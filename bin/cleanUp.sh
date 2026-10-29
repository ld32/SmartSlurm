#!/bin/bash

#set -x

# M1: sentinel marking the end of the step's own (tool) output and the start of
# the SmartSlurm job report / teardown. Emitted first so everything below --
# including the background "Cleanup counter" ticks -- is clearly SmartSlurm's,
# leaving the tool output above cleanly bracketed. Fires on success and failure
# because this script runs from the job's EXIT trap.
echo
echo "========== SmartSlurm job report: $1 =========="
echo

for i in {1..200}; do sleep 1; echo Cleanup counter: $i; done & 

# to call this:  0     1         2       3       4          5       6       7         8     9      10         11       12     13       14          15
#cleanUp.sh          "flag "software" "$ref" "$inputSize" $core   $memDef  $minDef   $mem  $time  $partition slurmAcc  inputs extraM extraTime userEmail

output="Running: $0"
for param in "$@"; do
    if [[ "$param" == *\ * ]]; then
        output="$output \"$param\""
    else
        output="$output $param"
    fi
done
echo "$output"

flag=$1 
software=$2
ref=$3
inputSize=$4 # input might not exist when job was submitted.
core=$5
memDef=$6
minDef=$7
totalM=$8

# # if jobs has --mem
# totalM=$SLURM_MEM_PER_NODE

# # if job has --mem-per-cpu and -c
# [ -z "$totalM" ] && totalM=$((SLURM_MEM_PER_CPU * SLURM_JOB_CPUS_PER_NODE))

# # if job has --mem-per-cpu and -n
# [ -z "$totalM" ] &&  totalM=$((SLURM_MEM_PER_CPU * SLURM_CPUS_PER_TASK))

totalT=${9}
partition=${10}
slurmAcc=${11}
inputs=${12}
extraMemC=${13}
extraTime=${14}
userEmail=${15}
excludeFailedNodes=${16}
pipeline=${17:-unknown}          # v2: the pipeline this step belongs to (ssbatch -L)

# v2 jobRecord: schema helpers + provenance. The record= column ORDER below must
# match JR_COLUMNS_V2 in jobRecordSchema.sh (name-addressed on read).
. "$(dirname "$0")/jobRecordSchema.sh" 2>/dev/null
jrNode=${SLURMD_NODENAME:-$(hostname 2>/dev/null)}
jrSample=$(jrSampleFromFlag "$flag")

# Defense in depth: make sure the store is v2 before we read/prune/append. ssbatch
# migrates at submit, but a v1/mixed file reaching cleanUp (older job, direct use)
# must be migrated here too, or the v1-positional prune below would delete the v2
# rows we append. Idempotent + locked, so concurrent jobs are safe.
[ -n "${smartSlurmJobRecordDir:-}" ] && jrMigrateStore "$smartSlurmJobRecordDir/jobRecord.txt"

# statKey names the per-group fitted .stat/.mem/.time/.memTime files, exactly as
# estimateResource.sh and calculateMemTime.sh write them: pipeline.step.refGroup
# (refGroup pools same-size refs via jrRefClass). When this job over-reserves or
# fails, those group files are stale and must be dropped so they regenerate from
# the updated jobRecord -- so the invalidation globs below MUST use this key, not
# the old software.ref. (extraMem.* stays keyed software.ref, matching its reader.)
refGroup=$(jrRefClass "$ref" 2>/dev/null); [ -n "$refGroup" ] || refGroup="$ref"
statKey="$pipeline.$software.$refGroup"

# jrFailPrune <step> <ref> <backupFile> : rewrite jobRecord.txt removing failed and
# failure-disproven records for one step/ref (keeps validated COMPLETED history).
# Version-branched: v1 backup uses the original positional two-pass unchanged; v2
# uses the same logic name-addressed (and preserves the header). Reads the backup
# twice (threshold pass, then filter pass), writing the pruned store.
jrFailPrune() {
    local a=$1 b=$2 backup=$3 out="$smartSlurmJobRecordDir/jobRecord.txt"
    # Pass 1 (read-only, name-addressed via the shared shim): the failure envelope
    # for this step/ref -- the smallest failing input size and largest failing
    # allocation, for OOM and OOT separately.
    local th mN fsM fmM tN fsT ftM
    th=$(jrEmitV2 "$backup" | awk -F, -v a="$a" -v b="$b" "$JR_AWK_COLMAP"'
            $(col["step"])==a && $(col["ref"])==b && $(col["status"])!="COMPLETED"{
              if($(col["status"])~/OOM/){ if(!mN++ || $(col["inputSize"])<fsM)fsM=$(col["inputSize"]); if($(col["memAllocated"])>fmM)fmM=$(col["memAllocated"]) }
              if($(col["status"])~/OOT/){ if(!tN++ || $(col["inputSize"])<fsT)fsT=$(col["inputSize"]); if($(col["timeAllocated"])>ftM)ftM=$(col["timeAllocated"]) }
            }
            END{ print mN+0, fsM+0, fmM+0, tN+0, fsT+0, ftM+0 }')
    read -r mN fsM fmM tN fsT ftM <<< "$th"
    # Pass 2 (atomic, locked, header-preserved): drop this step/ref's failed rows
    # and its "failure-disproven" COMPLETED rows (a success at >= a failing size
    # using <= the failing allocation -- an under-reporter); keep everything else.
    # jrRewrite reads from the live store; seed it from the backup first so the
    # pruned store is built from the same snapshot the thresholds came from.
    cp -p "$backup" "$out" 2>/dev/null
    jrRewrite "$out" '
        keep=1
        if ($(col["step"])==a && $(col["ref"])==b) {
            if ($(col["status"])!="COMPLETED") keep=0
            else {
                s=0
                if (mN && $(col["inputSize"])>=fsM && $(col["memUsed"])<=fmM) s=1
                if (tN && $(col["inputSize"])>=fsT && $(col["timeUsed"])<=ftM) s=1
                if (s) keep=0
            }
        }
        if (keep) print' \
        -v a="$a" -v b="$b" -v mN="$mN" -v fsM="$fsM" -v fmM="$fmM" -v tN="$tN" -v fsT="$fsT" -v ftM="$ftM"
}

#skipEstimate=${15} # todo. remove it
#smartSlurmJobRecordDir=${15}

# if not successful, delete cromwell error file, so that the job not labeled fail
#ls -l execution
[ -f $smartSlurmLogDir/$flag.success ] || rm -r execution/rc 2>/dev/null

out=$smartSlurmLogDir/"$flag.out"; err=$smartSlurmLogDir/$flag.err; script=$smartSlurmLogDir/$flag.sh; succFile=$smartSlurmLogDir/$flag.success; failFile=$smartSlurmLogDir/$flag.failed; checkpointDir=$smartSlurmLogDir/$flag

[ -f .exitcode ] && touch $succFile

# wait for slurm database update
sleep 5

sacct=`sacct --format=JobID,Submit,Start,End,MaxRSS,State,NodeList%30,Partition,ReqTRES%30,TotalCPU,Elapsed%14,Timelimit%14,ExitCode --units=M -j $SLURM_JOBID`

#sacct=`cat ~/fakeSacct.txt`

#sacct report is not accurate, let us use the total memory
#totalM=${sacct#*,mem=}; totalM=${totalM%%M,n*}

echo -e "\nJob summary:\n$sacct"
# echo *Notice the sacct report above: while the main job is still running for sacct command, user task is completed.

# record job for future estimating mem and time
jobStat=`echo -e "$sacct" | tail -n 1`

#from: "sacct --format=JobID,Submit,Start,End,MaxRSS,State,NodeList%30,Partition,ReqTRES%30,TotalCPU,Elapsed%14,Timelimit%14 --units=M -j $SLURM_JOBID"

START=`head -n 1 $smartSlurmLogDir/job_$SLURM_JOBID.memCPU.txt | cut -d' ' -f6`

FINISH=`date +%s`

[ -z "$START" ] && START=$((FINISH - 600))

echo  start: $START fisnish: $FINISH

# time in minutes
min=$((($FINISH - $START + 59)/60))

# memory in M
memSacct=`echo $jobStat | cut -d" " -f5`

[[ "$memSacct" == "RUNNING" ]] && memSacct=NA

node=`echo $jobStat | cut -d" " -f7`

if [ -f $succFile ] ; then # or nextflow successful job
    [ -z "$snakemakeSuccFlag" ] || touch "$snakemakeSuccFlag"
    jobStatus=COMPLETED
elif [[ "$sacct" == *TIMEOUT* ]] && [[ "$sacct" == *OUT_OF_ME* ]]; then
    jobStatus=OOMOOT    
elif [[ "$sacct" == *TIMEOUT* ]]; then
    jobStatus=OOT
elif [[ "$sacct" == *OUT_OF_ME* ]]; then
    jobStatus=OOM
elif [[ "$sacct" == *FAILED* ]]; then
    jobStatus=Fail
    [ -z "$snakemakeFailFlag" ] || touch "$snakemakeFailFlag"
elif [[ "$sacct" == *CANCEL* ]]; then
    jobStatus=Canceled
    [ -z "$snakemakeFailFlag" ] || touch "$snakemakeFailFlag"
else
    tLog=`tail -n 22 $out | grep ^srun`
    if [[ "$tLog" == *"task 0: Out Of Memory"* ]] || [[ "$tLog" == *"Cannot allocate memory"* ]]; then
        jobStatus="OOM"
        echo The job is actually out-of-memory by according to the log:
        echo $tLog
        scontrol show jobid -dd $SLURM_JOBID
    else
       jobStatus="Unknown"
       [ -z "$snakemakeFailFlag" ] || touch "$snakemakeFailFlag"
   fi
fi

echo jobStatus: $jobStatus

#echo -e  "Last row of job summary: $jobStat"
echo start: $START finish: $FINISH mem: $memSacct min: $min

# sacct actually give very not accurate result. Let's use cgroup report
#mem=`cat /sys/fs/cgroup/memory/slurm/uid_*/job_$SLURM_JOBID/memory.max_usage_in_bytes`

srunM=`cut -d' ' -f2 $smartSlurmLogDir/job_$SLURM_JOBID.memCPU.txt | sort -n | tail -n1`
#srunM=$((srunM / 1024 / 1024 ))

[ -z "$srunM" ] && srunM=0

echo jobStatus: $jobStatus cgroupMaxMem: $srunM

memSacct=${memSacct%M}; memSacct=${memSacct%.*} #remove M and decimals

# Not sure if this is needed.
[[ "$memSacct" != "NA" ]] && [ "$memSacct" -gt "$srunM" ] && srunM=$memSacct

if [[ "$inputs" != "none" ]] && [[ "$inputSize" == "0" ]]; then
    inputSize=`{ du --apparent-size -c -L ${inputs//,/ } 2>/dev/null || echo notExist; } | tail -n 1 | cut -f 1`

    if [[ "$inputSize" == "notExist" ]]; then
        echo Some or all input files not exist: .$inputs.
        echo Error! missingInputFile: .${inputs//,/ }.
        echo The input size is used for job records to estimame later job resournce needs.
        #exit
    fi
fi

# three possiblities to have this file:
# 1 job was OOT
# 2 job was oom
# 3 job resournce was adjusted by upsteam jobs
# 4 job checkpointed but OOT or OOM # this is special, because need to survive re-run from ssbatch
if [ -f ${out%.out}.adjust ]; then
   #tText=`cat ${out%.out}.adjust`
   #totalM=`echo $tText | cut -d' ' -f1`
   #totalT=`echo $tText | cut -d' ' -f2`
   #extraMC=`echo $tText | cut -d' ' -f3`
   IFS=' ' read -r inputSize totalM totalT extraMemC  <<< `cat $smartSlurmLogDir/$flag.adjust`
fi 

#[ -f $smartSlurmJobRecordDir/stats/extraMem.$software.$ref ] && extraMem=`sort -nr $smartSlurmJobRecordDir/stats/extraMem.$software.$ref | head -n1`


# totalM=200; totalT=200; srunM=1000; min=200;
overEstimate=""; overText=""; ratioM=""; ratioT=""
if [ "$memDef" -ne "$totalM" ] || [ "$minDef" -ne "$totalT" ]; then
  if [[ "$jobStatus" == COMPLETED ]]; then
    if [ "$totalM" -gt $((srunM * 2)) ] && [ $srunM -ge 1024 ]; then
        overEstimate="O:"
        overText="$SLURM_JOBID over-reserved resounce M: $srunM/$totalM T: $min/$totalT"
	#echo -e "`pwd`\n$out\n$SLURM_JOBID over-reserved resounce M: $srunM/$totalM T: $min/$totalT" | mail -s "Over reserved $SLURM_JOBID" $USER
        echo $overText

        # delete the .stat files, so that next time will re-generate
        mv $smartSlurmJobRecordDir/stats/$statKey.* $smartSlurmJobRecordDir/stats/back 2>/dev/null

        # get current extraMem
        currentExtraMem=0
        [ -f $smartSlurmJobRecordDir/stats/extraMem.$software.$ref ] && currentExtraMem=`sort -n $smartSlurmJobRecordDir/stats/extraMem.$software.$ref | tail -n1 | cut -d' ' -f1`

    fi
  fi
fi
ratioM=`echo "scale=2;$srunM/$totalM"|bc`; ratioT=`echo "scale=2;$min/$totalT"|bc`

# todo: move this part to main job, so that when release job, this job record can be used for the statics
# v2 column order == JR_COLUMNS_V2: pipeline,step,ref,sample,jobID,inputSize,memUsed,
#   memAllocated,timeUsed,timeAllocated,cores,status,user,node,epoch,isoDate,flag
record=$(jrRow pipeline="$pipeline" step="$2" ref="$ref" sample="$jrSample" jobID="$SLURM_JOBID" inputSize="$inputSize" memUsed="$srunM" memAllocated="$totalM" timeUsed="$min" timeAllocated="$totalT" cores="$core" status="$jobStatus" user="$USER" node="$jrNode" epoch="$(date +%s)" isoDate="$(date +%Y-%m-%dT%H:%M:%S)" flag="$flag")
# (dataToPlot mirrors the record; its only readers are the default-off ssMakePlots
#  plot parsers below, which feat/retire-unused-plots removes.)
echo dataToPlot,$record

# --- per-run efficiency row -------------------------------------------------
# One tab-separated row per job, written UNCONDITIONALLY (unlike the jobRecord
# append, which is gated), to a per-run file that checkRun reads to show
# time/mem efficiency. Units live in the header names so the file reads on its
# own: allocMemMB/usedMemMB (=$totalM/$srunM, megabytes) and
# allocTimeMin/usedTimeMin (=$totalT/$min, minutes). A human-readable ISO `date`
# trails the machine epoch (kept for sorting). The sample key is the first
# input's basename (falling back to the program name), matching how
# runAsPipeline groups the A+ per-sample rows. checkRunGrid reads by column
# position and skips this header row, so the header names are free to carry
# units and the trailing date column leaves every used index unchanged.
effFile="$smartSlurmLogDir/efficiency.tsv"
[ -f "$effFile" ] || printf 'epoch\tflag\tsample\tprogram\tref\tstatus\tallocMemMB\tusedMemMB\tallocTimeMin\tusedTimeMin\tinputs\tdate\n' > "$effFile"
effSample="${inputs#,}"; effSample="${effSample%%,*}"; effSample="${effSample%% *}"; effSample=`basename "$effSample" 2>/dev/null`
{ [ -z "$effSample" ] || [ "$effSample" = none ]; } && effSample="$2"
effEpoch=`date +%s`; effDate=`date '+%Y-%m-%dT%H:%M:%S'`
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$effEpoch" "$flag" "$effSample" "$2" "$ref" "$jobStatus" \
    "$totalM" "$srunM" "$totalT" "$min" "$inputs" "$effDate" >> "$effFile"

# 7-day prune: keep rows that are recent (epoch within a week) OR COMPLETED. One
# name-addressed jrRewrite call -- migrates first, locks, rewrites atomically,
# preserves the header. (Replaces the old version-branched inline awk, including
# the v1-positional fallback that deleted v2 rows read through the wrong columns.)
jrRewrite "$smartSlurmJobRecordDir/jobRecord.txt" \
    'if (($(col["epoch"])+0 > t) || ($(col["status"])=="COMPLETED")) print' \
    -v t="$(date -d '7 days ago' +%s)"


records=`jrEmitV2 "$smartSlurmJobRecordDir/jobRecord.txt" | awk -F, -v a="$2" -v b="$3" "$JR_AWK_COLMAP"'$(col["step"])==a && $(col["ref"])==b { print $(col["inputSize"]), $(col["memUsed"]) }' | sort -u -n`

#   less than 500 records or current one is larger than all old data and not exist already
if [[ "$jobStatus" == "COMPLETED" ]] && [[ "$(printf '%s\n' "$records" | grep -c .)" -lt 500 || "$(printf '%s\n' "$records" | tail -n1 | cut -d' ' -f1)" -lt "$inputSize" ]] && ! printf '%s\n' "$records" | grep -qxF "$inputSize $srunM"; then
    #if [ ! -f ${out%.out}.startFromCheckpoint ]; then
        [ -z "$START" ] || jrAppendLine "$smartSlurmJobRecordDir/jobRecord.txt" "$record"
        echo -e "Added this line to $smartSlurmJobRecordDir/jobRecord.txt:\n$record"
    #else
    #    echo -e "Has ${out%.out}.startFromCheckpoint. Did not add this line to $smartSlurmJobRecordDir/jobRecord.txt:\n$record"
    #fi
    mv $smartSlurmJobRecordDir/stats/$statKey.* $smartSlurmJobRecordDir/stats/back 2>/dev/null
else
    echo Did not add this record to $smartSlurmJobRecordDir/jobRecord.txt
fi


# If >10 records for this step/ref landed since the fitted .stat was built, drop
# the group .stat files so the next estimate regenerates them from fresh data.
# v2: compare the numeric `epoch` column against the stat file's mtime, name-
# addressed via jrEmitV2 (the old code read v1 positional fields $21/$12/$13 and
# split the human date, which no longer exist -- it silently counted 0).
if [ -f $smartSlurmJobRecordDir/stats/$statKey.mem.stat ]; then
    tFile=`stat -c %Y $smartSlurmJobRecordDir/stats/$statKey.mem.stat`

    newRecords=`jrEmitV2 "$smartSlurmJobRecordDir/jobRecord.txt" | awk -F, -v prog="$2" -v refx="$3" -v t="$tFile" "$JR_AWK_COLMAP"'($(col["epoch"])+0) > t && $(col["step"])==prog && $(col["ref"])==refx' | wc -l`
    if [ $newRecords -gt 10 ]; then
        mv $smartSlurmJobRecordDir/stats/$statKey.*.stat $smartSlurmJobRecordDir/stats/back 2>/dev/null
    fi
fi

echo Final mem: $srunM M, time: $min mins

# for nextflow log
cat .command.log >> $out 2>/dev/null 

# for testing
#rm $succFile; jobStatus=OOM
#rm $succFile
#jobStatus=OOM
if [ ! -f $succFile ]; then

    if  [ -f "${out%.out}.likelyCheckpointDoNotWork" ] && [ $((srunM/totalM*100)) -lt 10 ]; then
        echo -e "This step might not be good to checkpoint?\nOr maybe you did not give enough memory? Please test run it again with higher memory:\n$flag" | mail -s "!!!!$SLURM_JOBID:$flag" $USER
    fi
    touch $failFile
    # if checkpoint due to low memory or checkpoint failed due to memory, or actually out of memory
    [ -f "${out%.out}.likelyCheckpointOOM" ] && jobStatus=OOM 


    # Handle OOM, OOT, and OOMOOT cases
    if [[ "$jobStatus" == "OOM" ]] || [[ "$jobStatus" == "OOT" ]] || [[ "$jobStatus" == "OOMOOT" ]]; then

        record=$(jrRow pipeline="$pipeline" step="$2" ref="$ref" sample="$jrSample" jobID="$SLURM_JOBID" inputSize="$inputSize" memUsed="$srunM" memAllocated="$totalM" timeUsed="$min" timeAllocated="$totalT" cores="$core" status="$jobStatus" user="$USER" node="$jrNode" epoch="$(date +%s)" isoDate="$(date +%Y-%m-%dT%H:%M:%S)" flag="$flag")

        jrAppendLine "$smartSlurmJobRecordDir/jobRecord.txt" "$record"

        recentFails=`jrEmitV2 "$smartSlurmJobRecordDir/jobRecord.txt" | awk -F, -v a="$2" -v b="$3" -v t="$(date -d '1 hour ago' +%s)" "$JR_AWK_COLMAP"'$(col["status"])!="COMPLETED" && ($(col["epoch"])+0)>t && $(col["step"])==a && $(col["ref"])==b' | wc -l`
        echo recentFails in the last hour for $software $ref: $recentFails



        # delete all record about the step is thare more than 3 failed jobs in last hour
        if [ $recentFails -gt $numberOfFailedJobsLastHourToRecalculateStats ]; then
            backup=$smartSlurmJobRecordDir/jobRecord.back.$(date +%Y-%m-%d_%H-%M-%S).txt 
            if [ ! -f $backup ]; then 
                mv $smartSlurmJobRecordDir/jobRecord.txt $backup 2>/dev/null
                echo "Too many fails in the last hour for $software $ref, backing jobRecord.txt to# $backup and removing failed + failure-disproven records for $software $ref (validated history kept)." 2>/dev/null   
                jrFailPrune "$2" "$3" "$backup" 
                #ls -l $smartSlurmJobRecordDir/ 
                #ls -l $smartSlurmJobRecordDir/stats/
                rm $smartSlurmJobRecordDir/stats/$statKey.* 2>/dev/null
            fi
            # email user to review jobRecord.txt and related stat files, because there are too many fails in the last hour for this software and ref, which might indicate some systemic issue with the job or the cluster
            echo -e "There are $recentFails fails in the last hour for $software $ref, which might indicate some systemic issue with the job or the cluster. All job records for this software and ref  have been backed up to $backup; failed and failure-disproven records were removed, validated COMPLETED history kept in jobRecord.txt." | mail -s "!!!Too many oom/oot/over-estimates for $software $ref" $USER

        fi


        echo "Will re-queue after sending email... Job Status: $jobStatus"

        # Handle memory adjustment for OOM or OOMOOT
        if [[ "$jobStatus" == "OOM" ]] || [[ "$jobStatus" == "OOMOOT" ]]; then
            echo "Processing memory adjustment for $jobStatus"
            echo "old extraMem:"
            cat $smartSlurmJobRecordDir/stats/extraMem.$software.$ref

            extraMemN=$(( totalM - srunM + 1 ))
            maxExtra=`sort -n $smartSlurmJobRecordDir/stats/extraMem.$software.$ref | tail -n1 | cut -d' ' -f1`
            [ -z "$maxExtra" ] && maxExtra=5
            [ $extraMemN -gt $maxExtra ] && echo $extraMemN $totalM $inputSize $SLURM_JOBID >> $smartSlurmJobRecordDir/stats/extraMem.$software.$ref && maxExtra=$extraMemN

            echo "new extraMem:"
            cat $smartSlurmJobRecordDir/stats/extraMem.$software.$ref

            mv ${out%.out}.likelyCheckpointOOM ${out%.out}.likelyCheckpointOOM.old 2>/dev/null

            # Calculate new memory via the shared escalation helper (config.txt) so a
            # later manual rerun's failure floor escalates identically. (+2*maxExtra
            # margin preserved.)
            mem=`nextMemAfterFail $totalM $recentFails $((maxExtra*2))`
            echo "Calculated new memory: $mem M"
        else
            # For OOT only, keep existing memory
            mem=$totalM
        fi

        # Handle time adjustment for OOT or OOMOOT
        if [[ "$jobStatus" == "OOT" ]] || [[ "$jobStatus" == "OOMOOT" ]]; then
            echo "Processing time adjustment for $jobStatus"
            # Calculate new time via the shared escalation helper (config.txt).
            min=`nextTimeAfterFail $min $recentFails`
            echo "Calculated new time: $min minutes"
        fi

        # Calculate partition and time format
        #hours=$((($min + 59) / 60))
        #partition=`scontrol show job $id | grep Partition= | awk '{for(i=1;i<=NF;i++) if ($i ~ /Partition=/) {split($i,a,"="); print a[2];}}'`

        partition=`grep "^#myMem=" $script | awk -F'myPartition=' '{print $2}' | awk '{print $1}'`


        adjustPartition $mem $min $partition
        
        seconds=$(($min * 60))
        time=`eval "echo $(date -ud "@$seconds" +'$((%s/3600/24))-%H:%M:%S')"`

        echo $inputSize $mem $min $maxExtra > ${out%.out}.adjust

        # Set up environment variables for requeue command
        export myPartition=$partition
        export myTime=$time
        export myMem=${mem}M

        requeueCmd=`grep "Command used to submit the job:" $script | tail -n 1`
        requeueCmd=${requeueCmd#*submit the job: }
        requeueCmd=${requeueCmd/ -H/}
        requeueCmd=${requeueCmd//\$myPartition/$myPartition}
        requeueCmd=${requeueCmd//\$myTime/$myTime}
        requeueCmd=${requeueCmd//\$myMem/$myMem}
        requeueCmd=$( echo $requeueCmd | sed -E 's/-d afterok:[0-9]+//g' )

        # Attempt to requeue with retry logic
        for try in {1..3}; do
            # Create different flag files for different failure types
            #flagFile="$failFile.requeued.$try"
            case "$jobStatus" in
                "OOM") flagFile="$failFile.requeued.$try.mem" ;;
                "OOT") flagFile="$failFile.requeued.$try.time" ;;
                "OOMOOT") flagFile="$failFile.requeued.$try.oomoot" ;;
            esac

            if [ ! -f "$flagFile" ]; then
                touch "$flagFile"
                
                echo "Attempt $try - cmd: $requeueCmd"
                newJobID=`$requeueCmd`

                if [[ "$newJobID" =~ ^[0-9]+$ ]]; then
                    echo "# mem=$myMem time=$myTime reason=$jobStatus" >> $script
                    
                    # Update dependent jobs
                    IFS=$'\n'
                    for line in `grep $SLURM_JOBID $smartSlurmLogDir/allJobs.txt | grep -v ^$SLURM_JOBID`; do
                        job=${line%% *}
                        deps=`echo $line | awk '{print $2}'`

                        if [[ $deps == null ]]; then
                            deps=""
                        elif [[ $deps == ${deps/:/} ]]; then
                            deps="Dependency=afterok:$newJobID"
                        else
                            tmp=""; IFS=$' '
                            for t in ${deps//:/ }; do
                                [ "$SLURM_JOBID" == "$t" ] && tmp=$tmp:$newJobID || tmp=$tmp:$t
                            done
                            [ -z "$tmp" ] && deps="" || deps="Dependency=afterok$tmp"
                        fi
                        
                        if [ ! -z "$deps" ]; then
                            scontrol update jobid=$job $deps
                            flg=`echo $line | awk '{print $3}'`
                            echo `grep "Command used to submit the job:" $smartSlurmLogDir/$flg.sh | tail -n 1 | sed "s/$SLURM_JOBID/$newJobID/" ` >> $smartSlurmLogDir/$flg.sh
                        fi
                    done

                    echo "Requeued successfully"
                    [ -f $failFile ] && rm $failFile

                    # Send notification email with specific failure type
                    case "$jobStatus" in
                        "OOM") s="OOM.Requeued:$SLURM_JOBID-$newJobID:$SLURM_JOB_NAME" ;;
                        "OOT") s="OOT.Requeued:$SLURM_JOBID-$newJobID:$SLURM_JOB_NAME" ;;
                        "OOMOOT") s="OOMOOT.Requeued:$SLURM_JOBID-$newJobID:$SLURM_JOB_NAME" ;;
                    esac
                    
                    # Auto-requeue is silent now: it is in-progress recovery, not a
                    # break or a completion. Outcomes are reported by the per-sample
                    # status email, the per-run efficiency email, and checkRun.

                    # (no longer touch $out.$newJobID: 0-byte job-id marker, only
                    #  read cosmetically by checkRun; requeue is recorded below.)
                    cp $smartSlurmLogDir/allJobs.txt $smartSlurmLogDir/allJobs.requeue.$SLURM_JOBID.as.$newJobID

                    # Lock mechanism for job ID adjustments
                    while true; do
                        if mkdir $smartSlurmLogDir/job.adjusting.lock  2>/dev/null; then
                            break
                        else 
                            folder_mtime=$(stat -c %Y $smartSlurmLogDir/job.adjusting.lock )
                            current_time=$(date +%s)
                            if [[ $((current_time - folder_mtime)) -gt 10 ]]; then
                                touch $smartSlurmLogDir/job.adjusting.lock 
                                break
                            fi    
                        fi
                        echo "waiting for lock to adjust job ids"
                        sleep 1
                    done  

                    sed -i "s/$SLURM_JOBID/$newJobID/" $smartSlurmLogDir/allJobs.txt
                    rm -r $smartSlurmLogDir/job.adjusting.lock

                    cp $smartSlurmLogDir/job_$SLURM_JOBID.memCPU.txt $smartSlurmLogDir/job_$newJobID.memCPU.txt
                    echo 0 0 0 0 0 0 0 >> $smartSlurmLogDir/job_$newJobID.memCPU.txt
                    break

                else
                    echo "re-submit failed on attempt $try"
                fi
            fi
        done

        # Clean up statistics files
        if [ "$inputs" == "none" ]; then
            mv $smartSlurmJobRecordDir/stats/$statKey.* $smartSlurmJobRecordDir/stats/back  2>/dev/null
        else
            if [ -f $smartSlurmJobRecordDir/stats/$statKey.mem.stat ] || [ -f $smartSlurmJobRecordDir/stats/$statKey.time.stat ]; then
                mv $smartSlurmJobRecordDir/stats/$statKey.* $smartSlurmJobRecordDir/stats/back  2>/dev/null
            fi
        fi

    else
        echo Not sure why this job failed. Did not run out of time or memory. Please check for yourself.
        
        if [ ! -z "$excludeFailedNodes" ]; then 
            echo -n ,$SLURMD_NODENAME >> $smartSlurmJobRecordDir/stats/badNodes.$software.$ref
            
            badNodeList="-x `cat $smartSlurmJobRecordDir/stats/badNodes.$software.$ref`"
   
            export myPartition=$partition
            export myTime=$totalT
            export myMem=${totalM}M
            requeueCmd=`grep "Command used to submit the job:" $script | tail -n 1`
            requeueCmd=${requeueCmd#*submit the job: }
            requeueCmd=${requeueCmd//\$myPartition/$myPartition}
            requeueCmd=${requeueCmd//\$myTime/$myTime}
            requeueCmd=${requeueCmd//\$myMem/$myMem $badNodeList}
            requeueCmd=${requeueCmd/ -H/}
            requeueCmd=$( echo $requeueCmd | sed -E 's/-d afterok:[0-9]+//g' ) #-d afterok:38023271

            for try in {1..3}; do
                if [ ! -f $failFile.requeued.$try.fail ]; then
                    touch $failFile.requeued.$try.fail
                    
                    echo cmd: $requeueCmd

                    newJobID=`$requeueCmd`

                    if [[ "$newJobID" =~ ^[0-9]+$ ]]; then
                        echo "# mem=$myMem time=$myTime " >> $script
                        IFS=$'\n'
                        for line in `grep $SLURM_JOBID $smartSlurmLogDir/allJobs.txt | grep -v ^$SLURM_JOBID`; do
                            job=${line%% *}
                            deps=`echo $line | awk '{print $2}'`

                            if [[ $deps == null ]]; then
                                deps=""
                            elif [[ $deps == ${deps/:/} ]]; then
                                deps="Dependency=afterok:$newJobID"
                            else
                                tmp=""; IFS=$' '
                                for t in ${deps//:/ }; do
                                    [ "$SLURM_JOBID" == "$t" ] && tmp=$tmp:$newJobID || tmp=$tmp:$t
                                done
                                [ -z "$tmp" ] && deps="" || deps="Dependency=afterok$tmp"
                            fi
                            if [ ! -z "$deps" ]; then
                                scontrol update jobid=$job $deps
                                flg=`echo $line | awk '{print $3}'`
                                echo `grep "Command used to submit the job:" $smartSlurmLogDir/$flg.sh | tail -n 1 | sed "s/$SLURM_JOBID/$newJobID/" ` >> $smartSlurmLogDir/$flg.sh
                            fi
                        done

                        #if `sh $PWD/$smartSlurmLogDir/$flag.requeueCMD; rm $PWD/$smartSlurmLogDir/$flag.requeueCMD;`; then
                        #rm $smartSlurmLogDir/$flag.requeueCMD #;  then
                        #if `srun --jobid $SLURM_JOBID $acc "pwd"`; then
                        #if `scontrol requeue $SLURM_JOBID; sleep 2; scontrol update JobId=$SLURM_JOBID MinMemoryNode=$mem`; then

                        echo Re-queued successfully
                        [ -f $failFile ] && rm $failFile

                        # Auto-requeue (bad node) is silent now; see note above.

                        # (no longer touch $out.$newJobID: see note at the other
                        #  requeue site -- 0-byte cosmetic job-id marker.)
                        cp $smartSlurmLogDir/allJobs.txt $smartSlurmLogDir/allJobs.requeue.$SLURM_JOBID.as.$newJobID
                        while true; do
                            if mkdir $smartSlurmLogDir/job.adjusting.lock  2>/dev/null; then
                                break
                            else 
                                folder_mtime=$(stat -c %Y $smartSlurmLogDir/job.adjusting.lock )
                                current_time=$(date +%s)

                                if [[ $((current_time - folder_mtime)) -gt 10 ]]; then
                                    touch $smartSlurmLogDir/job.adjusting.lock 
                                    break
                                fi    
                            fi
                            echo waiting for lock to adjust job ids 
                            sleep 1
                        done                     

                        sed -i "s/$SLURM_JOBID/$newJobID/" $smartSlurmLogDir/allJobs.txt
                        rm -r $smartSlurmLogDir/job.adjusting.lock
                        cp $smartSlurmLogDir/job_$SLURM_JOBID.memCPU.txt $smartSlurmLogDir/job_$newJobID.memCPU.txt
                        echo 0 0 0 0 0 0 0 >> $smartSlurmLogDir/job_$newJobID.memCPU.txt
                        break
                    
                    else 
                        echo re-submit failed;
                    fi    
                fi
            done
        fi 
    fi
    #if [ "$min" -ge 20 ] && [ ! -z "$alwaysRequeueIfFail" ] && [ "$jobStatus" != "Cancelled" ]; then
    #    ( sleep 2;  scontrol requeue $SLURM_JOBID; ) &
    #fi

# job was successful
elif [ ! -z "$overEstimate" ]; then
        record=$(jrRow pipeline="$pipeline" step="$2" ref="$ref" sample="$jrSample" jobID="$SLURM_JOBID" inputSize="$inputSize" memUsed="$srunM" memAllocated="$totalM" timeUsed="$min" timeAllocated="$totalT" cores="$core" status="OverEstimate" user="$USER" node="$jrNode" epoch="$(date +%s)" isoDate="$(date +%Y-%m-%dT%H:%M:%S)" flag="$flag")
        jrAppendLine "$smartSlurmJobRecordDir/jobRecord.oom.oot.over.txt" "$record"

    
        recentFails=`jrEmitV2 "$smartSlurmJobRecordDir/jobRecord.txt" | awk -F, -v a="$2" -v b="$3" -v t="$(date -d '1 hour ago' +%s)" "$JR_AWK_COLMAP"'$(col["status"])!="COMPLETED" && ($(col["epoch"])+0)>t && $(col["step"])==a && $(col["ref"])==b' | wc -l`
        echo recentFails in the last hour for $software $ref: $recentFails

        # delete all record about the step is thare more than 3 failed jobs in last hour
        if [ $recentFails -gt $numberOfFailedJobsLastHourToRecalculateStats ]; then
            backup=$smartSlurmJobRecordDir/jobRecord.back.$(date +%Y-%m-%d_%H-%M-%S).txt 
            if [ ! -f $backup ]; then 
                mv $smartSlurmJobRecordDir/jobRecord.txt $backup 2>/dev/null
                echo Too many fails in the last hour for $software $ref, backing jobRecord.txt to# $backup and removing failed + failure-disproven records for $software $ref, validated history kept. 2>/dev/null   
                jrFailPrune "$2" "$3" "$backup" 
                #ls -l $smartSlurmJobRecordDir/ 
                #ls -l $smartSlurmJobRecordDir/stats/
                rm $smartSlurmJobRecordDir/stats/$statKey.* 2>/dev/null
            fi
            # email user to review jobRecord.txt and related stat files, because there are too many fails in the last hour for this software and ref, which might indicate some systemic issue with the job or the cluster
            echo -e "There are $recentFails fails in the last hour for $software $ref, which might indicate some systemic issue with the job or the cluster. All job records for this software and ref  have been backed up to $backup; failed and failure-disproven records were removed, validated COMPLETED history kept in jobRecord.txt." | mail -s "!!!Too many oom/oot/over-estimates for $software $ref" $USER

        fi

fi

# NOTE: the per-job/barchart gnuplot pngs, their .dataMem/.dataTime/.output CSVs,
# and the on-demand DAG png render used to live here (gated behind ssMakePlots).
# They are gone: efficiency now lives in efficiency.tsv + checkRun (press e /
# checkRunGrid --efficiency), and checkRun's g key renders the workflow graph
# live, so nothing in the pipeline needs to write those files. The per-run dollar
# savings figure, which was derived from the now-removed plot CSV, is dropped with
# it; the per-job savings below (from memCPU.txt) is unchanged.

rate=0.0013 # $0.0013 per G per hour
maxSaved=`cut -d' ' -f4 $smartSlurmLogDir/job_$SLURM_JOBID.memCPU.txt | sort -n | tail -n1`
jMin=`wc -l $smartSlurmLogDir/job_$SLURM_JOBID.memCPU.txt | cut -d ' ' -f 1`
savedDollar=$(echo "$rate * $maxSaved / 1024 * $jMin / 60" | bc -l)

maxUsed=`cut -d' ' -f2 $smartSlurmLogDir/job_$SLURM_JOBID.memCPU.txt | sort -n | tail -n1`
#jMin=`wc -l $smartSlurmLogDir/job_$SLURM_JOBID.memCPU.txt | cut -d ' ' -f 1`
savedDollar1=$(echo "$rate * ($defaultMem - $maxUsed) / 1024 * $jMin / 60" | bc -l)



minimumsize=9000
actualsize=`wc -c < "$out" || echo 0`

[ -f $succFile ] && s="${overEstimate}Succ:$SLURM_JOBID:$SLURM_JOB_NAME" || s="$jobStatus:$SLURM_JOBID:$SLURM_JOB_NAME"

if [ "$actualsize" -ge "$minimumsize" ]; then
   #toSend=`echo Job script content:; cat $script;`
   toSend="$overText\nLog path: $out"
   toSend="$toSend\nOutput is too big for email. Please find output from log mentioned above."
   toSend="$toSend\n...\nFirst 20 row of output:\n`head -n 20 $out`"
   toSend="$toSend\n...\nLast 20 row of output:\n`tail -n 20 $out`"
else
   #toSend=`echo Job script content:; cat $script; echo; echo Job output:; cat $out;`
   toSend="$overText\nLog path: $out"
   toSend="$toSend\nJob log:\n`cat $out`"
   #toSend="$s\n$toSend"
fi

toSend="SmartSlurm Saved $savedDollar by dynamic allocation of memory for this job.\n$toSend"
[ -z "$snakemakeSuccFlag" ] && toSend="SmartSlurm Saved $savedDollar1 by split workflow into steps from this job.\n$toSend"


if [ -f "$err" ]; then
    actualsize=`wc -c < "$err"`
    if [ "$actualsize" -ge "$minimumsize" ]; then
        toSend="$toSend\nError file is too big for email. Please find output in: $err"
        toSend="$toSend\n...\nLast 10 rows of error file:\n`tail -n 10 $err`"
    else
        toSend="$toSend\n\nError output:\n`cat  $err`"
    fi
fi

#cp /tmp/job_$SLURM_JOBID.mem.txt $smartSlurmLogDir/

# move to inside of job
#summarizeRun.sh $smartSlurmLogDir $flag 

[ -f $smartSlurmLogDir/summary.$SLURM_JOBID ] && toSend="`cat $smartSlurmLogDir/summary.$SLURM_JOBID`\n$toSend" && s="${toSend%% *} $s"

#echo -e "tosend:\n$toSend"
#echo -e "$toSend" >> ${err%.err}.email



[ -z "$userEmail" ] || USER=$userEmail


# ---------------------------------------------------------------------------
# Consolidated, text-only notifications (no png attachments). The old code sent
# one email per job with up to 7 pngs (and a copy to a hard-coded dev address);
# this replaces that with two separate, de-duplicated goals:
#   status      -> one per-SAMPLE alert: WARN on break, CONFIRM on completion
#   efficiency  -> one per-RUN email carrying the per-sample efficiency table,
#                  pointing at `checkRun` (press e) for the per-step breakdown
# Each fires once, guarded by a marker file, so volume is ~one per sample plus
# one per run instead of one per job. Honors --email noEmail (everything off)
# and noSuccEmail (failures only). The sample key is the first input's basename
# (fallback: program name), matching efficiency.tsv and allJobs.txt column 6.
# ---------------------------------------------------------------------------
_ss_sendmail() {   # <subject> <body>
    [ "$lessEmail" = noEmail ] && { echo "noEmail: suppressed [$1]"; return 0; }
    echo -e "$2" | mail -s "$1" "$USER" 2>/dev/null && echo "email sent: $1" || echo "email not sent: $1"
}

alljobs="$smartSlurmLogDir/allJobs.txt"

if [ ! -f "$succFile" ]; then
    # ---- sample BREAK: one warning per sample (first failure wins) ----
    brkMark="$smartSlurmLogDir/.emailed.break.$effSample"
    if [ ! -f "$brkMark" ]; then
        : > "$brkMark"
        _ss_sendmail "SmartSlurm FAIL: sample $effSample ($jobStatus)" \
            "Sample '$effSample' FAILED at step $flag ($jobStatus).\nRe-run the same command to continue -- SmartSlurm will resize and resubmit.\n\n$toSend"
    else
        echo "break already emailed for sample $effSample"
    fi
else
    # ---- sample COMPLETE: fire from the terminal (highest-step) job of the
    #      sample, so one confirmation per sample. Sample + step are read from
    #      allJobs.txt (col 3 = flag -> step; col 6 = inputs -> sample) the same
    #      way effSample is derived, so both sides agree.
    if [ "$lessEmail" != noSuccEmail ] && [ -f "$alljobs" ]; then
        myStep="${flag%%.*}"
        maxStep=`awk -v s="$effSample" '
            function samp(inp,   x,p,n){ sub(/^,/,"",inp); sub(/ .*/,"",inp); split(inp,x,","); n=split(x[1],p,"/"); return p[n] }
            NR>1 && NF>=6 { st=$3; sub(/\..*/,"",st); smp=samp($6); if(smp==""||smp=="none") smp=$4; if(smp==s && st+0>m) m=st+0 }
            END { print m+0 }' "$alljobs"`
        cmpMark="$smartSlurmLogDir/.emailed.done.$effSample"
        if [ "${myStep:-0}" -ge "${maxStep:-0}" ] && [ ! -f "$cmpMark" ]; then
            : > "$cmpMark"
            _ss_sendmail "SmartSlurm done: sample $effSample completed" \
                "Sample '$effSample' completed successfully (through step $myStep).\n\nSee 'checkRun' (press e) for resource efficiency."
        fi
    fi

    # ---- per-RUN efficiency email: fire once, when EVERY job has succeeded ----
    # (A run with any failed sample does not get a "complete" email; its break
    # warning covers it. Only a fully-successful run reports efficiency.)
    if [ "$lessEmail" != noSuccEmail ] && [ "$lessEmail" != noEmail ] && [ -f "$alljobs" ]; then
        total=`awk 'NR>1 && NF>=3' "$alljobs" | wc -l`
        succN=0
        while read -r fl; do
            [ -z "$fl" ] && continue
            [ -e "$smartSlurmLogDir/$fl.success" ] && succN=$((succN+1))
        done < <(awk 'NR>1 && NF>=3{print $3}' "$alljobs")
        # Fire only from a DOWNSTREAM/terminal step -- one that depends on an
        # earlier step (the flag's 2nd field, the dependency, is non-zero). A
        # source producer (dep field 0) is excluded, so it cannot prematurely
        # trigger the run email when allJobs.txt is still being populated (as in
        # the synchronous harness, where submit == run). On real Slurm allJobs is
        # complete up front, so succN>=total already pins the true last job; this
        # guard just makes the harness agree.
        myDep=`echo "$flag" | cut -d. -f2`
        runMark="$smartSlurmLogDir/.emailed.effrun"
        if [ "${total:-0}" -ge 1 ] && [ "${succN:-0}" -ge "${total:-0}" ] && [ "${myDep:-0}" != 0 ] && [ ! -f "$runMark" ]; then
            : > "$runMark"
            effTable=`"$(dirname "$0")/checkRunGrid" --efficiency --no-color --log "$smartSlurmLogDir" 2>/dev/null`
            _ss_sendmail "SmartSlurm run complete: resource efficiency" \
                "All jobs in this run finished successfully.\n\n$effTable\n\nSee 'checkRun' (press e) for the per-step breakdown."
        fi
    fi
fi
sleep 5
