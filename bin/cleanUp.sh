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
        mv $smartSlurmJobRecordDir/stats/$software.$ref.* $smartSlurmJobRecordDir/stats/back 2>/dev/null

        # get current extraMem
        currentExtraMem=0
        [ -f $smartSlurmJobRecordDir/stats/extraMem.$software.$ref ] && currentExtraMem=`sort -n $smartSlurmJobRecordDir/stats/extraMem.$software.$ref | tail -n1 | cut -d' ' -f1`

    fi
  fi
fi
ratioM=`echo "scale=2;$srunM/$totalM"|bc`; ratioT=`echo "scale=2;$min/$totalT"|bc`

# todo: move this part to main job, so that when release job, this job record can be used for the statics 
                                #3defult,  5given,  7cGroupUsed                  sacct used
record="$SLURM_JOBID,$inputSize,$memDef,$minDef,$totalM,$totalT,$srunM,$min,$jobStatus,$USER,$memSacct,$2,$ref,$flag,$core,$extraMemC,$defaultExtraTime,0$ratioM,0$ratioT,`date`,`date +%s`"
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

# delete rows other than 1 week and not COMPLETED
#awk -F',' -v t="$(date -d '7 days ago' +%s)" '{d=$21; gsub(/^"|"$/,"",d); cmd="date -d \"" d "\" +%s"; cmd|getline ts; close(cmd); if(ts>t && $9!="COMPLETED") print}' "$smartSlurmJobRecordDir/jobRecord.txt" > "$smartSlurmJobRecordDir/jobRecord.txt.new" && mv "$smartSlurmJobRecordDir/jobRecord.txt.new" "$smartSlurmJobRecordDir/jobRecord.txt"  
awk -F"," -v t="$(date -d '7 days ago' +%s)" '{ if($21 > t || $9 == "COMPLETED") print $0 }' $smartSlurmJobRecordDir/jobRecord.txt > $smartSlurmJobRecordDir/jobRecord.txt.new && mv $smartSlurmJobRecordDir/jobRecord.txt.new $smartSlurmJobRecordDir/jobRecord.txt


records=`awk -F"," -v a=$2 -v b=$3 '{ if($12 == a && $13 == b) {print $2, $7 }}' $smartSlurmJobRecordDir/jobRecord.txt | sort -u -n`

#   less than 500 records or current one is larger than all old data and not exist already
if [[ "$jobStatus" == "COMPLETED" ]] && [[ "$(echo $records | wc -l)" -lt 500 || "`echo $records | tail -n1 | cut -d' ' -f1`" -lt "$inputSize" ]] && ! echo $records | grep "$inputSize $srunM"; then
    #if [ ! -f ${out%.out}.startFromCheckpoint ]; then
        [ -z "$START" ] || echo $record >> $smartSlurmJobRecordDir/jobRecord.txt
        echo -e "Added this line to $smartSlurmJobRecordDir/jobRecord.txt:\n$record"
    #else
    #    echo -e "Has ${out%.out}.startFromCheckpoint. Did not add this line to $smartSlurmJobRecordDir/jobRecord.txt:\n$record"
    #fi
    mv $smartSlurmJobRecordDir/stats/$software.$ref.* $smartSlurmJobRecordDir/stats/back 2>/dev/null
else
    echo Did not add this record to $smartSlurmJobRecordDir/jobRecord.txt
fi


# the last column is date. first get the mofification time for $smartSlurmJobRecordDir/stats/$software.$ref.mem.stat file and see how many records in jobrecord.txt have rows with time new than that 
# if more than 10, then delete the stat file ( will re-generate next time to estimate)
if [ -f $smartSlurmJobRecordDir/stats/$software.$ref.mem.stat ]; then
    tFile=`stat -c %Y $smartSlurmJobRecordDir/stats/$software.$ref.mem.stat`

    # todo: time format in jobRecord.txt is not consistent. some are in seconds, some are in date format. need to fix it.  
    newRecords=`awk -F"," -v prog=$2 -v refx=$3 -v t="$tFile" '{ split($21, d, " "); if( d[2] > t && $12 == prog && $13 == refx) print $0}' $smartSlurmJobRecordDir/jobRecord.txt | wc -l`
    if [ $newRecords -gt 10 ]; then
        mv $smartSlurmJobRecordDir/stats/$software.$ref.*.stat $smartSlurmJobRecordDir/stats/back 2>/dev/null
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

        record="$SLURM_JOBID,$inputSize,$memDef,$minDef,$totalM,$totalT,$srunM,$min,$jobStatus,$USER,$memSacct,$2,$ref,$flag,$core,$extraMemC,$defaultExtraTime,0$ratioM,0$ratioT,`date`,`date +%s`" 

        echo $record >> $smartSlurmJobRecordDir/jobRecord.txt

        recentFails=`awk -F"," -v a=$2 -v b=$3 -v t="$(date -d '1 hour ago' +%s)" '{ if(  $9 != "COMPLETED" && $21 > t && $12 == a && $13 == b ) print $0}' $smartSlurmJobRecordDir/jobRecord.txt | wc -l`
        echo recentFails in the last hour for $software $ref: $recentFails



        # delete all record about the step is thare more than 3 failed jobs in last hour
        if [ $recentFails -gt $numberOfFailedJobsLastHourToRecalculateStats ]; then
            backup=$smartSlurmJobRecordDir/jobRecord.back.$(date +%Y-%m-%d_%H-%M-%S).txt 
            if [ ! -f $backup ]; then 
                mv $smartSlurmJobRecordDir/jobRecord.txt $backup 2>/dev/null
                echo "Too many fails in the last hour for $software $ref, backing jobRecord.txt to# $backup and removing failed + failure-disproven records for $software $ref (validated history kept)." 2>/dev/null   
                awk -F"," -v a=$2 -v b=$3 'NR==FNR{ if($12==a && $13==b && $9!="COMPLETED"){ if($9 ~ /OOM/){ if(!mN++ || $2<fsM)fsM=$2; if($5>fmM)fmM=$5 } if($9 ~ /OOT/){ if(!tN++ || $2<fsT)fsT=$2; if($6>ftM)ftM=$6 } } next } { if(!($12==a && $13==b)){ print; next } if($9!="COMPLETED") next; s=0; if(mN && $2>=fsM && $7<=fmM) s=1; if(tN && $2>=fsT && $8<=ftM) s=1; if(!s) print }' $backup $backup > $smartSlurmJobRecordDir/jobRecord.txt 
                #ls -l $smartSlurmJobRecordDir/ 
                #ls -l $smartSlurmJobRecordDir/stats/
                rm $smartSlurmJobRecordDir/stats/$software.$ref.* 2>/dev/null
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

                    touch $out.$newJobID
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
            mv $smartSlurmJobRecordDir/stats/$software.$ref.* $smartSlurmJobRecordDir/stats/back  2>/dev/null
        else
            if [ -f $smartSlurmJobRecordDir/stats/$software.$ref.mem.stat ] || [ -f $smartSlurmJobRecordDir/stats/$software.$ref.time.stat ]; then
                mv $smartSlurmJobRecordDir/stats/$software.$ref.* $smartSlurmJobRecordDir/stats/back  2>/dev/null
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

                        touch $out.$newJobID

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
        record="$SLURM_JOBID,$inputSize,$memDef,$minDef,$totalM,$totalT,$srunM,$min,OverEstimate,$USER,$memSacct,$2,$ref,$flag,$core,$extraMemC,$defaultExtraTime,0$ratioM,0$ratioT,`date`,`date +%s`" 
        echo $record >> $smartSlurmJobRecordDir/jobRecord.oom.oot.over.txt

    
        recentFails=`awk -F"," -v a=$2 -v b=$3 -v t="$(date -d '1 hour ago' +%s)" '{ if( $9 != "COMPLETED" && $21 > t && $12 == a && $13 == b ) print $0}' $smartSlurmJobRecordDir/jobRecord.txt | wc -l`
        echo recentFails in the last hour for $software $ref: $recentFails

        # delete all record about the step is thare more than 3 failed jobs in last hour
        if [ $recentFails -gt $numberOfFailedJobsLastHourToRecalculateStats ]; then
            backup=$smartSlurmJobRecordDir/jobRecord.back.$(date +%Y-%m-%d_%H-%M-%S).txt 
            if [ ! -f $backup ]; then 
                mv $smartSlurmJobRecordDir/jobRecord.txt $backup 2>/dev/null
                echo Too many fails in the last hour for $software $ref, backing jobRecord.txt to# $backup and removing failed + failure-disproven records for $software $ref, validated history kept. 2>/dev/null   
                awk -F"," -v a=$2 -v b=$3 'NR==FNR{ if($12==a && $13==b && $9!="COMPLETED"){ if($9 ~ /OOM/){ if(!mN++ || $2<fsM)fsM=$2; if($5>fmM)fmM=$5 } if($9 ~ /OOT/){ if(!tN++ || $2<fsT)fsT=$2; if($6>ftM)ftM=$6 } } next } { if(!($12==a && $13==b)){ print; next } if($9!="COMPLETED") next; s=0; if(mN && $2>=fsM && $7<=fmM) s=1; if(tN && $2>=fsT && $8<=ftM) s=1; if(!s) print }' $backup $backup > $smartSlurmJobRecordDir/jobRecord.txt 
                #ls -l $smartSlurmJobRecordDir/ 
                #ls -l $smartSlurmJobRecordDir/stats/
                rm $smartSlurmJobRecordDir/stats/$software.$ref.* 2>/dev/null
            fi
            # email user to review jobRecord.txt and related stat files, because there are too many fails in the last hour for this software and ref, which might indicate some systemic issue with the job or the cluster
            echo -e "There are $recentFails fails in the last hour for $software $ref, which might indicate some systemic issue with the job or the cluster. All job records for this software and ref  have been backed up to $backup; failed and failure-disproven records were removed, validated COMPLETED history kept in jobRecord.txt." | mail -s "!!!Too many oom/oot/over-estimates for $software $ref" $USER

        fi

fi

#if [[ x == y ]]; then 

# all job plots -- gated OFF by default. Set ssMakePlots=yes to restore the
# gnuplot barcharts and per-job mem/cpu pngs. Efficiency now lives in
# efficiency.tsv + `checkRun` (the `e` key / `checkRunGrid --efficiency`), so by
# default no per-job pngs are generated and nothing crowds smartSlurmLog.
if [ "${ssMakePlots:-}" = yes ]; then
tm=$SLURM_JOBID #  `mktemp XXXXXXXX --dry-run`

#rm $smartSlurmLogDir/barchartMem.png  $smartSlurmLogDir/barchartTime.png 2>/dev/null
echo Category,Used,Wasted,Saved2,default,Saved1 > $smartSlurmLogDir/$tm.dataMem.csv

if [ -f .command.sh ] && [ -f .command.run ]; then 
    ls -cr $smartSlurmLogDir/../../../*/*/*/*.out | xargs -d '\n' grep ^dataToPlot | awk -F, '{printf "%s-%s-%s,%s,%s,%s,%s\n", substr($15,1,index($15,".")-1), substr($2, length($2)-3), substr($10,1,3),  $8 + $17 *2,   $6-$8-$17 *2, $4-$6, $4}' | sed s/-COM//g | sed s/-OO/-/g >> $smartSlurmLogDir/$tm.dataMem.csv
else 
    ls $smartSlurmLogDir/*.out | sort -n | xargs -d '\n' grep ^dataToPlot | awk -F, '{printf "%s-%s-%s,%s,%s,%s,%s\n", substr($15,1,index($15,".")-1), substr($2, length($2)-3), substr($10,1,3),  $8 + $17 *2,   $6-$8-$17 *2, $4-$6, $4}' | sed s/-COM//g | sed s/-OO/-/g >> $smartSlurmLogDir/$tm.dataMem.csv
fi    

# if less than 0, change to zeor
awk -F',' -v OFS=',' '{ for (i=1; i<=NF; i++) if ($i < 0) $i = 0; print }' $smartSlurmLogDir/$tm.dataMem.csv > $smartSlurmLogDir/$tm.output.csv

awk -F, -v OFS=',' -v max=$(awk -F, 'BEGIN {max=0} {if (NR!=1 && $5>max) max=$5} END {print max}' $smartSlurmLogDir/$tm.output.csv) '{if(NR==1) print $0; else {diff=max-$5; print $0 "," diff "," max}}' $smartSlurmLogDir/$tm.output.csv > $smartSlurmLogDir/$tm.dataMem.csv

width=`wc -l $smartSlurmLogDir/$tm.dataMem.csv | cut -d' ' -f1`; width=$((width*16)); [ $width -le 800 ] && width=800; 
gnuplot -e "set key outside; set key reverse; set key invert; set datafile separator ','; set style data histogram; set style histogram rowstacked gap 2; set style fill solid border rgb 'black'; set xtics rotate by -45; set terminal png size $width,600; set output '$smartSlurmLogDir/$tm.barchartMem.png'; set title 'Job vs. Memmory'; set ylabel 'Memory (MegaBytes)'; plot '$smartSlurmLogDir/$tm.dataMem.csv' using 2:xtic(1) title 'Used' lc rgb 'green', '' using 3:xtic(1) title 'Wasted' lc rgb 'red', '' using 4:xtic(1) title 'Saved2' lc rgb 'yellow', '' using 6:xtic(1) title 'Saved1' lc rgb 'pink'"

echo To see the plot:
echo display $smartSlurmLogDir/$tm.barchartMem.png

echo Category,Used,Wasted,default,Saved > $smartSlurmLogDir/$tm.dataTime.csv

ls $smartSlurmLogDir/*.out | sort -n | xargs -d '\n' grep ^dataToPlot | awk -F, '{printf "%s-%s-%s,%s,%s,%s,%s\n", substr($15,1,index($15,".")-1), substr($2, length($2)-3), substr($10,1,3),  $9*$16,   ($7-$9)*$16, ($5-$7)*$16, $5*$16}' | sed s/-COM//g | sed s/-OO/-/g >> $smartSlurmLogDir/$tm.dataTime.csv

awk -F',' -v OFS=',' '{ for (i=1; i<=NF; i++) if ($i < 0) $i = 0; print }' $smartSlurmLogDir/$tm.dataTime.csv > $smartSlurmLogDir/$tm.output.csv 

awk -F, -v OFS=',' -v max=$(awk -F, 'BEGIN {max=0} {if (NR!=1 && $5>max) max=$5} END {print max}' $smartSlurmLogDir/$tm.output.csv) '{if(NR==1) print $0; else {diff=max-$5; print $0 "," diff "," max}}' $smartSlurmLogDir/$tm.output.csv > $smartSlurmLogDir/$tm.dataTime.csv

#width=`wc -l $smartSlurmLogDir/$tm.dataTime.csv | cut -d' ' -f1`; width=$((width*16)); [ $width -le 800 ] && width=800; 
# all job time
gnuplot -e "set key outside; set key reverse; set key invert; set datafile separator ','; set style data histogram; set style histogram rowstacked gap 2; set style fill solid border rgb 'black'; set xtics rotate by -45; set terminal png size $width,600; set output '$smartSlurmLogDir/$tm.barchartTime.png'; set title 'Job vs. Time'; set ylabel 'Time (Mins)'; plot '$smartSlurmLogDir/$tm.dataTime.csv' using 2:xtic(1) title 'Used' lc rgb 'green', '' using 3:xtic(1) title 'Wasted' lc rgb 'red', '' using 4:xtic(1) title 'Saved' lc rgb 'yellow'" #", '' using 6:xtic(1) title 'Saved' lc rgb 'pink'"

echo To see the plot:
echo display $smartSlurmLogDir/$tm.barchartTime.png

#rm $smartSlurmLogDir/$tm.*.csv


# todo: should make the plot wider instead of shink it: 
# https://stackoverflow.com/questions/13869439/gnuplot-how-to-increase-the-width-of-my-graph

# rowTotal=`wc -l $smartSlurmLogDir/job_$SLURM_JOBID.memCPU.txt | cut -d' ' -f1`
# if [ "$rowTotal" -gt 50 ]; then 
#     maxMem=0; maxCpu=0; 
#     rate=`echo "scale=2;$rowTotal/50"|bc`
#     IFS=$'\n'; rowCount1=0; rowCount2=0
#     echo > $smartSlurmLogDir/job_$SLURM_JOBID.memCPU1.txt
#     for t in `cat $smartSlurmLogDir/job_$SLURM_JOBID.memCPU.txt`; do
#         mem=`echo $t | cut -d' ' -f2`
#         cpu=`echo $t | cut -d' ' -f5`
#         [ "$mem" -gt $maxMem ] && maxMem=$mem && mem1=`echo $t | cut -d' ' -f3,4`
#         [ "$cpu" -gt $maxCpu ] && maxCpu=$cpu
#         rowCount1=$((rowCount1 + 1))
#         rowMax=`echo "scale=2;$rowCount2*$rate"|bc`
#         rowMax=${rowMax%.*}; [ -z "$rowMax" ] && rowMax=1; 
#         if [ "$rowMax" -le "$rowCount1" ]; then 
#             rowCount2=$((rowCount2 + 1))
#             echo $rowCount2 $maxMem $mem1 $maxCpu >> $smartSlurmLogDir/job_$SLURM_JOBID.memCPU1.txt
#             maxMem=0; maxCpu=0;
#         fi 
#     done
# else 
#     cp $smartSlurmLogDir/job_$SLURM_JOBID.memCPU.txt $smartSlurmLogDir/job_$SLURM_JOBID.memCPU1.txt
# fi 

# todo: need adjust reserved memory if the job was ajusted by upstream job!!!!!!!
# time vs. memory for current job
width=`wc -l $smartSlurmLogDir/job_$SLURM_JOBID.memCPU.txt | cut -d' ' -f1`; width=$((width*16)); [ $width -le 800 ] && width=800
gnuplot -e "set key outside; set key reverse; set key invert; set datafile separator ' '; set style data histogram; set style histogram rowstacked gap 2; set style fill solid border rgb 'black'; set xtics rotate by -45; set terminal png size $width,600; set output '$smartSlurmLogDir/job_$SLURM_JOBID.mem.png'; set title 'Time vs. Mem for job $SLURM_JOBID'; set xlabel 'Time (min)'; set ylabel 'Mem (M)'; plot '$smartSlurmLogDir/job_$SLURM_JOBID.memCPU.txt' using 2:xtic(1) title 'Used' lc rgb 'green', '' using 3:xtic(1) title 'Wasted' lc rgb 'red', '' using 4:xtic(1) title 'Saved' lc rgb 'yellow'"

# time vs. CPU usage for current job
width=`wc -l $smartSlurmLogDir/$tm.dataMem.csv | cut -d' ' -f1`; width=$((width*16)); [ $width -le 800 ] && width=800
gnuplot -e "set key outside; set key reverse; set key invert; set datafile separator ' '; set style data histogram; set style histogram rowstacked gap 2; set style fill solid border rgb 'black'; set xtics rotate by -45; set terminal png size $width,600; set output '$smartSlurmLogDir/job_$SLURM_JOBID.cpu.png'; set title 'Time vs. CPU Usage for job $SLURM_JOBID'; set xlabel 'Time (min)'; set ylabel 'CPU Usage (%)'; plot '$smartSlurmLogDir/job_$SLURM_JOBID.memCPU.txt' using 5:xtic(1) title 'Used' lc rgb 'green'"

#fi


# echo "$counter $total_memory_usage $(($reservedMem - $total_memory_usage)) $saved ${total_cpu_usage%.*}" >> job_$SLURM_JOBID.memCPU.txt

fi   # end ssMakePlots (per-job/barchart pngs gated off by default)

rate=0.0013 # $0.0013 per G per hour
maxSaved=`cut -d' ' -f4 $smartSlurmLogDir/job_$SLURM_JOBID.memCPU.txt | sort -n | tail -n1`
jMin=`wc -l $smartSlurmLogDir/job_$SLURM_JOBID.memCPU.txt | cut -d ' ' -f 1`
savedDollar=$(echo "$rate * $maxSaved / 1024 * $jMin / 60" | bc -l)

maxUsed=`cut -d' ' -f2 $smartSlurmLogDir/job_$SLURM_JOBID.memCPU.txt | sort -n | tail -n1`
#jMin=`wc -l $smartSlurmLogDir/job_$SLURM_JOBID.memCPU.txt | cut -d ' ' -f 1`
savedDollar1=$(echo "$rate * ($defaultMem - $maxUsed) / 1024 * $jMin / 60" | bc -l)

#echo Category,Used,Wasted,Saved2,default,Saved1 > $smartSlurmLogDir/$tm.dataMem.csv
# The per-run savings aggregate is read from the plot CSV, which is only built
# when ssMakePlots=yes. With plots off (the default) that file does not exist, so
# compute these only when it is present -- otherwise savedMem is empty and bc/awk
# would error and blank the "SmartSlurm Saved ... for this run" line. The per-run
# efficiency email already summarizes used-vs-allocated for the whole run.
if [ -f "$smartSlurmLogDir/$tm.dataMem.csv" ]; then
    savedMem=`awk -F',' '{sum += $4} END {print sum}' $smartSlurmLogDir/$tm.dataMem.csv`
    savedMem1=`awk -F',' '{sum += $6} END {print sum}' $smartSlurmLogDir/$tm.dataMem.csv`
    savedDollar3=$(echo "$rate * $savedMem / 1024 * $jMin / 60" | bc -l)
    savedDollar4=$(echo "$rate * $savedMem1 / 1024 * $jMin / 60" | bc -l)
fi


minimumsize=9000
actualsize=`wc -c $out || echo 0`

[ -f $succFile ] && s="${overEstimate}Succ:$SLURM_JOBID:$SLURM_JOB_NAME" || s="$jobStatus:$SLURM_JOBID:$SLURM_JOB_NAME"

if [ "${actualsize% *}" -ge "$minimumsize" ]; then
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

[ -n "$savedDollar3" ] && toSend="So far, SmartSlurm Saved $savedDollar3 by dynamic allocation of memory for this run.\n$toSend"
[ -n "$savedDollar4" ] && [ -z "$snakemakeSuccFlag" ] && toSend="So far, SmartSlurm Saved $savedDollar4 by split workflow into steps from this run. \n$toSend"

if [ -f "$err" ]; then
    actualsize=`wc -c $err`
    if [ "${actualsize% *}" -ge "$minimumsize" ]; then
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
[ -z "$userEmail" ] || USER=$userEmail

# Optional DAG render -- ON DEMAND only. checkRun's `g` key renders the workflow
# graph live, so by default we no longer auto-render a png here. Set
# ssMakePlots=yes to restore the old behavior.
if [ "${ssMakePlots:-}" = yes ] && [ -z "$snakemakeSuccFlag" ] && [ -n "$condaEnvPath" ]; then
    module load conda/miniforge3/24.11.3-0 2>/dev/null
    conda activate "$condaEnvPath" 2>/dev/null
    workflowPlot.py "$smartSlurmLogDir" && dot -Tsvg "$smartSlurmLogDir/jobs.dot" -o "$smartSlurmLogDir/dag.svg" && convert "$smartSlurmLogDir/dag.svg" "$smartSlurmLogDir/dag.png"
fi

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
