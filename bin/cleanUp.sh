#!/bin/bash

#set -x

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
sleep 1

sacct=`sacct --format=JobID,Submit,Start,End,MaxRSS,State,NodeList%30,Partition,ReqTRES%30,TotalCPU,Elapsed%14,Timelimit%14,ExitCode --units=M -j $SLURM_JOBID`

#sacct=`cat ~/fakeSacct.txt`

#sacct report is not accurate, let us use the total memory
#totalM=${sacct#*,mem=}; totalM=${totalM%%M,n*}

echo -e "\nJob summary:\n$sacct"
# echo *Notice the sacct report above: while the main job is still running for sacct command, user task is completed.

# record job for future estimating mem and time
jobStat=`echo -e "$sacct" | tail -n 1`

#from: "sacct --format=JobID,Submit,Start,End,MaxRSS,State,NodeList%30,Partition,ReqTRES%30,TotalCPU,Elapsed%14,Timelimit%14 --units=M -j $SLURM_JOBID"

START=`head -n 1 $smartSlurmLogDir/job_$SLURM_JOB_ID.memCPU.txt | cut -d' ' -f6`
FINISH=`date +%s`

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
        scontrol show jobid -dd $SLURM_JOB_ID
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

#set -x
# totalM=200; totalT=200; srunM=1000; min=200;
overReserved=""; overText=""; ratioM=""; ratioT=""
if [ "$memDef" -ne "$totalM" ] || [ "$minDef" -ne "$totalT" ]; then
  if [[ "$jobStatus" == COMPLETED ]]; then
    if [ "$totalM" -gt $((srunM * 2)) ] && [ $srunM -ge 1024 ]; then
        overReserved="O:"
        overText="$SLURM_JOBID over-reserved resounce M: $srunM/$totalM T: $min/$totalT"
	#echo -e "`pwd`\n$out\n$SLURM_JOBID over-reserved resounce M: $srunM/$totalM T: $min/$totalT" | mail -s "Over reserved $SLURM_JOBID" $USER
        echo $overText
    fi
  fi
fi
ratioM=`echo "scale=2;$srunM/$totalM"|bc`; ratioT=`echo "scale=2;$min/$totalT"|bc`

# todo: move this part to main job, so that when release job, this job record can be used for the statics 
                                #3defult,  5given,  7cGroupUsed                  sacct used
record="$SLURM_JOB_ID,$inputSize,$memDef,$minDef,$totalM,$totalT,$srunM,$min,$jobStatus,$USER,$memSacct,$2,$ref,$flag,$core,$extraMemC,$defaultExtraTime,0$ratioM,0$ratioT,`date`" 
echo dataToPlot,$record

if [[ $jobStatus == "COMPLETED" ]]; then # && [[ "${skipEstimate}" == n ]]; then
        records=`grep COMPLETED $smartSlurmJobRecordDir/jobRecord.txt 2>/dev/null | awk -F"," -v a=$2 -v b=$3 '{ if($12 == a && $13 == b) {print $2, $7 }}' | sort -u -n`
        #timeRecords=`grep COMPLETED $smartSlurmJobRecordDir/jobRecord.txt 2>/dev/null | awk -F"," -v a=$2 -v b=$3 '{ if($12 == a && $13 == b) {print $8 }}' | sort -u -nr | tr '\n' ' '`
        #memRecords=`grep COMPLETED $smartSlurmJobRecordDir/jobRecord.txt 2>/dev/null | awk -F"," -v a=$2 -v b=$3 '{ if($12 == a && $13 == b) {print $7 }}' | sort -u -nr | tr '\n' ' '`
        #timeRecords=`grep COMPLETED $smartSlurmJobRecordDir/jobRecord.txt 2>/dev/null | awk -F"," -v a=$2 -v b=$3 '{ if($12 == a && $13 == b) {print $8 }}' | sort -u -nr | tr '\n' ' '`

        #maxMem=`cat $smartSlurmJobRecordDir/stats/$software.$ref.mem.stat.noInput  2>/dev/null | sort -nr | tr '\n' ' ' | cut -f 1 -d " "`
        #maxTime=`cat $smartSlurmJobRecordDir/stats/$software.$ref.time.stat.noInput  2>/dev/null | sort -nr | tr '\n' ' ' | cut -f 1 -d " "
        #if [ "$(echo $memRecords | wc -w)" -lt 20 ] || [ "${memRecords%% *}" -lt "${srunM%.*}" ] || [ "${timeRecords%% *}" -lt "$min" ]; then

        #   less than 20 records                  # or current one is larger than all old data   # and not exist already
        if [ "$(echo $records | wc -l)" -lt 200 ] || [ "`echo $records | tail -n1 | cut -d' ' -f1`" -lt "$inputSize" ] && ! echo $records | grep "$inputSize $srunM"; then
            #if [ ! -f ${out%.out}.startFromCheckpoint ]; then
                [ -z "$START" ] || echo $record >> $smartSlurmJobRecordDir/jobRecord.txt
                echo -e "Added this line to $smartSlurmJobRecordDir/jobRecord.txt:\n$record"
            #else
            #    echo -e "Has ${out%.out}.startFromCheckpoint. Did not added this line to $smartSlurmJobRecordDir/jobRecord.txt:\n$record"
            #fi
            mv $smartSlurmJobRecordDir/stats/$software.$ref.* $smartSlurmJobRecordDir/stats/back 2>/dev/null
        else
            echo Did not add this record to $smartSlurmJobRecordDir/jobRecord.txt
        fi
        #cat $smartSlurmJobRecordDir/jobRecord.txt
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
        echo -e "This step might not good to checkpoint?\nIt might because you did not give enouther memory. Please test run it with higher memmory:\n$flag" | mail -s "!!!!$SLURM_JOBID:$flag" $USER
    fi
    touch $failFile
    # if checkpoint due to low memory or checkpoint failed due to memory, or actually out of memory
    if  [ -f "${out%.out}.likelyCheckpointOOM" ] || [[ "$jobStatus" == "OOM" ]]; then

        jobStatus=OOM

        #set -x
        echo Will re-queue after sending email...

        #if [[ "${skipEstimate}" == n ]]; then
            echo old extraMem:
            cat $smartSlurmJobRecordDir/stats/extraMem.$software.$ref

            extraMemN=$(( totalM - srunM + 1 ))

            maxExtra=`sort -n $smartSlurmJobRecordDir/stats/extraMem.$software.$ref | tail -n1 | cut -d' ' -f1`
            [ -z "$maxExtra" ] && maxExtra=5
            [ $extraMemN -gt $maxExtra ] && echo $extraMemN $totalM $inputSize $SLURM_JOBID >> $smartSlurmJobRecordDir/stats/extraMem.$software.$ref && maxExtra=$extraMemN

            echo new extraMem:
            cat $smartSlurmJobRecordDir/stats/extraMem.$software.$ref
        #else
        #    maxExtra=5
        #fi

        mv ${out%.out}.likelyCheckpointOOM ${out%.out}.likelyCheckpointOOM.old 2>/dev/null
        (
        #set -x;
        if [ $totalM -lt 100 ]; then
            totalM=100
            newFactor=8
        elif [ $totalM -lt 200 ]; then
            totalM=200
            newFactor=5
        elif [ $totalM -lt 512 ]; then
            totalM=512
            newFactor=4
        elif [ $totalM -lt 1024 ]; then
            totalM=1024
            newFactor=3
        elif [ $totalM -lt 10240 ]; then
            newFactor=2
        elif [ $totalM -lt 51200 ]; then
            newFactor=1.5
        else
            newFactor=1.2
        fi



        mem=`echo "($totalM*$newFactor*1.2+$maxExtra*2)/1" | bc`
        
        # testing here
        #mem=1000

        echo Trying to requeue with $mem M
        echo $inputSize $mem $totalT $maxExtra > ${out%.out}.adjust

        #[[ "$USER" == ld32 ]] && hostName=login00 || hostName=o2.hms.harvard.edu
        #set -x
        #if `ssh $hostName "scontrol requeue $SLURM_JOBID; scontrol update JobId=$SLURM_JOBID MinMemoryNode=$mem;"`; then
        #echo "scontrol requeue $SLURM_JOBID; sleep 2; scontrol update JobId=$SLURM_JOBID MinMemoryNode=$mem;" > $smartSlurmLogDir/$flag.requeueCMD

        hours=$((($totalT + 59) / 60))
        adjustPartition $hours $partition

        export myPartition=$partition
        export myTime=$totalT
        export myMem=${mem}M
        requeueCmd=`grep "Command used to submit the job:" $script | tail -n 1`
        requeueCmd=${requeueCmd#*submit the job: }
        requeueCmd=${requeueCmd/ -H/}
        requeueCmd=${requeueCmd//\$myPartition/$myPartition}
        requeueCmd=${requeueCmd//\$myTime/$myTime}
        requeueCmd=${requeueCmd//\$myMem/$myMem}

		requeueCmd=$( echo $requeueCmd | sed -E 's/-d afterok:[0-9]+//g' ) #-d afterok:38023271
        for try in {1..5}; do
            if [ ! -f $failFile.requeued.$try.mem ]; then
                #sleep 2
                touch $failFile.requeued.$try.mem
                
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

                    echo Requeued successfully
                    [ -f $failFile ] && rm $failFile

                    s="OOM.Requeued:$SLURM_JOBID-$newJobID:$SLURM_JOB_NAME"
                    echo -e "" | mail -s "$s" $USER
                    [[ $USER != ld32 ]] && echo -e "" | mail -s "$s" ld32

                    touch $out.$newJobID

                    cp $smartSlurmLogDir/allJobs.txt $smartSlurmLogDir/allJobs.requeue.$SLURM_JOB_ID.as.$newJobID

                    while true; do
                        if `mkdir $smartSlurmLogDir/job.adjusting.lock  2>/dev/null`; then
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
                set +x
            fi
        done;
        set +x;
        ) &

        # delete stats and redo them
        if [ "$inputs" == "none" ]; then
            mv $smartSlurmJobRecordDir/stats/$software.$ref.* $smartSlurmJobRecordDir/stats/back  2>/dev/null
        else
            # remove bad records
            # todo: all need to delete some resords form jobRecords.txt!!!
            if [ -f $smartSlurmJobRecordDir/stats/$software.$ref.mem.stat ]; then
                # .  $smartSlurmJobRecordDir/stats/$software.$ref.mem.stat
                # Finala=`printf "%.15f\n" $Finala`
                # Finalb=`printf "%.15f\n" $Finalb`
                # Maximum=`printf "%.15f\n" $Maximum`
                # echo Finala: $Finala Finalb: $Finalb Maximum: $Maximum

                # awk -F"," -v a=$2 -v b=$3 -v c=$Finala -v d=$Finalb '{ if ( ! ($12 == a && $13 == b && $2 * c + d > $7)) {print}}' $smartSlurmJobRecordDir/jobRecord.txt > $smartSlurmJobRecordDir/jobRecord.txt1
                # echo diff output:
                # diff $smartSlurmJobRecordDir/jobRecord.txt $smartSlurmJobRecordDir/jobRecord.txt1
                # mv $smartSlurmJobRecordDir/jobRecord.txt1 $smartSlurmJobRecordDir/jobRecord.txt
                mv $smartSlurmJobRecordDir/stats/$software.$ref.* $smartSlurmJobRecordDir/stats/back  2>/dev/null
            fi
        fi
    elif [[ "$jobStatus" == "OOT" ]]; then
#set -x
        if [ $min -lt 20 ]; then
            min=20
            newFactor=2
        elif [ $min -lt 30 ]; then
            min=30
            newFactor=4
        elif [ $min -lt 60 ]; then
            min=60
            newFactor=3
        elif [ $min -lt 120 ]; then
            min=120
            newFactor=2
        elif [ $min -lt 480 ]; then
            min=480
            newFactor=1.5
        else
            newFactor=1.2
        fi
 
        #newFactor=2
        min=`echo "($min*$newFactor*1.2)/1" | bc`
        echo Trying to requeue $try with $min minutes
        echo $inputSize $totalM $min $maxExtra > ${out%.out}.adjust

        hours=$((($min + 59) / 60))
        adjustPartition $hours $partition
        seconds=$(($min * 60))

        # https://chat.openai.com/chat/80e28ff1-4885-4fe3-8f21-3556d221d7c6

        time=`eval "echo $(date -ud "@$seconds" +'$((%s/3600/24))-%H:%M:%S')"`
        export myPartition=$partition
        export myTime=$time
        export myMem=${totalM}M
        requeueCmd=`grep "Command used to submit the job:" $script | tail -n 1`
        requeueCmd=${requeueCmd#*submit the job: }
        requeueCmd=${requeueCmd//\$myPartition/$myPartition}
        requeueCmd=${requeueCmd//\$myTime/$myTime}
        requeueCmd=${requeueCmd//\$myMem/$myMem}
        requeueCmd=${requeueCmd/ -H/}
        requeueCmd=$( echo $requeueCmd | sed -E 's/-d afterok:[0-9]+//g' ) #-d afterok:38023271

        for try in {1..5}; do
            if [ ! -f $failFile.requeued.$try.time ]; then
                touch $failFile.requeued.$try.time
                
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

                    echo Requeued successfully
                    [ -f $failFile ] && rm $failFile

                    s="OOT.Requeued:$SLURM_JOBID-$newJobID:$SLURM_JOB_NAME"
                    echo -e "" | mail -s "$s" $USER
                    [[ $USER != ld32 ]] && echo -e "" | mail -s "$s" ld32

                    touch $out.$newJobID

                    cp $smartSlurmLogDir/allJobs.txt $smartSlurmLogDir/allJobs.requeue.$SLURM_JOB_ID.as.$newJobID
                    while true; do
                        if `mkdir $smartSlurmLogDir/job.adjusting.lock  2>/dev/null`; then
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
        set +x 

         # delete stats and redo them
        if [[ "$inputs" == "none" ]]; then
            mv $smartSlurmJobRecordDir/stats/$software.$ref.* $smartSlurmJobRecordDir/stats/back  2>/dev/null

            #[ -f $smartSlurmJobRecordDir/stats/$software.$ref.mem.stat.noInput ] && mv $smartSlurmJobRecordDir/stats/$software.$ref.mem.stat.noInput $smartSlurmJobRecordDir/stats/$software.$ref.mem.stat.noInput.$(stat -c %y $smartSlurmJobRecordDir/stats/$software.$ref.mem.stat.noInput | tr " " ".")
        else
             # remove bad records
             # todo: all need to delete some resords form jobRecords.txt!!!
            if [ -f $smartSlurmJobRecordDir/stats/$software.$ref.time.stat ]; then
                # .  $smartSlurmJobRecordDir/stats/$software.$ref.time.stat
                # Finala=`printf "%.15f\n" $Finala`
                # Finalb=`printf "%.15f\n" $Finalb`
                # Maximum=`printf "%.15f\n" $Maximum`
                # echo Finala: $Finala Finalb: $Finalb Maximum: $Maximum

                # awk -F"," -v a=$2 -v b=$3 -v c=$Finala -v d=$Finalb '{ if ( ! ($12 == a && $13 == b && $2 * c + d > $8 ) ) {print}}' $smartSlurmJobRecordDir/jobRecord.txt > $smartSlurmJobRecordDir/jobRecord.txt1
                # echo diff output:
                # diff $smartSlurmJobRecordDir/jobRecord.txt $smartSlurmJobRecordDir/jobRecord.txt1
                # mv $smartSlurmJobRecordDir/jobRecord.txt1 $smartSlurmJobRecordDir/jobRecord.txt
                mv $smartSlurmJobRecordDir/stats/$software.$ref.* $smartSlurmJobRecordDir/stats/back  2>/dev/null
            fi
        fi

    else
        echo Not sure why job failed. Not run out of time or memory. Pelase check youself.
        
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

                        echo Requeued successfully
                        [ -f $failFile ] && rm $failFile

                        s="FailBadNode.Requeued:$SLURM_JOBID-$newJobID:$SLURM_JOB_NAME"
                        echo -e "" | mail -s "$s" $USER
                        [[ $USER != ld32 ]] && echo -e "" | mail -s "$s" ld32

                        touch $out.$newJobID

                        cp $smartSlurmLogDir/allJobs.txt $smartSlurmLogDir/allJobs.requeue.$SLURM_JOB_ID.as.$newJobID
                        while true; do
                            if `mkdir $smartSlurmLogDir/job.adjusting.lock  2>/dev/null`; then
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

fi

#if [[ x == y ]]; then 

# all job plots
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


# echo "$counter $total_memory_usage $(($reservedMem - $total_memory_usage)) $saved ${total_cpu_usage%.*}" >> job_$SLURM_JOB_ID.memCPU.txt

rate=0.0013 # $0.0013 per G per hour 
maxSaved=`cut -d' ' -f4 $smartSlurmLogDir/job_$SLURM_JOBID.memCPU.txt | sort -n | tail -n1`
jMin=`wc -l $smartSlurmLogDir/job_$SLURM_JOBID.memCPU.txt | cut -d ' ' -f 1`
savedDollar=$(echo "$rate * $maxSaved / 1024 * $jMin / 60" | bc -l)

maxUsed=`cut -d' ' -f2 $smartSlurmLogDir/job_$SLURM_JOBID.memCPU.txt | sort -n | tail -n1`
#jMin=`wc -l $smartSlurmLogDir/job_$SLURM_JOBID.memCPU.txt | cut -d ' ' -f 1`
savedDollar1=$(echo "$rate * ($defaultMem - $maxUsed) / 1024 * $jMin / 60" | bc -l)

#echo Category,Used,Wasted,Saved2,default,Saved1 > $smartSlurmLogDir/$tm.dataMem.csv
savedMem=`awk -F',' '{sum += $4} END {print sum}' $smartSlurmLogDir/$tm.dataMem.csv`
savedMem1=`awk -F',' '{sum += $6} END {print sum}' $smartSlurmLogDir/$tm.dataMem.csv`
savedDollar3=$(echo "$rate * $savedMem / 1024 * $jMin / 60" | bc -l)
savedDollar4=$(echo "$rate * $savedMem1 / 1024 * $jMin / 60" | bc -l)


minimumsize=9000
actualsize=`wc -c $out || echo 0`

[ -f $succFile ] && s="${overReserved}Succ:$SLURM_JOBID:$SLURM_JOB_NAME" || s="$jobStatus:$SLURM_JOBID:$SLURM_JOB_NAME"

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

toSend="So far, SmartSlurm Saved $savedDollar3 by dynamic allocation of memory for this run.\n$toSend"
[ -z "$snakemakeSuccFlag" ] && toSend="So far, SmartSlurm Saved $savedDollar4 by split workflow into steps from this run. \n$toSend"

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

[ -f $smartSlurmLogDir/summary.$SLURMJOB_ID ] && toSend="`cat $smartSlurmLogDir/summary.$SLURMJOB_ID`\n$toSend" && s="${toSend%% *} $s"

#echo -e "tosend:\n$toSend"
#echo -e "$toSend" >> ${err%.err}.email



if [[ $USER != ld32 ]]; then
    if [ -f $smartSlurmJobRecordDir/stats/$software.$ref.mem.png ]; then
        echo -e "$toSend" | mail -s "$s" -a $smartSlurmLogDir/job_$SLURM_JOBID.mem.png -a $smartSlurmLogDir/job_$SLURM_JOBID.cpu.png -a $smartSlurmLogDir/$tm.barchartMem.png -a $smartSlurmLogDir/$tm.barchartTime.png -a $smartSlurmJobRecordDir/stats/$software.$ref.mem.png -a $smartSlurmJobRecordDir/stats/$software.$ref.time.png ld32
    elif [ -f $smartSlurmJobRecordDir/stats/back/$software.$ref.time.png ]; then
        echo -e "$toSend" | mail -s "$s" -a $smartSlurmLogDir/job_$SLURM_JOBID.mem.png -a $smartSlurmLogDir/job_$SLURM_JOBID.cpu.png -a $smartSlurmLogDir/$tm.barchartMem.png -a $smartSlurmLogDir/$tm.barchartTime.png -a $smartSlurmJobRecordDir/stats/back/$software.$ref.mem.png -a $smartSlurmJobRecordDir/stats/back/$software.$ref.time.png ld32
    else
        echo -e "$toSend" | mail -s "$s" -a $smartSlurmLogDir/job_$SLURM_JOBID.mem.png -a $smartSlurmLogDir/job_$SLURM_JOBID.cpu.png -a $smartSlurmLogDir/$tm.barchartMem.png -a $smartSlurmLogDir/$tm.barchartTime.png ld32
    fi
fi

[ -z "$userEmail" ] || USER=$userEmail

if [ -z "$snakemakeSuccFlag" ]; then    
    #set -x 
    module load conda/miniforge3/24.11.3-0
    conda activate smartSlurmEnv
    workflowPlot.py $smartSlurmLogDir &&  dot -Tsvg $smartSlurmLogDir/jobs.dot -o $smartSlurmLogDir/dag.svg && convert $smartSlurmLogDir/dag.svg $smartSlurmLogDir/dag.png
fi 

#set +x 

if [[ "$lessEmail" == "noEmail" ]]; then 
    echo noEmail
elif [[ "$lessEmail" == "noSuccEmail" ]] && [[ "$s" == *Succ* ]]; then
    echo lessEmail: not sending email because job successful
 
elif [ -f $smartSlurmLogDir/dag.png ]; then
    #echo -e "$toSend" | sendmail `head -n 1 ~/.forward`
    if [ -f $smartSlurmJobRecordDir/stats/$software.$ref.mem.png ]; then
        echo -e "$toSend" | mail -s "$s" -a $smartSlurmLogDir/dag.png -a $smartSlurmLogDir/job_$SLURM_JOBID.mem.png -a $smartSlurmLogDir/job_$SLURM_JOBID.cpu.png -a $smartSlurmLogDir/$tm.barchartMem.png -a $smartSlurmLogDir/$tm.barchartTime.png -a $smartSlurmJobRecordDir/stats/$software.$ref.mem.png -a $smartSlurmJobRecordDir/stats/$software.$ref.time.png $USER && echo email sent1 || \
            { echo Email not sent1.; echo -e "Subject: $s\n$toSend" | sendmail `head -n 1 ~/.forward` && echo Email sent11. || echo Email still not sent11; }
    elif [ -f $smartSlurmJobRecordDir/stats/back/$software.$ref.time.png ]; then
        echo -e "$toSend" | mail -s "$s" -a $smartSlurmLogDir/dag.png -a $smartSlurmLogDir/job_$SLURM_JOBID.mem.png -a $smartSlurmLogDir/job_$SLURM_JOBID.cpu.png -a $smartSlurmLogDir/$tm.barchartMem.png -a $smartSlurmLogDir/$tm.barchartTime.png -a $smartSlurmJobRecordDir/stats/back/$software.$ref.mem.png -a $smartSlurmJobRecordDir/stats/back/$software.$ref.time.png $USER && echo email sent2 || \
            { echo Email not sent2.; echo -e "Subject: $s\n$toSend" | sendmail `head -n 1 ~/.forward` && echo Email sent21. || echo Email still not sent21; }
    elif [ -f $smartSlurmLogDir/job_$SLURM_JOBID.mem.png ]; then
        echo -e "$toSend" | mail -s "$s" -a $smartSlurmLogDir/dag.png -a $smartSlurmLogDir/job_$SLURM_JOBID.mem.png -a $smartSlurmLogDir/job_$SLURM_JOBID.cpu.png -a $smartSlurmLogDir/$tm.barchartMem.png -a $smartSlurmLogDir/$tm.barchartTime.png $USER && echo email sent3 || \
        { echo Email not sent3.; echo -e "Subject: $s\n$toSend" | sendmail `head -n 1 ~/.forward` && echo Email sent31. || echo Email still not sent31; }
    else 
        echo -e "$toSend" | mail -s "$s" $USER && echo email sent4 || \
        { echo Email not sent4.; echo -e "Subject: $s\n$toSend" | sendmail `head -n 1 ~/.forward` && echo Email sent by second try41. || echo Email still not sent41; }
    fi
else 
    #echo -e "$toSend" | sendmail `head -n 1 ~/.forward`
    if [ -f $smartSlurmJobRecordDir/stats/$software.$ref.mem.png ]; then
        echo -e "$toSend" | mail -s "$s" -a $smartSlurmLogDir/job_$SLURM_JOBID.mem.png -a $smartSlurmLogDir/job_$SLURM_JOBID.cpu.png -a $smartSlurmLogDir/$tm.barchartMem.png -a $smartSlurmLogDir/$tm.barchartTime.png -a $smartSlurmJobRecordDir/stats/$software.$ref.mem.png -a $smartSlurmJobRecordDir/stats/$software.$ref.time.png $USER && echo email sent1 || \
            { echo Email not sent1.; echo -e "Subject: $s\n$toSend" | sendmail `head -n 1 ~/.forward` && echo Email sent11. || echo Email still not sent11; }
    elif [ -f $smartSlurmJobRecordDir/stats/back/$software.$ref.time.png ]; then
        echo -e "$toSend" | mail -s "$s" -a $smartSlurmLogDir/job_$SLURM_JOBID.mem.png -a $smartSlurmLogDir/job_$SLURM_JOBID.cpu.png -a $smartSlurmLogDir/$tm.barchartMem.png -a $smartSlurmLogDir/$tm.barchartTime.png -a $smartSlurmJobRecordDir/stats/back/$software.$ref.mem.png -a $smartSlurmJobRecordDir/stats/back/$software.$ref.time.png $USER && echo email sent2 || \
            { echo Email not sent2.; echo -e "Subject: $s\n$toSend" | sendmail `head -n 1 ~/.forward` && echo Email sent21. || echo Email still not sent21; }
    elif [ -f $smartSlurmLogDir/job_$SLURM_JOBID.mem.png ]; then
        echo -e "$toSend" | mail -s "$s" -a $smartSlurmLogDir/job_$SLURM_JOBID.mem.png -a $smartSlurmLogDir/job_$SLURM_JOBID.cpu.png -a $smartSlurmLogDir/$tm.barchartMem.png -a $smartSlurmLogDir/$tm.barchartTime.png $USER && echo email sent3 || \
        { echo Email not sent3.; echo -e "Subject: $s\n$toSend" | sendmail `head -n 1 ~/.forward` && echo Email sent31. || echo Email still not sent31; }
    else 
        echo -e "$toSend" | mail -s "$s" $USER && echo email sent4 || \
        { echo Email not sent4.; echo -e "Subject: $s\n$toSend" | sendmail `head -n 1 ~/.forward` && echo Email sent by second try41. || echo Email still not sent41; }
    fi
fi 
sleep 5
