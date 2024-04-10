#!/bin/bash

#set -x

for i in {1..200}; do sleep 1; echo cleanup $i; done & 

# to call this:  0     1         2       3       4          5       6       7      8     9      10         11       12     13     14
#cleanUp.sh          "flag "software" "$ref" "$inputSize" $core   $memO  $timeO   $mem  $time  $partition slurmAcc  inputs extraM extraTime smartSlurmJobRecordDir

#Running /home/ld32/smartSlurm/bin/cleanUp.sh 16.15.cliper.sample1.treatment cliper none 0 12 10240 15 10240 15 short rccg none 5 5

#smartSlurmLogDir=`dirname $1`
flag=$1 #`basename $1`
software=$2
ref=$3
inputSize=$4 # input does not exist when job was submitted.
core=$5
memO=$6
timeO=$7
#totalM=$8

# if jobs has --mem
totalM=$SLURM_MEM_PER_NODE

# if job has --mem-per-cpu and -c
[ -z "$totalM" ] && totalM=$((SLURM_MEM_PER_CPU * SLURM_JOB_CPUS_PER_NODE))

# if job has --mem-per-cpu and -n
[ -z "$totalM" ] &&  totalM=$((SLURM_MEM_PER_CPU * SLURM_CPUS_PER_TASK))

totalT=${9}
partition=${10}
slurmAcc=${11}
inputs=${12}
extraMemC=${13}
extraTime=${14}
#skipEstimate=${15} # todo. remove it
#smartSlurmJobRecordDir=${15}

echo Running $0 $@
#echo pwd: `pwd`

# if not successful, delete cromwell error file, so that the job not labeled fail
#ls -l execution
[ -f $smartSlurmLogDir/$flag.success ] || rm -r execution/rc 2>/dev/null

#cd ${1%log}

# if [ -z "$smartSlurmsmartSlurmJobRecordDir" ]; then
#     if [ -f ~/.smartSlurm/config/config.txt ]; then
#         source ~/.smartSlurm/config/config.txt
#     else
#         source $(dirname $0)/../config/config.txt || { echoerr Config list file not found: config.txt; exit 1; }
#     fi
# fi

#touch /tmp/job_$SLURM_JOB_ID.done
#if [[ -z "$1" ]]; then

   #out=$flag.out; out=${out/\%jerr=${4##* }; err=${err/\%j/$SLURM_JOB_ID}; script=${4% *}; script=${script#* }; succFile=${script/\.sh/}.success;      failFile=${script/\.sh/}.failed;
#    out=slurm-$SLURM_JOBID.out; err=slurm-$SLURM_JOBID.err; script=$flag.sh; succFile=$flag.success; failFile=$flag.failed;
#else
    out=$smartSlurmLogDir/"$flag.out"; err=$smartSlurmLogDir/$flag.err; script=$smartSlurmLogDir/$flag.sh; succFile=$smartSlurmLogDir/$flag.success; failFile=$smartSlurmLogDir/$flag.failed; checkpointDir=$smartSlurmLogDir/$flag
#fi

# if [ -f $succFile ]; then

#     adjustDownStreamJobs.sh $smartSlurmLogDir # todo, make sure to igore if running single job
#     #rm $failFile 2>/dev/null
# fi 


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

#START=`echo $jobStat | cut -d" " -f3`

#START=`date -d "$START" +%s`

START=`head -n 1 $smartSlurmLogDir/job_$SLURM_JOB_ID.memCPU.txt | cut -d' ' -f6`

#FINISH=`echo $jobStat | cut -d" " -f4`

#FINISH=`date -d "$FINISH" +%s`

#[ -z "$FINISH" ] &&
FINISH=`date +%s`

echo  start: $START fisnish: $FINISH
#seconds=$((FINISH - START))

# todo: might use real time directly?
# time is too short, might be wrong result
#[ $seconds -lt 10 ] && echo time is so short. is it right? Let us use real time... && FINISH=`date +%s`

# time in minutes
min=$((($FINISH - $START + 59)/60))

#[[ "$mim" == 0 ]] && min=1

# memory in M
memSacct=`echo $jobStat | cut -d" " -f5`

[[ "$memSacct" == "RUNNING" ]] && memSacct=NA

# node
node=`echo $jobStat | cut -d" " -f7`


# todo: might directly delete this part
#case "$jobStat" in
# jobRecord.txt header
#1user 2software 3ref 4inputName 5inputSizeInK 6CPUNumber 7memoryO 8timeO 9readMem 10RequestedTime 11jobID 12memoryM 13minRun 14Node 15 finalStatus
#*COMPLETED* )  jobStatus="COMPLETED" && echo *Notice the sacct report above: while the main job is still running for sacct command, user task is completed.;;

#*TIMEOUT*   )  jobStatus="OOT";;

#*OUT_OF_ME*   ) jobStatus="OOM";;

#*CANCELLED*	) jobStatus="Cancelled";;

#*FAILED*	) jobStatus="Fail";;

#*          )  jobStatus="Unknown";;


#esac

# for testing
#jobStatus="OOM"

#echo jobStatus: $jobStatus

# sacct might give wrong resules
#[[ $jobStatus != "OOM" ]] && [[ $jobStatus != "OOT" ]] && [[ $jobStatus != "Cancelled" ]] && [ -f $succFile ] && jobStatus="COMPLETED"

#if [ -f $smartSlurmLogDir/$flag/checkpointDoneMem ]; then
#    jobStatus=OOM
#elif [ -f $smartSlurmLogDir/$flag/checkpointDoneTime ]; then
#    jobStatus=OOT
if [ -f $succFile ]; then
    jobStatus=COMPLETED
elif [[ "$sacct" == *TIMEOUT* ]]; then
    jobStatus=OOT
elif [[ "$sacct" == *OUT_OF_ME* ]]; then
    jobStatus=OOM
elif [[ "$sacct" == *FAILED* ]]; then
    jobStatus=Fail
elif [[ "$sacct" == *CANCEL* ]]; then
    jobStatus=Canceled
else
    tLog=`tail -n 22 $out | grep ^srun`
    if [[ "$tLog" == *"task 0: Out Of Memory"* ]]; then
        jobStatus="OOM"
        echo The job is actually out-of-memory by according to the log:
        echo $tLog
        scontrol show jobid -dd $SLURM_JOB_ID
    else
       jobStatus="Unknown"
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
   tText=`cat ${out%.out}.adjust`
   #totalM=`echo $tText | cut -d' ' -f1`
   totalT=`echo $tText | cut -d' ' -f2`
   extraMC=`echo $tText | cut -d' ' -f3`
fi

[ -f $smartSlurmJobRecordDir/stats/extraMem.$software.$ref ] && extraMem=`sort -nr $smartSlurmJobRecordDir/stats/extraMem.$software.$ref | head -n1`

#set -x
# totalM=200; totalT=200; srunM=1000; min=200;
overReserved=""; overText=""; ratioM=""; ratioT=""
if [ "$memO" -ne "$totalM" ] && [ "$totalT" -ne "$totalT" ]; then
  if [[ "$jobStatus" == COMPLETED ]]; then
    if [ "$totalM" -gt $((srunM * 2)) ] && [ $srunM -ge 1024 ]; then
        overReserved="O:"
        overText="$SLURM_JOBID over-reserved resounce M: $srunM/$totalM T: $min/$totalT"
	#echo -e "`pwd`\n$out\n$SLURM_JOBID over-reserved resounce M: $srunM/$totalM T: $min/$totalT" | mail -s "Over reserved $SLURM_JOBID" $USER
        echo $overText
    fi
  fi
  ratioM=`echo "scale=2;$srunM/$totalM"|bc`; ratioT=`echo "scale=2;$min/$totalT"|bc`

fi
#set +x
                                #3defult,  5given,  7cGroupUsed                  sacct used
record="$SLURM_JOB_ID,$inputSize,$memO,$totalT,$totalM,$totalT,$srunM,$min,$jobStatus,$USER,$memSacct,$2,$ref,$flag,$core,$extraMemC,$defaultExtraTime,$ratioM,$ratioT,`date`"  # 16 extraM
echo dataToPlot,$record


#if [[ ! -f $smartSlurmJobRecordDir/stats/$software.$ref.mem.stat || "$2" == "regularSbatch" ]]; then

if [[ $jobStatus == "COMPLETED" ]]; then # && [[ "${skipEstimate}" == n ]]; then

    #if [[ "$inputSize" == 0 ]]; then # || "$2" == "regularSbatch" ]] ; then
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
                echo $record >> $smartSlurmJobRecordDir/jobRecord.txt
                echo -e "Added this line to $smartSlurmJobRecordDir/jobRecord.txt:\n$record"
            #else
            #    echo -e "Has ${out%.out}.startFromCheckpoint. Did not added this line to $smartSlurmJobRecordDir/jobRecord.txt:\n$record"
            #fi
            mv $smartSlurmJobRecordDir/stats/$software.$ref.* $smartSlurmJobRecordDir/stats/back 2>/dev/null
        else
            echo Did not add this record to $smartSlurmJobRecordDir/jobRecord.txt

        fi
    # else
    #     maxMem=`cat $smartSlurmJobRecordDir/stats/$software.$ref.mem.txt  2>/dev/null | cut -f 2 -d " " | sort -nr | tr '\n' ' ' | cut -f 1 -d ' '`

    #     maxTime=`cat $smartSlurmJobRecordDir/stats/$software.$ref.time.txt  2>/dev/null | cut -f 2 -d " " | sort -nr | tr '\n' ' ' | cut -f 1 -d ' '`

    #     if [ -z "$maxMem" ] || [ "${maxMem%.*}" -lt "${srunM%.*}" ] || [ -z "$maxTime" ] || [ "$maxTime" -lt "$min" ]; then
    #         echo $record >> $smartSlurmJobRecordDir/jobRecord.txt

    #         echo -e "Added this line to $smartSlurmJobRecordDir/jobRecord.txt:\n$record"
    #         rm $smartSlurmJobRecordDir/stats/$software.$ref.*  2>/dev/null

    #     else
    #         echo Did not add this record to $smartSlurmJobRecordDir/jobRecord.txt
    #     fi
    # fi
#else
    # todo: may not need failed job records?
#    echo $record >> $smartSlurmJobRecordDir/jobRecord.txt
    #echo -e "Added this line to $smartSlurmJobRecordDir/jobRecord.txt:\n$record"
fi
echo Final mem: $srunM M, time: $min mins


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

            #extraMemN=$(( ( totalM - srunM ) *2 ))
            extraMemN=$(( totalM - srunM + 1 ))
            #[[ "$extraMemN" == 0 ]] && extraMemN=1
            #[ ! -f "${out%.out}.likelyCheckpointOOM" ] &&
            [ $extraMemN -gt 0 ] && echo $extraMemN $totalM $inputSize $SLURM_JOBID >> $smartSlurmJobRecordDir/stats/extraMem.$software.$ref

            #oomCount=`wc -l $smartSlurmJobRecordDir/stats/extraMem.$software.$ref | cut -d' ' -f1`
            maxExtra=`sort -n $smartSlurmJobRecordDir/stats/extraMem.$software.$ref | tail -n1 | cut -d' ' -f1`
            [ -z "$maxExtra" ] && maxExtra=5
            echo new extraMem:
            cat $smartSlurmJobRecordDir/stats/extraMem.$software.$ref
        #else
        #    maxExtra=5
        #fi

        mv ${out%.out}.likelyCheckpointOOM ${out%.out}.likelyCheckpointOOM.old 2>/dev/null
        (
        #set -x;
        for try in {1..5}; do
            if [ ! -f $failFile.requeued.$try.mem ]; then
                #sleep 2
                touch $failFile.requeued.$try.mem
                if [ $totalM -lt 200 ]; then
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

                #newFactor=2
                mem=`echo "($totalM*$newFactor+$maxExtra*2)/1" | bc`
                echo trying to requeue $try with $mem M
                echo $mem $totalT $maxExtra > ${out%.out}.adjust

                #[[ "$USER" == ld32 ]] && hostName=login00 || hostName=o2.hms.harvard.edu
                #set -x
                #if `ssh $hostName "scontrol requeue $SLURM_JOBID; scontrol update JobId=$SLURM_JOBID MinMemoryNode=$mem;"`; then
                #echo "scontrol requeue $SLURM_JOBID; sleep 2; scontrol update JobId=$SLURM_JOBID MinMemoryNode=$mem;" > $smartSlurmLogDir/$flag.requeueCMD

                hours=$((($totalT + 59) / 60))
                adjustPartition $hours $partition

                export  myPartition=$partition
                export myTime=$totalT
                export myMem=${mem}M
                requeueCmd=`grep "Command used to submit the job:" $script | tail -n 1`
                requeueCmd=${requeueCmd#*submit the job: }
                requeueCmd=${requeueCmd//\$myPartition/$myPartition}
                requeueCmd=${requeueCmd//\$myTime/$myTime}
                requeueCmd=${requeueCmd//\$myMem/$myMem}
                newJobID=`$requeueCmd`

                if [[ "$newJobID" =~ ^[0-9]+$ ]]; then
                    echo "# mem=$myMem time=$myTime " >> $script
                    IFS=$'\n'
                    for line in `grep $SLURM_JOBID $smartSlurmLogDir/allJobs.txt | grep -v ^$SLURM_JOBID`; do
                        job=${line%% *}
                        deps=`echo $line | awk '{print $2}'`

                        if [[ $deps == null ]]; then
                            deps=""
                        elif [[ $deps == ${deps/\./} ]]; then
                            deps="Dependency=afterok:$newJobID"
                        else
                            tmp=""; IFS=$' '
                            for t in ${deps//\./ }; do
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
        #set +x;
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
        for try in {1..5}; do
            if [ ! -f $failFile.requeued.$try.time ]; then
                touch $failFile.requeued.$try.time
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
                    min
                    newFactor=2
                elif [ $min -lt 480 ]; then
                    min=480
                    newFactor=1.5
                else
                    newFactor=1.2
                fi

                #newFactor=2
                min=`echo "($min*$newFactor)/1" | bc`
                echo trying to requeue $try with $min minutes
                echo $totalM $min $maxExtra > ${out%.out}.adjust


                hours=$((($min + 59) / 60))
                adjustPartition $hours $partition

                seconds=$(($min * 60))

                # https://chat.openai.com/chat/80e28ff1-4885-4fe3-8f21-3556d221d7c6

                time=`eval "echo $(date -ud "@$seconds" +'$((%s/3600/24))-%H:%M:%S')"`
                export  myPartition=$partition
                export myTime=$time
                export myMem=${totalM}M
                requeueCmd=`grep "Command used to submit the job:" $script | tail -n 1`
                requeueCmd=${requeueCmd#*submit the job: }
                requeueCmd=${requeueCmd//\$myPartition/$myPartition}
                requeueCmd=${requeueCmd//\$myTime/$myTime}
                requeueCmd=${requeueCmd//\$myMem/$myMem}
                newJobID=`$requeueCmd`

                if [[ "$newJobID" =~ ^[0-9]+$ ]]; then
                    echo "# mem=$myMem time=$myTime " >> $script
                    IFS=$'\n'
                    for line in `grep $SLURM_JOBID $smartSlurmLogDir/allJobs.txt | grep -v ^$SLURM_JOBID`; do
                        job=${line%% *}
                        deps=`echo $line | awk '{print $2}'`

                        if [[ $deps == null ]]; then
                            deps=""
                        elif [[ $deps == ${deps/\./} ]]; then
                            deps="Dependency=afterok:$newJobID"
                        else
                            tmp=""; IFS=$' '
                            for t in ${deps//\./ }; do
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
        
                
                    
                
                # sleep 1
                # echo trying to requeue $try
                # touch $failFile.requeued.$try.time

                # # 80G memory
                # #[ "${mem%M}" -gt 81920 ] && [ "$try" -gt 2 ] && break

                # scontrol requeue $SLURM_JOBID && echo job re-submitted || echo job not re-submitted.

                # # time=${10}
                # # [[ "$time" == *-* ]] && { day=${time%-*}; tem=${time#*-}; hour=${tem%%:*}; min=${tem#*:}; min=${min%%:*}; sec=${tem#$hour:$min}; sec=${sec#:}; } || { [[ "$time" =~ ^[0-9]+$ ]] && min=$time || { sec=${time##*:}; min=${time%:*}; min=${min##*:}; hour=${time%$min:$sec}; hour=${hour%:}; day=0;} }

                # # [ -z "$day" ] && day=0; [ -z "$hour" ] && hour=0; [ -z "$min" ] && min=0;[ -z "$sec" ] && sec=0

                # # echo $day $day,  $hour hour,  $min min,  $sec sec

                # factor=2 #$((1 + 1/e(0.1 * $try)))

                # # # how many hours for sbatch command if we double the time
                # # hours=$(($day * 2 * 24 + $hour * 2 + ($min * 2 + 59 + ($sec * 2 + 59) / 60 ) / 60))

                # #min=${10} # the orignal estimated time

                # #[ "$min" -lt 20 ] && min=20 # at least 20 minutes

                # hours=$((($min * $factor + 59) / 60))

                # adjustPartition $hours $partition

                
                # if [[ "$partition" != "${partition}" ]]; then
                #     scontrol update jobid=$SLURM_JOBID Partition=$partition TimeLimit=$time
                # else
                #     scontrol update jobid=$SLURM_JOBID TimeLimit=$time
                #     #scontrol update jobstep=123456.2 TimeLimit=02:00:00
                # fi

                # echo $totalM $(( min * factor )) $extraMemC > ${out%.out}.adjust
                # echo job resubmitted: $SLURM_JOBID with time: $time partition: $partition, mem is not changed

                # [ -f $failFile ] && rm $failFile

                # echo 0 0 0 0 0 0 0 >> $smartSlurmLogDir/job_$SLURM_JOBID.memCPU.txt

                # break
            fi
        done
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
    fi
    if [ "$min" -ge 20 ] && [ ! -z "$alwaysRequeueIfFail" ] && [ "$jobStatus" != "Cancelled" ]; then
        ( sleep 2;  scontrol requeue $SLURM_JOBID; ) &
    fi

fi

echo "Sending email..."

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

summarizeRun.sh $smartSlurmLogDir

toSend="`cat $smartSlurmLogDir/summary`\n$toSend"

s="${toSend%% *} $s"

#echo -e "tosend:\n$toSend"
#echo -e "$toSend" >> ${err%.err}.email

#echo -e "$toSend" | sendmail `head -n 1 ~/.forward`
if [ -f $smartSlurmJobRecordDir/stats/$software.$ref.mem.png ]; then
    echo -e "$toSend" | mail -s "$s" -a $smartSlurmLogDir/job_$SLURM_JOBID.mem.png -a $smartSlurmLogDir/job_$SLURM_JOBID.cpu.png -a $smartSlurmLogDir/barchartMem.png -a $smartSlurmLogDir/barchartTime.png -a $smartSlurmJobRecordDir/stats/$software.$ref.mem.png -a $smartSlurmJobRecordDir/stats/$software.$ref.time.png $USER && echo email sent || \
        { echo Email not sent.; echo -e "$s\n$toSend" | sendmail `head -n 1 ~/.forward` && echo Email sent by second try. || echo Email still not sent!!; }
elif [ -f $smartSlurmJobRecordDir/stats/back/$software.$ref.time.png ]; then
    echo -e "$toSend" | mail -s "$s" -a $smartSlurmLogDir/job_$SLURM_JOBID.mem.png -a $smartSlurmLogDir/job_$SLURM_JOBID.cpu.png -a $smartSlurmLogDir/barchartMem.png -a $smartSlurmLogDir/barchartTime.png -a $smartSlurmJobRecordDir/stats/back/$software.$ref.mem.png -a $smartSlurmJobRecordDir/stats/back/$software.$ref.time.png $USER && echo email sent || \
        { echo Email not sent.; echo -e "$s\n$toSend" | sendmail `head -n 1 ~/.forward` && echo Email sent by second try. || echo Email still not sent!!; }

else
    echo -e "$toSend" | mail -s "$s" -a $smartSlurmLogDir/job_$SLURM_JOBID.mem.png -a $smartSlurmLogDir/job_$SLURM_JOBID.cpu.png -a $smartSlurmLogDir/barchartMem.png -a $smartSlurmLogDir/barchartTime.png $USER && echo email sent || \
    { echo Email not sent.; echo -e "$s\n$toSend" | sendmail `head -n 1 ~/.forward` && echo Email sent by second try. || echo Email still not sent!!; }

fi

if [[ $USER != ld32 ]]; then
    if [ -f $smartSlurmJobRecordDir/stats/$software.$ref.mem.png ]; then
        echo -e "$toSend" | mail -s "$s" -a $smartSlurmLogDir/job_$SLURM_JOBID.mem.png -a $smartSlurmLogDir/job_$SLURM_JOBID.cpu.png -a $smartSlurmLogDir/barchartMem.png -a $smartSlurmLogDir/barchartTime.png -a $smartSlurmJobRecordDir/stats/$software.$ref.mem.png -a $smartSlurmJobRecordDir/stats/$software.$ref.time.png ld32
    elif [ -f $smartSlurmJobRecordDir/stats/back/$software.$ref.time.png ]; then
        echo -e "$toSend" | mail -s "$s" -a $smartSlurmLogDir/job_$SLURM_JOBID.mem.png -a $smartSlurmLogDir/job_$SLURM_JOBID.cpu.png -a $smartSlurmLogDir/barchartMem.png -a $smartSlurmLogDir/barchartTime.png -a $smartSlurmJobRecordDir/stats/back/$software.$ref.mem.png -a $smartSlurmJobRecordDir/stats/back/$software.$ref.time.png ld32

    else
        echo -e "$toSend" | mail -s "$s" -a $smartSlurmLogDir/job_$SLURM_JOBID.mem.png -a $smartSlurmLogDir/job_$SLURM_JOBID.cpu.png -a $smartSlurmLogDir/barchartMem.png -a $smartSlurmLogDir/barchartTime.png ld32
    fi
fi

echo

# create an empty file so that it is easier to match job name to job ID
#touch $out.$SLURM_JOB_ID
# move this to job sbumission time,

# wait for email to be sent
sleep 5
