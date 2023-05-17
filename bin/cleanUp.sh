#!/bin/bash

#set -x 

# to call this:  0     1           2           3       4         5          6       7        8     9     10      11       12           13 
#cleanUp.sh       "projectDir"  "$software" "$ref" "$flag" "$inputSize"   $core   $memO  $timeO   $mem  $time  $partition slurmAcc  inputs 


echo Running $0 $@

if [ -f ~/.smartSlurm/config/config.txt ]; then 
    source ~/.smartSlurm/config/config.txt
else     
    source $(dirname $0)/../config/config.txt || { echo Config list file not found: config.txt; exit 1; }
fi

#touch /tmp/job_$SLURM_JOB_ID.done
if [[ -z "$1" ]]; then 

   #out=$4.out; out=${out/\%jerr=${4##* }; err=${err/\%j/$SLURM_JOB_ID}; script=${4% *}; script=${script#* }; succFile=${script/\.sh/}.success;      failFile=${script/\.sh/}.failed; 
    out=slurm-$SLURM_JOBID.out; err=slurm-$SLURM_JOBID.err; script=$4.sh; succFile=$4.success; failFile=$4.failed;
else 
    out=$1/log/"${4}.out"; err=$1/log/${4}.err; script=$1/log/${4}.sh; succFile=$1/log/${4}.success; failFile=$1/log/${4}.failed; checkpointDir=$1/log/${4} 
fi 

# wait for slurm database update
sleep 10

sacct=`sacct --format=JobID,Submit,Start,End,MaxRSS,State,NodeList%30,Partition,ReqTRES%30,TotalCPU,Elapsed%14,Timelimit%14 --units=M -j $SLURM_JOBID` 

#sacct=`cat ~/fakeSacct.txt`

#sacct report is not accurate, let us use the total memory
#totalM=${sacct#*,mem=}; totalM=${totalM%%M,n*}

echo -e "\nJob summary:\n$sacct"
# echo *Notice the sacct report above: while the main job is still running for sacct command, user task is completed.

# record job for future estimating mem and time
jobStat=`echo -e "$sacct" | tail -n 1`

#from: "sacct --format=JobID,Submit,Start,End,MaxRSS,State,NodeList%30,Partition,ReqTRES%30,TotalCPU,Elapsed%14,Timelimit%14 --units=M -j $SLURM_JOBID" 

START=`echo $jobStat | cut -d" " -f3`

START=`date -d "$START" +%s`

FINISH=`echo $jobStat | cut -d" " -f4`      

FINISH=`date -d "$FINISH" +%s`

[ -z "$FINISH" ] && FINISH=`date +%s`

echo  start: $START fisnish: $FINISH

# time in minutes
min=$((($FINISH - $START + 59)/60))

#[[ "$mim" == 0 ]] && min=1

# memory in M
memSacct=`echo $jobStat | cut -d" " -f5`

[[ "$memSacct" == "RUNNING" ]] && memSacct=NA

# node
node=`echo $jobStat | cut -d" " -f7`

case "$jobStat" in
# jobRecord.txt header
#1user 2software 3ref 4inputName 5inputSizeInK 6CPUNumber 7memoryO 8timeO 9readMem 10RequestedTime 11jobID 12memoryM 13minRun 14Node 15 finalStatus    
*COMPLETED* )  jobStatus="COMPLETED" && echo *Notice the sacct report above: while the main job is still running for sacct command, user task is completed.;;

*TIMEOUT*   )  jobStatus="OOT";;

*OUT_OF_ME*   ) jobStatus="OOM";;
        
*CANCELLED*	) jobStatus="Cancelled";;

*FAILED*	) jobStatus="Fail";;

*          )  jobStatus="Unknown";;


esac

# for testing
#jobStatus="OOM"

echo jobStatus: $jobStatus 

# sacct might give wrong resules
[[ $jobStatus != "OOM" ]] && [[ $jobStatus != "OOT" ]] && [[ $jobStatus != "Cancelled" ]] && [ -f $succFile ] && jobStatus="COMPLETED"

echo -e  "Last row of job summary: $jobStat" 
echo start: $START finish: $FINISH mem: $memSacct min: $min



# sacct actually give very not accurate result. Let's use cgroup report
#mem=`cat /sys/fs/cgroup/memory/slurm/uid_*/job_$SLURM_JOBID/memory.max_usage_in_bytes`

srunM=`sort -n log/job_$SLURM_JOBID.mem.txt | tail -n1 | cut -d' ' -f2`
#srunM=$((srunM / 1024 / 1024 ))

echo jobStatus: $jobStatus cgroupMaxMem: $srunM 


memSacct=${memSacct%M}; memSacct=${memSacct%.*} #remove M and decimals

 
# Not sure if this is needed.
[[ "$memSacct" != "NA" ]] && [ "$memSacct" -gt "$srunM" ] && srunM=$memSacct

# sometimes, the last row say srun cancelled, but the job is actually out of memory or out
if [[ $jobStatus == "Cancelled" ]]; then  
    if [[ "$sacct" == *TIMEOUT* ]]; then 
        echo The job is actually timeout
        jobStatus="OOT"
    elif [[ "$sacct" == *OUT_OF_ME* ]]; then 
        echo The job is actually out-of-memory
        jobStatus="OOM"
    fi
fi
#set -x

#some times, it reports unknown, but is it is oom
if [[ $jobStatus == Unknown ]]; then
    tLog=`tail -n 22 $out | grep ^srun`
    [[ "$tLog" == *"task 0: Out Of Memory"* ]] && jobStatus="OOM" && echo The job is actually out-of-memory by according to the log: && echo $tLog
fi  
set +x 

inputs=${13}
inputSize=$5     # input does not exist when job was submitted.
if [ "$inputs" != "none" ] && [ "$5" == "0" ]; then
    inputSize=`{ du --apparent-size -c -L ${inputs//,/ } 2>/dev/null || echo notExist; } | tail -n 1 | cut -f 1`

    if [[ "$inputSize" == "notExist" ]]; then 
        echo Some or all input files not exist: $inputs
        echo Error! missingInputFile: ${inputs//,/ } 
        exit 
    fi
fi    

[ -f ${out%.out}.adjust ] && totalM=`cat ${out%.out}.adjust | cut -d' ' -f1` || totalM=$9

[ -f ${out%.out}.adjust ] && totalT=`cat ${out%.out}.adjust | cut -d' ' -f2` || totalT=${10}

[ -f ${out%.out}.adjust ] && extraMemC=`cat ${out%.out}.adjust | cut -d' ' -f3` || extraMemC=${14}

[ -f $jobRecordDir/stats/extraMem.$2.$3 ] && extraMem=`sort -nr $jobRecordDir/stats/extraMem.$2.$3 | head -n1`

# email alert if only half of resource is used
# totalM=200; totalT=200; srunM=1000; min=200; 		
if [ "$7" -ne "$totalM" ] && [ "$8" -ne "$totalT" ]; then 
   if [ $((totalM-extraM)) -gt $((srunM * 2)) ] || [ $totalT -gt $(( min * 2 )) ]; then 
	echo Over-reserved! 
	echo -e "`pwd`\n$out\n$SLURM_JOBID over-reserved resounce M: $srunM/$totalM T: $min/$totalT" | mail -s "Over reserved $SLURM_JOBID" $USER   
   fi
fi    
                                #3defult,  5given,  7cGroupUsed                  sacct used      
 record="$SLURM_JOB_ID,$inputSize,$7,$8,$totalM,$totalT,$srunM,$min,$jobStatus,$USER,$memSacct,$2,$3,$4,$6,$extraMemC,$extraTime,`date`"  # 16 extraM 
echo dataToPlot,$record
    
#if [[ ! -f $jobRecordDir/stats/$2.$3.mem.stat || "$2" == "regularSbatch" ]]; then 
    
if [[ $jobStatus == "COMPLETED" ]]; then 
   
    #if [[ "$inputSize" == 0 ]]; then # || "$2" == "regularSbatch" ]] ; then
        records=`grep COMPLETED $jobRecordDir/jobRecord.txt 2>/dev/null | awk -F"," -v a=$2 -v b=$3 '{ if($12 == a && $13 == b) {print $2, $7 }}' | sort -u -n` 
        #timeRecords=`grep COMPLETED $jobRecordDir/jobRecord.txt 2>/dev/null | awk -F"," -v a=$2 -v b=$3 '{ if($12 == a && $13 == b) {print $8 }}' | sort -u -nr | tr '\n' ' '`
        #memRecords=`grep COMPLETED $jobRecordDir/jobRecord.txt 2>/dev/null | awk -F"," -v a=$2 -v b=$3 '{ if($12 == a && $13 == b) {print $7 }}' | sort -u -nr | tr '\n' ' '` 
        #timeRecords=`grep COMPLETED $jobRecordDir/jobRecord.txt 2>/dev/null | awk -F"," -v a=$2 -v b=$3 '{ if($12 == a && $13 == b) {print $8 }}' | sort -u -nr | tr '\n' ' '`
    
    
        #maxMem=`cat $jobRecordDir/stats/$2.$3.mem.stat.noInput  2>/dev/null | sort -nr | tr '\n' ' ' | cut -f 1 -d " "`
        #maxTime=`cat $jobRecordDir/stats/$2.$3.time.stat.noInput  2>/dev/null | sort -nr | tr '\n' ' ' | cut -f 1 -d " "
        #if [ "$(echo $memRecords | wc -w)" -lt 20 ] || [ "${memRecords%% *}" -lt "${srunM%.*}" ] || [ "${timeRecords%% *}" -lt "$min" ]; then
        
        #   less than 20 records                  # or current one is larger than all old data   # and not exist already
        if [ "$(echo $records | wc -l)" -lt 200 ] || [ "`echo $records | tail -n1 | cut -d' ' -f1`" -lt "$inputSize" ] && ! echo $records | grep "$inputSize $srunM"; then
            echo $record >> $jobRecordDir/jobRecord.txt
            echo -e "Added this line to $jobRecordDir/jobRecord.txt:\n$record"
            mv $jobRecordDir/stats/$2.$3.* $jobRecordDir/stats/back 2>/dev/null
        else 
            echo Did not add this record to $jobRecordDir/jobRecord.txt
        fi  
    # else         
    #     maxMem=`cat $jobRecordDir/stats/$2.$3.mem.txt  2>/dev/null | cut -f 2 -d " " | sort -nr | tr '\n' ' ' | cut -f 1 -d ' '`
        
    #     maxTime=`cat $jobRecordDir/stats/$2.$3.time.txt  2>/dev/null | cut -f 2 -d " " | sort -nr | tr '\n' ' ' | cut -f 1 -d ' '`

    #     if [ -z "$maxMem" ] || [ "${maxMem%.*}" -lt "${srunM%.*}" ] || [ -z "$maxTime" ] || [ "$maxTime" -lt "$min" ]; then
    #         echo $record >> $jobRecordDir/jobRecord.txt
    #         echo -e "Added this line to $jobRecordDir/jobRecord.txt:\n$record"
    #         rm $jobRecordDir/stats/$2.$3.*  2>/dev/null
    #     else 
    #         echo Did not add this record to $jobRecordDir/jobRecord.txt
    #     fi  
    # fi
#else
    # todo: may not need failed job records?
#    echo $record >> $jobRecordDir/jobRecord.txt
    #echo -e "Added this line to $jobRecordDir/jobRecord.txt:\n$record"
fi
echo Final mem: $srunM, time: $min minutes
   

# for testing
#rm $succFile; jobStatus=OOM 

if [ ! -f $succFile ]; then
    touch $failFile

    # if checkpoint failed or out of memory
    if  [ -f ${out%.out}.likelyCheckpointOOM ] || [[ "$jobStatus" == "OOM" ]]; then
        
        mv ${out%.out}.likelyCheckpointOOM ${out%.out}.likelyCheckpointOOM.old

        jobStatus=OOM

        #set -x
        echo Will re-queue after sending email...
    
        echo old extraMem: 
        cat $jobRecordDir/stats/extraMem.$2.$3

        #extraMemN=$(( ( totalM - srunM ) *2 ))
        extraMemN=$(( totalM - srunM + 1 ))
        #[[ "$extraMemN" == 0 ]] && extraMemN=1
        [ $extraMemN -gt 0 ] && echo $extraMemN $totalM $inputSize $SLURM_JOBID >> $jobRecordDir/stats/extraMem.$2.$3
        #oomCount=`wc -l $jobRecordDir/stats/extraMem.$2.$3 | cut -d' ' -f1`
        maxExtra=`sort -n $jobRecordDir/stats/extraMem.$2.$3 | tail -n1 | cut -d' ' -f1`
        [ -z "$maxExtra" ] && maxExtra=5
        extraMem=$(( $maxExtra * 2 ))
        
        echo new extraMem: 
        cat $jobRecordDir/stats/extraMem.$2.$3
        #set +x        
        ( sleep 5; 
        for try in {1..8}; do
            if [ ! -f $failFile.requeued.$try.mem ]; then
                sleep 2
                #mem=$(( $mem * ( 2 ^ $try ) ))
                touch $failFile.requeued.$try.mem
                #alpha=1
                #newFactor=`echo "1.2+1/e($alpha*$mem/1000)" | bc -l | xargs printf "%.2f"`
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
                mem=`echo "($totalM*$newFactor+$extraMem)/1" | bc`
                # try=1, then factor is 5, try 2 factor is 3, try 3 is 2 ...
                #mem=$(( $mem * (1 + 1/e(0.1 * $try))))
                #mem=4000
                echo trying to requeue $try with $mem M
                echo $mem $totalT $extraMem > ${out%.out}.adjust
                
                # 80G memory
                #[ "$mem" -gt 81920 ] && [ "$try" -gt 2 ] && break

                hostName=`hostname`
                
                if `ssh $hostName "scontrol requeue $SLURM_JOBID; scontrol update JobId=$SLURM_JOBID MinMemoryNode=$mem;"`; then

                #if `scontrol requeue $SLURM_JOBID; scontrol update JobId=$SLURM_JOBID MinMemoryNode=$mem`; then 
                    echo Requeued successfully
                    [ -f $failFile ] && rm $failFile

                    s="Requeued:$SLURM_JOBID:$SLURM_JOB_NAME"
                    echo -e "" | mail -s "$s" $USER

                    break
                fi    
        
            fi    
        done ) &
        #else
           # echo Could not find the original mem value.
        #    echo Job failed of out-of-memory. Please resubmit with more memory check youself.
        #fi  
        
        # delete stats and redo them
        if [ "$inputs" == "none" ]; then
            mv $jobRecordDir/stats/$2.$3.* $jobRecordDir/stats/back  2>/dev/null
        else 
            # remove bad records
            if [ -f $jobRecordDir/stats/$2.$3.mem.stat ]; then 
                # .  $jobRecordDir/stats/$2.$3.mem.stat 
                # Finala=`printf "%.15f\n" $Finala`
                # Finalb=`printf "%.15f\n" $Finalb`
                # Maximum=`printf "%.15f\n" $Maximum`
                # echo Finala: $Finala Finalb: $Finalb Maximum: $Maximum

                # awk -F"," -v a=$2 -v b=$3 -v c=$Finala -v d=$Finalb '{ if ( ! ($12 == a && $13 == b && $2 * c + d > $7)) {print}}' $jobRecordDir/jobRecord.txt > $jobRecordDir/jobRecord.txt1 
                # echo diff output:
                # diff $jobRecordDir/jobRecord.txt $jobRecordDir/jobRecord.txt1
                # mv $jobRecordDir/jobRecord.txt1 $jobRecordDir/jobRecord.txt
                mv $jobRecordDir/stats/$2.$3.* $jobRecordDir/stats/back  2>/dev/null
            fi
        fi 

        # awk -F"," -v a=$s -v b=$r -v c=$Finala -v d=$Finalb '{ ! ($12 == a && $13 == b && $2 * c + d > $7) {print}}' /home/ld32/.smartSlurm/jobRecord.txt

        # awk -F"," -v a=$s -v b=$r -v c=$Finala -v d=$Finalb '{ ! ($12 == a && $13 == b && $2 * c + d > $7) {print}}' $jobRecordDir/jobRecord.txt > $jobRecordDir/jobRecord.txt1 


        #rm $jobRecordDir/stats/$2.$3.mem.stat* $jobRecordDir/stats/$2.$3.time.stat* 2>/dev/null
        #scontrol requeue $SLURM_JOBID && echo job re-submitted || echo job not re-submitted.
        #scontrol requeue 45937 && echo job re-submitted || echo job not re-submitted.

        # sleep 5 45937
        #scontrol hold $SLURM_JOBID
        #scontrol hold 45937
        #sleep 5
        # this doe not work due to cgroup out memory error 
        #sbatch -p priority -t 5 -A rccg --wrap "scontrol requeue $SLURM_JOBID; scontrol update JobId=$SLURM_JOBID MinMemoryCPU=50000;" 
        # sleep 5; scontrol release $SLURM_JOBID; "
        #which sbatch 

        #sbatch -p priority -t 5 -A rccg --wrap "hostname"

        #scontrol update JobId=$SLURM_JOBID MinMemoryCPU=8000
        # todo: submit new job and adjust down steam jobs' dependency

        #echo job resubmitted: $SLURM_JOBID with mem: $mem
    
    elif [[ "$jobStatus" == "OOT" ]]; then
        
        for try in {1..8}; do
            if [ ! -f $failFile.requeued.$try.time ]; then
                sleep 2 
                echo trying to requeue $try
                touch $failFile.requeued.$try.time 

                # 80G memory
                #[ "${mem%M}" -gt 81920 ] && [ "$try" -gt 2 ] && break

             
                scontrol requeue $SLURM_JOBID && echo job re-submitted || echo job not re-submitted.

                # time=${10}
                # [[ "$time" == *-* ]] && { day=${time%-*}; tem=${time#*-}; hour=${tem%%:*}; min=${tem#*:}; min=${min%%:*}; sec=${tem#$hour:$min}; sec=${sec#:}; } || { [[ "$time" =~ ^[0-9]+$ ]] && min=$time || { sec=${time##*:}; min=${time%:*}; min=${min##*:}; hour=${time%$min:$sec}; hour=${hour%:}; day=0;} }

                # [ -z "$day" ] && day=0; [ -z "$hour" ] && hour=0; [ -z "$min" ] && min=0;[ -z "$sec" ] && sec=0

                # echo $day $day,  $hour hour,  $min min,  $sec sec

                factor=2 #$((1 + 1/e(0.1 * $try)))

                # # how many hours for sbatch command if we double the time
                # hours=$(($day * 2 * 24 + $hour * 2 + ($min * 2 + 59 + ($sec * 2 + 59) / 60 ) / 60))

                #min=${10} # the orignal estimated time

                #[ "$min" -lt 20 ] && min=20 # at least 20 minutes

                hours=$((($min * $factor + 59) / 60))

                adjustPartition $hours $partition

                seconds=$(($min * $factor  * 60))

                #[ $seconds -le 60 ] && time=11:0 && seconds=60

                #echo srun seconds: $seconds 
                
                # https://chat.openai.com/chat/80e28ff1-4885-4fe3-8f21-3556d221d7c6

                time=`eval "echo $(date -ud "@$seconds" +'$((%s/3600/24))-%H:%M:%S')"`

                if [[ "$partition" != "${11}" ]]; then 
                    scontrol update jobid=$SLURM_JOBID Partition=$partition TimeLimit=$time
                else 
                    scontrol update jobid=$SLURM_JOBID TimeLimit=$time
                    #scontrol update jobstep=123456.2 TimeLimit=02:00:00

                fi 
                
                echo $totalM $(( min * factor )) $extraMemC > ${out%.out}.adjust
                echo job resubmitted: $SLURM_JOBID with time: $time partition: $partition, mem is not changed

                [ -f $failFile ] && rm $failFile


                s="Requeued:$SLURM_JOBID:$SLURM_JOB_NAME"
                echo -e "" | mail -s "$s" $USER

                break
            fi
        done 
         # delete stats and redo them
        if [[ "$inputs" == "none" ]]; then
            mv $jobRecordDir/stats/$2.$3.* $jobRecordDir/stats/back  2>/dev/null

            #[ -f $jobRecordDir/stats/$2.$3.mem.stat.noInput ] && mv $jobRecordDir/stats/$2.$3.mem.stat.noInput $jobRecordDir/stats/$2.$3.mem.stat.noInput.$(stat -c %y $jobRecordDir/stats/$2.$3.mem.stat.noInput | tr " " ".")
        else 
             # remove bad records
            if [ -f $jobRecordDir/stats/$2.$3.time.stat ]; then 
                # .  $jobRecordDir/stats/$2.$3.time.stat
                # Finala=`printf "%.15f\n" $Finala`
                # Finalb=`printf "%.15f\n" $Finalb`
                # Maximum=`printf "%.15f\n" $Maximum`
                # echo Finala: $Finala Finalb: $Finalb Maximum: $Maximum

                # awk -F"," -v a=$2 -v b=$3 -v c=$Finala -v d=$Finalb '{ if ( ! ($12 == a && $13 == b && $2 * c + d > $8 ) ) {print}}' $jobRecordDir/jobRecord.txt > $jobRecordDir/jobRecord.txt1 
                # echo diff output: 
                # diff $jobRecordDir/jobRecord.txt $jobRecordDir/jobRecord.txt1
                # mv $jobRecordDir/jobRecord.txt1 $jobRecordDir/jobRecord.txt
                mv $jobRecordDir/stats/$2.$3.* $jobRecordDir/stats/back  2>/dev/null
            fi
        fi 
    
    else 
        echo Not sure why job failed. Not run out of time or memory. Pelase check youself.
    fi
elif [ ! -z "$1" ]; then
    adjustDownStreamJobs.sh $1/log
    rm $failFile 2>/dev/null
else 
    rm $failFile 2>/dev/null
fi

echo "Sending email..."

minimumsize=9000

actualsize=`wc -c $out || echo 0`

[ -f $succFile ] && s="Succ:$SLURM_JOBID:$SLURM_JOB_NAME" || s="$jobStatus:$SLURM_JOBID:$SLURM_JOB_NAME" 

if [ "${actualsize% *}" -ge "$minimumsize" ]; then
   #toSend=`echo Job script content:; cat $script;`
   toSend="Log dir: `dirname $out`\n"
   toSend="$toSend\nOutput is too big for email. Please find output in: $out"  
   toSend="$toSend\n...\nFirst 20 row of output:\n`head -n 20 $out`"
   toSend="$toSend\n...\nLast 20 row of output:\n`tail -n 20 $out`"
else
   #toSend=`echo Job script content:; cat $script; echo; echo Job output:; cat $out;`
   toSend="Log: $out"
   toSend=`echo; echo Job log:; cat $out;`
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



#cp /tmp/job_$SLURM_JOBID.mem.txt log/

summarizeRun.sh log/allJobs.txt $flag





#echo -e "tosend:\n$toSend"
echo -e "$toSend" >> ${err%.err}.email

#echo -e "$toSend" | sendmail `head -n 1 ~/.forward`
if [ -f $jobRecordDir/stats/$2.$3.mem.png ]; then 
    echo -e "$toSend" | mail -s "$s" -a log/job_$SLURM_JOBID.mem.png -a log/barchartMem.png -a log/barchartTime.png -a $jobRecordDir/stats/$2.$3.mem.png -a $jobRecordDir/stats/$2.$3.time.png $USER && echo email sent || \
        { echo Email not sent.; echo -e "Subject: $s\n$toSend" | sendmail `head -n 1 ~/.forward` && echo Email sent by second try. || echo Email still not sent!!; }
elif [ -f $jobRecordDir/stats/back/$2.$3.time.png ]; then
    echo -e "$toSend" | mail -s "$s" -a log/job_$SLURM_JOBID.mem.png -a log/barchartMem.png -a log/barchartTime.png -a $jobRecordDir/stats/back/$2.$3.mem.png -a $jobRecordDir/stats/back/$2.$3.time.png $USER && echo email sent || \
        { echo Email not sent.; echo -e "Subject: $s\n$toSend" | sendmail `head -n 1 ~/.forward` && echo Email sent by second try. || echo Email still not sent!!; }

else 
    echo -e "$toSend" | mail -s "$s" -a log/job_$SLURM_JOBID.mem.png -a log/barchartMem.png -a log/barchartTime.png $USER && echo email sent || \
    { echo Email not sent.; echo -e "Subject: $s\n$toSend" | sendmail `head -n 1 ~/.forward` && echo Email sent by second try. || echo Email still not sent!!; }        

fi

if [[ $USER != ld32 ]]; then
    if [ -f $jobRecordDir/stats/$2.$3.mem.png ]; then 
        echo -e "$toSend" | mail -s "$s" -a log/job_$SLURM_JOBID.mem.png -a log/barchartMem.png -a log/barchartTime.png -a $jobRecordDir/stats/$2.$3.mem.png -a $jobRecordDir/stats/$2.$3.time.png ld32 
    elif [ -f $jobRecordDir/stats/back/$2.$3.time.png ]; then
        echo -e "$toSend" | mail -s "$s" -a log/job_$SLURM_JOBID.mem.png -a log/barchartMem.png -a log/barchartTime.png -a $jobRecordDir/stats/back/$2.$3.mem.png -a $jobRecordDir/stats/back/$2.$3.time.png ld32 

    else 
        echo -e "$toSend" | mail -s "$s" -a log/job_$SLURM_JOBID.mem.png -a log/barchartMem.png -a log/barchartTime.png ld32
    fi
fi

echo 

# create an empty file so that it is easier to match job name to job ID
#touch $out.$SLURM_JOB_ID
# move this to job sbumission time, so that we have that file as early as possible

# wait for email to be sent
sleep 10
