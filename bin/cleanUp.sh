#!/bin/bash

#set -x 

# to call this:  0     1           2           3       4         5          6       7        8       9    10      11       12           13 
#cleanUp.sh       "projectDir"  "$software" "$ref" "$flag" "$inputSize"   $core   $memO  $timeO    $mem  $time  $partition slurmAcc  original.sbatch.command

echo Running $0 $@

if [[ -z "$1" ]]; then 

   #out=$4.out; out=${out/\%jerr=${4##* }; err=${err/\%j/$SLURM_JOB_ID}; script=${4% *}; script=${script#* }; succFile=${script/\.sh/}.success;      failFile=${script/\.sh/}.failed; 
    out=slurm-$SLURM_JOBID.out; err=slurm-$SLURM_JOBID.err; script=$4.sh; succFile=$4.success; failFile=$4.failed;
else 
    out=$1/logs/"${4}.out"; err=$1/logs/${4}.err; script=$1/logs/${4}.sh; succFile=$1/logs/${4}.success; failFile=$1/logs/${4}.failed;   
fi 

sleep 5

sacct=`sacct --format=JobID,Submit,Start,End,MaxRSS,State,NodeList%30,Partition,ReqTRES%30,TotalCPU,Elapsed%14,Timelimit%14 --units=M -j $SLURM_JOBID` 

#sacct=`cat ~/fakeSacct.txt`

echo -e "\nJob summary:\n$sacct"
# echo *Notice the sacct report above: while the main job is still running for sacct command, user task is completed.

# record job for future estimating mem and time
jobStat=`echo -e "$sacct" | tail -n 1`

#from: "sacct --format=JobID,Submit,Start,End,MaxRSS,State,NodeList%30,Partition,ReqTRES%30,TotalCPU,Elapsed%14,Timelimit%14 --units=M -j $SLURM_JOBID" 

START=`echo $jobStat | cut -d" " -f3`

START=`date -d "$START" +%s`

FINISH=`echo $jobStat | cut -d" " -f4`      

FINISH=`date -d "$FINISH" +%s`

# time in minutes
mins=$((($FINISH - $START + 59)/60))

# memory in M
mem=`echo $jobStat | cut -d" " -f5`

# node
node=`echo $jobStat | cut -d" " -f7`

case "$jobStat" in
# jobRecord.txt header
#1user 2software 3ref 4inputName 5inputSizeInK 6CPUNumber 7memoryO 8timeO 9readMem 10RequestedTime 11jobID 12memoryM 13minRun 14Node 15 finalStatus    
*COMPLETED* )  jobStatus="COMPLETED" && echo *Notice the sacct report above: while the main job is still running for sacct command, user task is completed.;;

*TIMEOUT*   )  jobStatus="OOT";;

*OUT_OF_ME*   ) jobStatus="OOM";;
        
*CANCELLED*	) jobStatus="Cancelled";;

*          )  jobStatus="Unknown";;


esac

# for testing
#jobStatus="OOM"

[[ $jobStatus != "COMPLETED" ]] && [ -f $succFile ] && rm $succFile

echo -e  "Last row of job summary: $jobStat" 
echo start: $START finish: $FINISH mem: $mem mins: $mins
echo jobStatus: $jobStatus

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

record="$SLURM_JOB_ID,$5,$7,$8,$9,${10},${mem%M},${mins},$jobStatus,$USER,$1,$2,$3,$4,$6,${node},$err,`date`,\"${13}\""

    
if [ ! -z "$record" ]; then
    #if [[ ! -f ~/smartSlurm/stats/$2.$3.mem.stat || "$2" == "regularSbatch" ]]; then 
        
    if [[ $jobStatus == "COMPLETED" ]]; then 
        memm=${mem%M*}
        if [ "$5" == 0 ]; then # || "$2" == "regularSbatch" ]] ; then
            maxMem=`cat ~/smartSlurm/stats/$2.$3.mem.stat.noInput | sort -nr | tr '\n' ' ' | cut -f 1 -d " "`
            maxTime=`cat ~/smartSlurm/stats/$2.$3.time.stat.noInput | sort -nr | tr '\n' ' ' | cut -f 1 -d " "`
            
            if [ -z "$maxMem" ] || [ "${maxMem%.*}" -lt "${memm%.*}" ] || [ -z "$maxTime" ] || [ "$maxTime" -lt "$mins" ]; then
                echo $record >> ~/smartSlurm/jobRecord.txt
                echo -e "Added this line to ~/smartSlurm/jobRecord.txt:\n$record"
                rm ~/smartSlurm/stats/$2.$3.mem.stat.noInput ~/smartSlurm/stats/$2.$3.time.stat.noInput
            else 
                echo Did not add this record to ~/smartSlurm/stats/jobRecord.txt
            fi  
        else         
            maxMem=`cat ~/smartSlurm/stats/$2.$3.mem.txt | cut -f 2 -d " " | sort -nr | tr '\n' ' ' | cut -f 1 -d ' '`
            
            maxTime=`cat ~/smartSlurm/stats/$2.$3.time.txt | cut -f 2 -d " " | sort -nr | tr '\n' ' ' | cut -f 1 -d ' '`
    
            if [ -z "$maxMem" ] || [ "${maxMem%.*}" -lt "${memm%.*}" ] || [ -z "$maxTime" ] || [ "$maxTime" -lt "$mins" ]; then
                echo $record >> ~/smartSlurm/jobRecord.txt
                echo -e "Added this line to ~/smartSlurm/jobRecord.txt:\n$record"
                rm ~/smartSlurm/stats/$2.$3.mem.stat ~/smartSlurm/stats/$2.$3.time.stat
            else 
                echo Did not add this record to ~/smartSlurm/stats/jobRecord.txt
            fi  
        fi
    else
        # todo: may not need failed job records?
        echo # $record >> ~/smartSlurm/jobRecord.txt
        #echo -e "Added this line to ~/smartSlurm/jobRecord.txt:\n$record"
    fi
    echo Final mem usage: $mem, time usage: $mins minutes
fi    

# for testing
#rm $succFile; jobStatus=OOM 

if [ ! -f $succFile ]; then
    touch $failFile

    #echo looking partition for hour: $hours 
    x=`realpath $0` 
    . ${x%\/bin\/cleanUp.sh}/config/partitions.txt || { echo Partition list file not found: partition.txt; exit 1; }

    if [[ "$jobStatus" == "OOM" ]]; then
        
        

        jobStat=`echo -e "$sacct" | head -n 3 | tail -n 1`

        mem=${jobStat#*mem=}; mem=${mem%M*}


        [ "$mem" -lt 500 ] && mem=500 # at least 500M


        if [ -n "$mem" ] && [ "$mem" -eq "$mem" ] 2>/dev/null; then
            mem=$(( $mem * 2 ))
            
            echo Trying to re-queue the job with memory: $mem
            
            # submit a small job to requeue the job, because it can not requeue itself
#             p=`adjustPartition 1 short`
#             echo /usr/bin/sbatch --parsable -p $p -t 5 -A ${12} --mail-type=NONE --wrap "scontrol requeue $SLURM_JOBID; scontrol update JobId=$SLURM_JOBID MinMemoryNode=$mem;"
#             jobID=`/usr/bin/sbatch -o /dev/null -e /dev/null --parsable --mail-type=NONE -p $p -t 5 -A ${12} --wrap "scontrol requeue $SLURM_JOBID; scontrol update JobId=$SLURM_JOBID MinMemoryNode=$mem;"`
            
            # does work, can not update memory for job steps
            # ( sleep 5; 
            #     scontrol requeue $SLURM_JOBID 
            #     scontrol update JobId=$SLURM_JOBID MinMemoryNode=$mem
            #     scontrol update JobId=$SLURM_JOBID 1 Memory=$(( $mem - 10 )) 
            
            # ) &
            # this works but need ssh to login nodes to run scontrol
            # put this block of code in to background and sleep, so that email run first before requeue 
            ( sleep 5; 
            
            for try in {1..8}; do
                if [ ! -f $failFile.requeued.$try ]; then
                    echo trying to requeue $try
                    touch $failFile.requeued.$try 
                   
                    # 80G memory
                    #[ "$mem" -gt 81920 ] && [ "$try" -gt 2 ] && break
            
        
                    # let's try to requeue using login nodes
                    for nodeIdx in {1..2}; do
                        echo trying ssh to node $nodeIdx
                        if `ssh login00 "scontrol requeue $SLURM_JOBID; scontrol update JobId=$SLURM_JOBID MinMemoryNode=$mem;"`; then 
                            echo Requeued successfully from computer login0$nodeIdx
                            [ -f $failFile ] && rm $failFile
                            break
                        fi    
                    done  
                    #break
                fi
            done ) &
        else
            echo Could not find the original mem value.
            echo Job failed of out-of-memory. Please resubmit with more memory check youself.
        fi  
        
        # delete stats and redo them
        if [ "$5" == 0 ]; then
            rm ~/smartSlurm/stats/$2.$3.mem.stat.noInput ~/smartSlurm/stats/$2.$3.time.stat.noInput 2>/dev/null
        else 
            rm ~/smartSlurm/stats/$2.$3.mem.stat ~/smartSlurm/stats/$2.$3.time.stat 2>/dev/null
        fi 
        #rm ~/smartSlurm/stats/$2.$3.mem.stat* ~/smartSlurm/stats/$2.$3.time.stat* 2>/dev/null
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
            if [ ! -f $failFile.requeued.$try ]; then
                echo trying to requeue $try
                touch $failFile.requeued.$try 

                # 80G memory
                #[ "${mem%M}" -gt 81920 ] && [ "$try" -gt 2 ] && break

                touch $failFile.requeued.$try 
                scontrol requeue $SLURM_JOBID && echo job re-submitted || echo job not re-submitted.

                # time=${10}
                # [[ "$time" == *-* ]] && { day=${time%-*}; tem=${time#*-}; hour=${tem%%:*}; min=${tem#*:}; min=${min%%:*}; sec=${tem#$hour:$min}; sec=${sec#:}; } || { [[ "$time" =~ ^[0-9]+$ ]] && min=$time || { sec=${time##*:}; min=${time%:*}; min=${min##*:}; hour=${time%$min:$sec}; hour=${hour%:}; day=0;} }

                # [ -z "$day" ] && day=0; [ -z "$hour" ] && hour=0; [ -z "$min" ] && min=0;[ -z "$sec" ] && sec=0

                # echo $day $day,  $hour hour,  $min min,  $sec sec

                # # how many hours for sbatch command if we double the time
                # hours=$(($day * 2 * 24 + $hour * 2 + ($min * 2 + 59 + ($sec * 2 + 59) / 60 ) / 60))

                [ "$mins" -lt 20 ] && mins=20 # at least 20 minutes

                hours=$((($mins * 2 + 59) / 60))

                partition=`adjustPartition $hours $partition`

                seconds=$(($mins * 2 * 60))

                #[ $seconds -le 60 ] && time=11:0 && seconds=60

                #echo srun seconds: $seconds 
                
                # https://chat.openai.com/chat/80e28ff1-4885-4fe3-8f21-3556d221d7c6

                time=`eval "echo $(date -ud "@$seconds" +'$((%s/3600/24))-%H:%M:%S')"`

                if [[ "$partition" != "${11}" ]]; then 
                    scontrol update jobid=$SLURM_JOBID Partition=$partition TimeLimit=$time
                else 
                    scontrol update jobid=$SLURM_JOBID TimeLimit=$time
                    scontrol update jobstep=123456.2 TimeLimit=02:00:00

                fi 

                echo job resubmitted: $SLURM_JOBID with time: $time partition: $partition

                [ -f $failFile ] && rm $failFile

                break
            fi
        done 
         # delete stats and redo them
        if [ "$5" == 0 ]; then
            rm ~/smartSlurm/stats/$2.$3.mem.stat.noInput ~/smartSlurm/stats/$2.$3.time.stat.noInput 2>/dev/null
        else 
            rm ~/smartSlurm/stats/$2.$3.mem.stat ~/smartSlurm/stats/$2.$3.time.stat 2>/dev/null
        fi 
    else 
        echo Not sure why job failed. Not run out of time or memory. Pelase check youself.
    fi
elif [ ! -z "$1" ]; then
    adjustDownStreamJobs.sh $1/logs $4     
fi

echo "Sending email..."

minimumsize=9000

actualsize=`wc -c $out`

[ -f $succFile ] && s="Succ:$SLURM_JOBID:$SLURM_JOB_NAME" || s="$jobStatus:$SLURM_JOBID:$SLURM_JOB_NAME" 

if [ "${actualsize% *}" -ge "$minimumsize" ]; then
   toSend=`echo Job script content:; cat $script;`
   toSend="$toSend\nOutput is too big for email. Please find output in: $out"  
   toSend="$toSend\n...\nLast 10 row of output:\n`tail -n 10 $out`"
else
   toSend=`echo Job script content:; cat $script; echo; echo Job output:; cat $out;`
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

#echo -e "tosend:\n$toSend"
echo -e "$toSend" >> ${err%.err}.email

#echo -e "$toSend" | sendmail `head -n 1 ~/.forward`
echo -e "$toSend" | mail -s "$s" $USER && echo email sent || \
    { echo Email not sent.; echo -e "$toSend \nTry again..." | sendmail `head -n 1 ~/.forward` && echo Email sent by second try. || echo Email still not sent!!; }

echo 

# wait for email to be sent
sleep 20

#[ -f $succFile ] && exit 0  

#exit 1; 


