#!/bin/bash

#set -x 

# to call this:  0     1           2           3       4         5          6       7        8       9    10      11       12           13 
#cleanUp.sh       "projectDir"  "$software" "$ref" "$flag" "$inputSize"   $core   $memO  $timeO    $mem  $time  $partition slurmAcc  original.sbatch.command

echo Running $0 $@

if [ -z "$1" ]; then 

    out=${4%% *}; out=${out/\%j/$SLURM_JOB_ID}; err=${4##* }; err=${err/\%j/$SLURM_JOB_ID}; script=${4% *}; script=${script#* }; succFile=${script/\.sh/}.success; failFile=${script/\.sh/}.failed; 
    
else 
    out=$1/logs/"${4}.out"; err=$1/logs/${4}.err; script=$1/logs/${4}.sh; succFile=$1/logs/${4}.success; failFile=$1/logs/${4}.failed;   
fi 

sleep 5

sacct=`sacct --format=JobID,Submit,Start,End,MaxRSS,State,NodeList%30,Partition,ReqTRES%30,TotalCPU,Elapsed%14,Timelimit%14 --units=M -j $SLURM_JOBID` 

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

[[ $jobStatus != "COMPLETED" ]] && [ -f $succFile ] && rm $succFile

record="$SLURM_JOB_ID,$5,$7,$8,$9,${10},${mem%M},${mins},$jobStatus,$USER,$1,$2,$3,$4,$6,${node},$err,`date`,\"${13}\""

echo -e  "Last row of job summary: $jobStat" 
echo start: $START finish: $FINISH mem: $mem mins: $mins
echo jobStatus: $jobStatus
    
if [ ! -z "$record" ]; then
#    if [[ ! -f ~/smartSlurm/stats/$2.$3.mem.stat.final || "$2" == "regularSbatch" ]]; then 
        echo $record >> ~/smartSlurm/myJobRecord.txt
        echo -e "Added this line to ~/smartSlurm/myJobRecord.txt:\n$record"
#    else 
#        echo Did not add this record to ~/smartSlurm/stats/myJobRecord.txt
#    fi
#else 
#    echo "Job record:\n$record\n" 
#    echo Did not add this record to ~/smartSlurm/stats/myJobRecord.txt1
#    echo Because we already have ~/smartSlurm/stats/$1.$2.mem.stat.final
fi

if [ ! -f $succFile ]; then
    touch $failFile

    #echo looking partition for hour: $hours 
    x=`realpath $0` 
    . ${x%\/bin\/cleanUp.sh}/config/partitions.txt || { echo Partition list file not found: partition.txt; exit 1; }

    if [[ "$jobStatus" == "OOM" ]]; then

        jobStat=`echo -e "$sacct" | head -n 3 | tail -n 1`

        mem=${jobStat#*mem=}; mem=${mem%M*}

        if [ -n "$mem" ] && [ "$mem" -eq "$mem" ] 2>/dev/null; then
            echo Submitting a job to re-queue the job. 
            mem=$(( $mem * 2 ))
            p=`adjustPartition 1 short`
            echo /usr/bin/sbatch --parsable -p $p -t 5 -A ${12} --mail-type=NONE --wrap "scontrol requeue $SLURM_JOBID; scontrol update JobId=$SLURM_JOBID MinMemoryNode=$mem;"
            jobID=`/usr/bin/sbatch -o /dev/null -e /dev/null --parsable --mail-type=NONE -p $p -t 5 -A ${12} --wrap "scontrol requeue $SLURM_JOBID; scontrol update JobId=$SLURM_JOBID MinMemoryNode=$mem;"`
            #scontrol top $jobID   
        else
            echo Could not find the original mem value.
            echo Job failed of out-of-memory. Please resubmit with more memory check youself.
        fi

#        [ -f ~/smartSlurm/stats/$2.$3.mem.stat.final ] && rm ~/smartSlurm/stats/$2.$3.mem.stat.final ~/smartSlurm/stats/$2.$3.time.stat.final  
        if [ "$5" == 0 ]; then
            rm ~/smartSlurm/stats/$2.$3.mem.stat.noInput ~/smartSlurm/stats/$2.$3.time.stat.noInput 2>/dev/null
        else 
            rm ~/smartSlurm/stats/$2.$3.mem.stat.final ~/smartSlurm/stats/$2.$3.time.stat.final 2>/dev/null
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
        scontrol requeue $SLURM_JOBID && echo job re-submitted || echo job not re-submitted.
    
        # time=${10}
        # [[ "$time" == *-* ]] && { day=${time%-*}; tem=${time#*-}; hour=${tem%%:*}; min=${tem#*:}; min=${min%%:*}; sec=${tem#$hour:$min}; sec=${sec#:}; } || { [[ "$time" =~ ^[0-9]+$ ]] && min=$time || { sec=${time##*:}; min=${time%:*}; min=${min##*:}; hour=${time%$min:$sec}; hour=${hour%:}; day=0;} }

        # [ -z "$day" ] && day=0; [ -z "$hour" ] && hour=0; [ -z "$min" ] && min=0;[ -z "$sec" ] && sec=0

        # echo $day $day,  $hour hour,  $min min,  $sec sec

        # # how many hours for sbatch command if we double the time
        # hours=$(($day * 2 * 24 + $hour * 2 + ($min * 2 + 59 + ($sec * 2 + 59) / 60 ) / 60))

        hours=$((($mins * 2 + 59) / 60))

        partition=`adjustPartition $hours $partition`

        seconds=$(($mins * 2 * 60))

        #[ $seconds -le 60 ] && time=11:0 && seconds=60

        #echo srun seconds: $seconds

        time=`eval "echo $(date -ud "@$seconds" +'$((%s/3600/24))-%H:%M:%S')"`

        if [[ "$partition" != "${11}" ]]; then 
            scontrol update jobid=$SLURM_JOBID Partition=$partition TimeLimit=$time
        else 
            scontrol update jobid=$SLURM_JOBID TimeLimit=$time
            
        fi 

        echo job resubmitted: $SLURM_JOBID with time: $time partition: $partition
        
        if [ "$5" == 0 ]; then
            rm ~/smartSlurm/stats/$2.$3.mem.stat.noInput ~/smartSlurm/stats/$2.$3.time.stat.noInput 2>/dev/null
        else 
            rm ~/smartSlurm/stats/$2.$3.mem.stat.final ~/smartSlurm/stats/$2.$3.time.stat.final 2>/dev/null
        fi    
#        [ -f ~/smartSlurm/stats/$2.$3.mem.stat.final ] && rm ~/smartSlurm/stats/$2.$3.mem.stat.final ~/smartSlurm/stats/$2.$3.time.stat.final      
        
    else 
        echo Not sure why job failed. Not run out of time or memory. Pelase check youself.
    fi
elif [ ! -z "$1" ]; then
    adjustDownStreamJobs.sh $1/logs $4     
fi

minimumsize=9000

actualsize=`wc -c $out`

[ -f $succFile ] && s="Success: job id:$SLURM_JOBID name:$SLURM_JOB_NAME" || s="Failed($jobStatus): job id:$SLURM_JOBID name:$SLURM_JOB_NAME" 

if [ "${actualsize% *}" -ge "$minimumsize" ]; then
   toSend=`echo Job script content:; cat $script;`
   toSend="$toSend\nOutput is too big for email. Please find output in: $out"  
   toSend="$toSend\n...\nLast 6 row of output:\n`tail -n 6 $out`"
else
   toSend=`echo Job script content:; cat $script; echo; echo Job output:; cat $out;`
   #toSend="$s\n$toSend"
fi

if [ -f "$err" ]; then 
    actualsize=`wc -c $err`
    if [ "${actualsize% *}" -ge "$minimumsize" ]; then
        toSend="$toSend\nError file is too big for email. Please find output in: $err"  
        toSend="$toSend\n...\nLast 6 rows of error file:\n`tail -n 6 $err`"
    else
        toSend="$toSend\n\nError output:\n`cat  $err`"
    fi
fi    

#echo -e "tosend:\n $toSend"
echo "Sending email..."
#echo -e "$toSend" | sendmail `head -n 1 ~/.forward`
echo -e "$toSend" | mail -s "$s" $USER && echo email sent || echo email not sent

echo 

# wait for email to be sent
sleep 20

#[ -f $succFile ] && exit 0  

#exit 1; 


