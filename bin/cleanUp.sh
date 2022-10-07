#!/bin/bash

#set -x 

# to call this:  0     1           2           3       4         5          6       7        8       9    10      11       12
#cleanUp.sh       "projectDir"  "$software" "$ref" "$flag" "$inputSize"   $core   $memO  $timeO    $mem  $time  $partition slurmAcc

echo Running $0 $@
                
para="$USER $2 $3 $4 $5 $6 $7 $8 $9 ${10}";  out=$1/flag/"${4}.out"; err=$1/flag/${4}.err; script=$1/flag/${4}.sh; succFile=$1/flag/${4}.success; failFile=$1/flag/${4}.failed;   

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

failReason=""

case "$jobStat" in
# jobRecord.txt header
#1user 2software 3ref 4inputName 5inputSizeInK 6CPUNumber 7memoryO 8timeO 9readMem 10RequestedTime 11jobID 12memoryM 13minRun 14Node 15 finalStatus    
*COMPLETED* )  record="$para $SLURM_JOB_ID  ${mem%M} ${mins} ${node} COMPLETED `date`" && echo *Notice the sacct report above: while the main job is still running for sacct command, user task is completed.;;

*TIMEOUT*   )  record="$para $SLURM_JOB_ID  ${mem%M} ${mins} ${node} needMoreTime $errFlag `date`" && failReason="(needMoreTime)";;

*OUT_OF_ME*   ) record="$para $SLURM_JOB_ID ${mem%M} ${mins} ${node} needMoreMem $errFlag `date`" && failReason="(needMoreMem)";;
        
*CANCELLED*	) record="$para $SLURM_JOB_ID ${mem%M} ${mins} ${node} Cancelled $errFlag `date`" && failReason="(cancelled)";;

esac

echo -e  "Last row of job summary: $jobStat" 
echo start: $START finish: $FINISH mem: $mem mins: $mins
echo failReason: $failReason
    
if [[ "$5" != "0" && -z "$failReason" && "${mem%M}" != "0" && ! -z "$record" && ! -f ~/.smartSlurm/$2.$3.mem.stat.final ]]; then 
    echo $record >> ~/.smartSlurm/myJobRecord.txt
    echo -e "Added this line to ~/.smartSlurm/myJobRecord.txt:\n$record"
    
else 
#    echo "Job record:\n$record\n" 
    echo Did not add this record to ~/.smartSlurm/myJobRecord.txt
#    echo Because we already have ~/.smartSlurm/$1.$2.mem.stat.final
fi

if [ ! -f $succFile ]; then
    touch $failFile

    #echo looking partition for hour: $hours 
    x=`realpath $0` 
    . ${x%\/bin\/cleanUp.sh}/config/partitions.txt || { echo Partition list file not found: partition.txt; exit 1; }

    if [[ "$failReason" == "(needMoreTime)" ]]; then
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

    elif [[ "$failReason" == "(needMoreMem)" ]]; then

        jobStat=`echo -e "$sacct" | head -n 3 | tail -n 1`

        mem=${jobStat#*mem=}; mem=${mem%M*}

        if [ -n "$mem" ] && [ "$mem" -eq "$mem" ] 2>/dev/null; then
            echo Submitting a job to re-queue the job. 
            mem=$(( $mem * 2 ))
            p=`adjustPartition 1 short`
            echo sbatch --parsable -p $p -t 5 -A ${12} --wrap "scontrol requeue $SLURM_JOBID; scontrol update JobId=$SLURM_JOBID MinMemoryNode=$mem;"
            jobID=`sbatch --parsable --mail-type=ALL -p $p -t 5 -A ${12} --wrap "scontrol requeue $SLURM_JOBID; scontrol update JobId=$SLURM_JOBID MinMemoryNode=$mem;"`
            #scontrol top $jobID   
        else
            echo Could not find the original mem value.
            echo Job failed of out-of-memory. Please resubmit with more memory check youself.
        fi

    
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
        
    else 
        echo Not sure why job failed. Not run out of time or memory. Pelase check youself.
    fi
else
    adjustDownStreamJobs.sh $1/flag $4     
fi

minimumsize=9000

actualsize=`wc -c $out`

[ -f $succFile ] && s="Success: job id:$SLURM_JOBID name:$SLURM_JOB_NAME" || s="Failed$failReason: job id:$SLURM_JOBID name:$SLURM_JOB_NAME" 

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

[ -f $succFile ] && exit 0  

exit 1; 


