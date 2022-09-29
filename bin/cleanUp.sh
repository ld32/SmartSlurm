#!/bin/bash

set -x 

# to call this:  0     1           2           3       4         5          6      7        8       9  10   11
#emailAndRecord.sh "projectDir"  "$software" "$ref" "$flag" "$inputSize"   $core $memO  $timeO  $mem  $time  partition
                
para="$USER $2 $3 $4 $5 $6 $7 $8 $9 $10";  out=$1/flag/"${4}.out"; err=$1/flag/"${4}.err"; script=$1/flag/"${4}.sh"; succFile=$1/flag/"${4}.success";   

sacct=`sacct --format=JobID,Submit,Start,End,MaxRSS,State,NodeList%30,Partition,ReqTRES%30,TotalCPU,Elapsed%14,Timelimit%14 --units=M -j $SLURM_JOBID` 

echo -e "\nJob summary:\n$sacct"
# echo *Notice the sacct report above: while the main job is still running for sacct command, user task is completed.



# record job for future estimating mem and time
jobStat=`echo -e "$sacct" | tail -n 1`

#echo Running: $0  $@

echo -e  "Last row of job summary: $jobStat" 
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

#echo start: $START finish: $FINISH mem: $mem mins: $mins

failReason=""

case "$jobStat" in
# jobRecord.txt header
#1user 2software 3ref 4inputName 5inputSizeInK 6CPUNumber 7memoryO 8timeO 9readMem 10RequestedTime 11jobID 12memoryM 13minRun 14Node 15 finalStatus    
*COMPLETED* )  record="$para $SLURM_JOB_ID  ${mem%M} ${mins} ${node} COMPLETED `date`" && echo *Notice the sacct report above: while the main job is still running for sacct command, user task is completed.;;

*TIMEOUT*   )  record="$para $SLURM_JOB_ID  ${mem%M} ${mins} ${node} needMoreTime $errFlag `date`" && failReason="(needMoreTime)";;
        
*CANCELLED*	)  grep "Exceeded job memory limit" $out $err 2>&1 >/dev/null && failReason="(needMoreMem)" && record="$para $SLURM_JOB_ID ${mem%M} ${mins} ${node} needMoreMem $flag.out `date`" || record="$para $SLURM_JOB_ID ${mem%M} ${mins} ${node} UnknowReason $errFlag `date`";;

esac
    
if [[ "$5" != "0" && -z "$failReason" && "$mem" != "0" && ! -z "$record" && ! -f ~/.rcbio/$1.$2.mem.stat.final ]; then 
    echo $record >> ~/.smartSlurm/myJobRecord.txt
    echo -e "Added this line to ~/.smartSlurm/myJobRecord.txt:\n$record"
    
else 
#    echo "Job record:\n$record\n" 
    echo Did not add this record to ~/.rcbio/myJobRecord.txt
#    echo Because we already have ~/.rcbio/$1.$2.mem.stat.final
fi

# do we need calculate stats here???

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

#to=`cat ~/.forward`
#echo -e "$s\n$toSend" | sendmail $to && echo email sent || echo email not sent

#adjustDownStreamJobs.sh $1/flag $4 

[ -f $succFile ] && exit 0  

touch $failFile

scontrol requeue $SLURM_JOBID 

#todo: check if out of time? or out of memory

if [[ "$failReason" == "(needMoreTime)" ]]; then
 
    time=$10
    [[ "$time" == *-* ]] && { day=${time%-*}; tem=${time#*-}; hour=${tem%%:*}; min=${tem#*:}; min=${min%%:*}; sec=${tem#$hour:$min}; sec=${sec#:}; } || { [[ "$time" =~ ^[0-9]+$ ]] && min=$time || { sec=${time##*:}; min=${time%:*}; min=${min##*:}; hour=${time%$min:$sec}; hour=${hour%:}; day=0;} }

    [ -z "$day" ] && day=0; [ -z "$hour" ] && hour=0; [ -z "$min" ] && min=0;[ -z "$sec" ] && sec=0

    echoerr $day $day,  $hour hour,  $min min,  $sec sec

    # how many hours for sbatch command if we double the time
    hours=$(($day * 2 * 24 + $hour * 2 + ($min * 2 + 59 + ($sec * 2 + 59) / 60 ) / 60))

    #echoerr looking partition for hour: $hours 
    x=`realpath $0` 
    . ${x%\/bin\/smartSbatch}/config/partitions.txt || { echoerr Partition list file not found: partition.txt; exit 1; }

    partition=`adjustPartition $hours $partition`

    timeN=...

    if [[ "$partition" != "$11" ]]; then 
        scontrol update jobid=$SLURM_JOBID Partition=$partition
    fi 



elif [[ "$failReason" == "(needMoreMem)" ]]; then
    


scontrol update jobid=$SLURM_JOBID TimeLimit=0:10:0 MinMemoryNode=40

echo job resubmitted: $SLURM_JOBID 

exit 1; 


