#!/bin/bash

#set -x 

# to call this:  0     1         2        3         4          5     6     7       8    9    
#emailAndRecord.sh "$software" "$ref" "$flag" "$inputSize"   $core $memO  $timeO  $mem  $time 
                
para="$USER $1 $2 $3 $4 $5 $6 $7 $8 $9";  out=flag/"${3}.out"; err=flag/"${3}.err"; script=flag/"${3}.sh"; succFile=flag/"${3}.success";   

if [[ "$4" != "0" ]]; then

    # record job for future estimating mem and time
    jobStat=`grep "Job done. Summary:" -A 6 $out | tail -n 1`

    #echo Running: $0  $@

    #echo -e  "Last row of job summary: $jobStat" 
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
        
    case "$jobStat" in
    # jobRecord.txt header
    #1user 2software 3ref 4inputName 5inputSizeInK 6CPUNumber 7memoryO 8timeO 9readMem 10RequestedTime 11jobID 12memoryM 13minRun 14Node 15 finalStatus    
    *COMPLETED* )  record="$para $SLURM_JOB_ID  ${mem%M} ${mins} ${node} COMPLETED `date`" && echo *Notice the sacct report above: while the main job is still running for sacct command, user task is completed.;;
    
    *TIMEOUT*   )  record="$para $SLURM_JOB_ID  ${mem%M} ${mins} ${node} needMoreTime $errFlag `date`" && failReason="(needMoreTime)";;
            
    *CANCELLED*	)  grep "Exceeded job memory limit" $out $err 2>&1 >/dev/null && failReason="(needMoreMem)" && record="$para $SLURM_JOB_ID ${mem%M} ${mins} ${node} needMoreMem $flag.out `date`" || record="$para $SLURM_JOB_ID ${mem%M} ${mins} ${node} UnknowReason $errFlag `date`";;

    esac
        
    #if [ ! -f ~/.rcbio/$1.$2.mem.stat.final ]; then 
    
        [ ! -z "$record" ] && echo $record >> ~/.smartSlurm/myJobRecord.txt &&  echo -e "Added this line to ~/.smartSlurm/myJobRecord.txt:\n$record"
    #else 
    #    echo "Job record:\n$record\n" 
    #    echo Did not add this record to ~/.rcbio/myJobRecord.txt
    #    echo Because we already have ~/.rcbio/$1.$2.mem.stat.final
    #fi
fi
# do we need calculate stats here???

minimumsize=9000

actualsize=`wc -c $out`

[ -f $succFile ] && s="Subject: Success: job id:$SLURM_JOBID name:$SLURM_JOB_NAME" || s="Subject: Failed$failReason: job id:$SLURM_JOBID name:$SLURM_JOB_NAME" 

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

echo -e "tosend:\n $toSend"
echo "Sending email..."
#echo -e "$toSend" | sendmail `head -n 1 ~/.forward`
echo -e "$toSend" | mail -s "$s" $USER
