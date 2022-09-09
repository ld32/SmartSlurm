#!/bin/bash

flag=$1

minimumsize=9000

actualsize=`wc -c $flag.out`

[ ! -f $flag.success ] && s="Subject: Failed: job id:$SLURM_JOBID name:$SLURM_JOB_NAME" ||  s="Subject: Success: job id:$SLURM_JOBID name:$SLURM_JOB_NAME"

stat=`tail -n 1 $flag.out`
[[ "$stat" == *COMPLETED* ]] && echo *Notice the sacct report above: while the main job is still running for sacct command, user task is completed. >> $flag.out

if [ "${actualsize% *}" -ge "$minimumsize" ]; then
   toSend=`echo Job script content:; cat $flag.sh`
   toSend="$toSend\nOutput is too big for email. Please find output in: $flag.out"  
   toSend="$toSend\n...\n`tail -n 6 $flag.out`"
else
   toSend=`echo Job script content:; cat $flag.sh; echo Job output:; cat $flag.out`
   #toSend="$s\n$toSend"
fi

#echo -e "$toSend" | sendmail $USER 
echo -e "$toSend" | mail -s "$s" $USER

echo Job done. Summary:;

sacct --format=JobID,Submit,Start,End,State,Partition,ReqTRES%20,Timelimit,CPUTime,MaxRSS,NodeList%30 --units=M -j $SLURM_JOBID;

[ ! -f $flag.success ] && exit 1 || exit 0

