#!/bin/bash

#set -x 

usage() { echo -e "Usage: \n${0##*/} <Job outout file name><Job error output file name><Subject for the email>"; exit 1; }

echo Running: $0 $@

[ -z "$1" ] || [[ "-h" == "$1" ]] || [[ "--help" == "$1" ]] && usage

out=$1; err=$2; title="$3"

totalM=$SLURM_MEM_PER_NODE

# if job has --mem-per-cpu and -c
[ -z "$totalM" ] && totalM=$((SLURM_MEM_PER_CPU * SLURM_JOB_CPUS_PER_NODE))

# if job has --mem-per-cpu and -n
[ -z "$totalM" ] &&  totalM=$((SLURM_MEM_PER_CPU * SLURM_CPUS_PER_TASK))

# wait for slurm database update
sleep 5

sacct=`sacct --format=JobID,Submit,Start,End,MaxRSS,State,NodeList%30,Partition,ReqTRES%30,TotalCPU,Elapsed%14,Timelimit%14,ExitCode --units=M -j $SLURM_JOBID`

echo -e "\nJob summary:\n$sacct"
echo *Notice the sacct report above: while the main job is still running for sacct command, user task is completed.

# record job for future estimating mem and time
jobStat=`echo -e "$sacct" | tail -n 1`

START=`head -n 1 $smartSlurmLogDir/job_$SLURM_JOB_ID.memCPU.txt | cut -d' ' -f6`

FINISH=`date +%s`

echo  start: $START fisnish: $FINISH

# time in minutes
min=$((($FINISH - $START + 59)/60))

[[ "$mim" == 0 ]] && min=1

# node
node=`echo $jobStat | cut -d" " -f7`

if [[ "$sacct" == *TIMEOUT* ]]; then
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

[ -f $smartSlurmLogDir/job_$SLURM_JOBID.succ ] && jobStatus=COMPLETED

echo jobStatus: $jobStatus

srunM=`cut -d' ' -f2 $smartSlurmLogDir/job_$SLURM_JOBID.memCPU.txt | sort -n | tail -n1`

#srunM=$((srunM / 1024 / 1024 ))

[ -z "$srunM" ] && srunM=0

echo jobStatus: $jobStatus srunM: $srunM


overReserved=""; overText=""; ratioM=""; ratioT=""

if [[ "$jobStatus" == COMPLETED ]]; then
    if [ "$totalM" -gt $((srunM * 2)) ] && [ $srunM -ge 1024 ]; then
        overReserved="O:"
        overText="$SLURM_JOBID over-reserved resounce M: $srunM/$totalM T: $min/$totalT"
	#echo -e "`pwd`\n$out\n$SLURM_JOBID over-reserved resounce M: $srunM/$totalM T: $min/$totalT" | mail -s "Over reserved $SLURM_JOBID" $USER
        echo $overText
    fi
  
  #ratioM=`echo "scale=2;$srunM/$totalM"|bc`; ratioT=`echo "scale=2;$min/$totalT"|bc`

fi

echo "Sending email..."

minimumsize=9000

actualsize=`wc -c $out || echo 0`

[ -f $smartSlurmLogDir/job_$SLURM_JOBID.succ ] && s="${overReserved}Succ:$SLURM_JOBID: $title" || s="$jobStatus:$SLURM_JOBID:$title"

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

cd $smartSlurmLogDir



rowTotal=`wc -l job_$SLURM_JOBID.memCPU.txt | cut -d' ' -f1`
if [ "$rowTotal" -gt 50 ]; then 
    maxMem=0; maxCpu=0; 
    rate=`echo "scale=2;$rowTotal/50"|bc`
    IFS=$'\n'; rowCount1=0; rowCount2=0
    echo > job_$SLURM_JOBID.memCPU1.txt
    for t in `cat job_$SLURM_JOBID.memCPU.txt`; do
        mem=`echo $t | cut -d' ' -f2`
        cpu=`echo $t | cut -d' ' -f5`
        [ "$mem" -gt $maxMem ] && maxMem=$mem && mem1=`echo $t | cut -d' ' -f3,4`
        [ "$cpu" -gt $maxCpu ] && maxCpu=$cpu
        rowCount1=$((rowCount1 + 1))
        rowMax=`echo "scale=2;$rowCount2*$rate"|bc`
        rowMax=${rowMax%.*}; [ -z "$rowMax" ] && rowMax=1; 
        if [ "$rowMax" -le "$rowCount1" ]; then 
            rowCount2=$((rowCount2 + 1))
            echo $rowCount2 $maxMem $mem1 $maxCpu >> job_$SLURM_JOBID.memCPU1.txt
            maxMem=0; maxCpu=0;
        fi 
    done
else 
    ln -s job_$SLURM_JOBID.memCPU.txt job_$SLURM_JOBID.memCPU1.txt
fi 

# time vs. memory for current job
gnuplot -e "set key outside; set key reverse; set key invert; set datafile separator ' '; set style data histogram; set style histogram rowstacked gap 2; set style fill solid border rgb 'black'; set xtics rotate by -45; set terminal png size 800,600; set output 'job_$SLURM_JOBID.mem.png'; set title 'Time vs. Mem for job $SLURM_JOBID'; set xlabel 'Time'; set ylabel 'Mem (M)'; plot 'job_$SLURM_JOBID.memCPU1.txt' using 2:xtic(1) title 'Used' lc rgb 'green', '' using 3:xtic(1) title 'Wasted' lc rgb 'red'"

# time vs. CPU usage for current job
gnuplot -e "set key outside; set key reverse; set key invert; set datafile separator ' '; set style data histogram; set style histogram rowstacked gap 2; set style fill solid border rgb 'black'; set xtics rotate by -45; set terminal png size 800,600; set output 'job_$SLURM_JOBID.cpu.png'; set title 'Time vs. CPU Usage for job $SLURM_JOBID'; set xlabel 'Time'; set ylabel 'CPU Usage (%)'; plot 'job_$SLURM_JOBID.memCPU1.txt' using 5:xtic(1) title 'Used' lc rgb 'green'"

#s="${toSend%% *} $s"

echo -e "$toSend" | mail -s "$s" -a job_$SLURM_JOBID.mem.png -a job_$SLURM_JOBID.cpu.png $USER && echo email sent || \
{ echo Email not sent.; echo -e "$s\n$toSend" | sendmail `head -n 1 ~/.forward` && echo Email sent by second try. || echo Email still not sent!!; }

if [[ $USER != ld32 ]]; then
    echo -e "$toSend" | mail -s "$s" -a job_$SLURM_JOBID.mem.png -a job_$SLURM_JOBID.cpu.png ld32
fi

echo To see the plot:
echo display $smartSlurmLogDir/job_$SLURM_JOBID.mem.png
echo display $smartSlurmLogDir/job_$SLURM_JOBID.cpu.png

sleep 5