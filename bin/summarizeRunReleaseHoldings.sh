#!/bin/bash

Usage="Usage: $0 [ workDir/log, the log folder name. ]  \nThis script will go through job name list in allJobs.txt to see if the jobs finish successfully or not."

#set -x

echo Running: $0 $@

if [ -f $smartSlurmLogDir/allJobs.txt ]; then 
    lines=`tail -n +2 $smartSlurmLogDir/allJobs.txt` # | awk 'NF>2{print $1, $2, $3}'`
else 
    exit 1; 
    #lines="$SLURMJOB_ID $2" # tail -n +2 allJobs.txt | awk 'NF>2{print $1, $2, $3}'`
fi 
 
echo pwd: `pwd`

IFS=$'\n'

out=`squeue -u $USER -t PD,R --noheader -o "%.18i-%t"`

#echo squeue out: $out

current=0; succ=0; fail=0; running=0; pending=0; requeue=0; unknown=0; unholdCounter=0; 
toSend="Summery for jobs in allJobs.txt:"
for line in $lines; do
    if [ ! -z "${line/ /}" ]; then
        #id=${line%% *}; name=${line##* }
        IFS=' ' read -a arrIN <<< "$line"

        id=${arrIN[0]}; 
        [[ $id == $SLURM_JOBID ]] && currentName=${arrIN[2]} && currentStep=$(echo "$currentName" | cut -d '.' -f 1-2)

        deps=${arrIN[1]}
        name=${arrIN[2]}
        program=${arrIN[3]}
        ref=${arrIN[4]}  ; ref=${ref//\//-}
        inputs="${arrIN[5]}"

        if [ -f $smartSlurmLogDir/$name.success ]; then
            toSend="$toSend\n${line:0:60} Done"
            succ=$((succ + 1))
        elif [ -f $smartSlurmLogDir/$name.failed ]; then
            toSend="$toSend\n${line:0:60} Failed"
            fail=$((fail + 1))
        elif [[ "$out" == *$id-R* ]]; then # && [[ "$id" != "$SLURM_JOBID" ]]; then
            toSend="$toSend\n${line:0:60} Running"
            running=$((running + 1))
        elif [[ "$out" == *$id-P* ]]; then 
            toSend="$toSend\n${line:0:60} Pending"
            pending=$((pending + 1))
            if [ $unholdCounter -gt 0 ]; then 
                thisStep=$(echo "$name" | cut -d '.' -f 1-2)
                if [[ $thisStep == $currentStep ]] && [[ "$deps" == null ]] && grep " \-H " $smartSlurmLogDir/$name.sh && [ ! -f $smartSlurmLogDir/$name.adjust ]; then 
                    unholdCounter=$((unholdCounter - 1))
                    
                    echo Trying to release: $name 

                    echo Job is beging released by $SLURM_JOBID $currentName >> $smartSlurmLogDir/$name.out

                    IFS=' ' read -r inputSize mem min extraMem <<< $(estimateResource.sh $program ${ref//\//-} $inputs $name 0 0 adjust)
                    
                    if [[ $mem != 0 ]] && [[ $min != 0 ]]; then
                        #set -x

                        hours=$((( min + 59 ) / 60 ))

                        echo looking partition for hour: $hours

                        adjustPartition $hours partition

                        seconds=$(( min * 60 ))

                        time=`eval "echo $(date -ud "@$seconds" +'$((%s/3600/24))-%H:%M:%S')"`

                        #set -x 

                        #scontrol show job $id

                        echo running: scontrol update jobid=$id timelimit=$time partition=$partition MinMemoryNode=${mem}

                        echo running: scontrol update jobid=$id timelimit=$time partition=$partition MinMemoryNode=${mem} >> $smartSlurmLogDir/$name.out

                        scontrol update JobId=$id TimeLimit=$time Partition=$partition  MinMemoryNode=${mem}
                        #scontrol show job $id

                        scontrol release $id
                        echo -e "Job is beging released by $SLURM_JOBID $currentName with new mem and time\n" >> $smartSlurmLogDir/$name.out
                        set +x
                    # didn't get estimate, but already have 3 successful jobs, release one job anyway
                    # because jobRecord need to be unique by programName + reference + inputSize + memory
                    # If the first 5 jobs have same unique value, there is only one record in jobRecords.txt
                    else
                        unholdCounter=0; 
                        echo Fail to estimate new resource. Directly release one job anyway. 
                        echo -e "Job is beging direcly released by $SLURM_JOBID $currentName\n" >> $smartSlurmLogDir/$name.out
                        scontrol release $id 
                    fi 
                    
                fi
            fi            
        elif [ -f $smartSlurmLogDir/$name.failed.requeued.1.time ]; then 
            toSend="$toSend\n${line:0:40} Requeued"
            requeue=$((requeue + 1))    
        else
            toSend="$toSend\n${line:0:40} Unknow"
            unknown=$((unknown + 1))
        fi
        if [ "$id" == "$SLURM_JOBID" ]; then 
            #check if statics is available for new ten hoding jobs
            unholdCounter=5 # only take care of jobs without dependency
        fi 
    fi
done

current=$((succ + fail + requeue))
total=$((succ + fail + running + pending + +requeue + unknown))
s="$current/$total Succ:$succ/$total Requeue:$requeue/$total Running:$running/$total Pending:$pending/$total Fail:$fail/$total Unknown:$unknown/$total"

[ -f $smartSlurmLogDir/allJobs.txt ] && echo -e "$s\n$toSend" > $smartSlurmLogDir/summary.$SLURMJOB_ID

# if [ $((running + pending)) -le 5 ]; then
#     echo -e "$toSend" | mail -s $s $USER
#     [ "$USER" != ld32 ] && echo -e "$toSend" | mail -s $s ld32
# fi


