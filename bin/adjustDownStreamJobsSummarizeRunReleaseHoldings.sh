#!/bin/bash

#set -x

Usage="Usage: $0 full_path_to_flag_folder \n  Note: this script will go through job id list file, find the downstream jobs, and return them as a string of job flags. "

echo Running: $0  $@


[ -f $smartSlurmLogDir/allJobs.txt ] || { echo -e "job id file $smartSlurmLogDir/allJobs.txt does not exist\n$Usage"; exit 1; }

# jobid, deps, flag, software, ref, input, inputSize
if [ -f $smartSlurmLogDir/allJobs.txt ]; then 
    lines=`grep -v ^job_id $smartSlurmLogDir/allJobs.txt` # | awk 'NF>2{print $1, $2, $3}'`
else 
    exit 1; 
fi 

IFS=$' ';

# directly get id, deps, software, ref, and input here, if input is none, directly skip this job
output=`echo $lines | awk '{if ($2 ~ /'"$SLURM_JOBID/"') print $1, $2, $3, $4, $5, $6;}'`

echo -e "Jobs on the same dependency level with current job:\n$output"
if [ -z "$output" ]; then 
    echo -e "Downstream job ids not found for $SLURM_JOBID"; #summarizeRunReleaseHoldings.sh; exit; }

else 

    IFS=$'\n';
    for i in $output; do
        echo 1working on $i
        IFS=' ' read -a arrIN <<< "$i"

        id=${arrIN[0]}
        deps=${arrIN[1]}
        name=${arrIN[2]}
        program=${arrIN[3]}
        ref=${arrIN[4]}  ; ref=${ref//\//-}
        inputs=${arrIN[5]}

        allDone=""
        IFS=$' ';
        for j in ${deps//:/ }; do
            echo 2working on $j
            [[ "$j" == "$SLURM_JOBID" ]] && continue; # ignore current job
            #echo look for the job flag for $j
            job=`echo $lines | awk '{if ($1 ~ /'"$j/"') print $3;}'`
            #[ -z "$job" ] && { echo -e "job name not found!"; exit; }
            [ -f "$smartSlurmLogDir/$job.success" ] && echo This job was done! $job || { allDone=no; break;}
        done

        if [ -z "$allDone" ]; then
            date
            #ls -lrt $smartSlurmLogDir
            echo Dependants for $name are all done. Ready to adjust mem/runtime...

            echo -e "Re-adjust resource by upsteam job job $SLURM_JOB_ID:" >> $smartSlurmLogDir/$name.out
            grep ^$SLURM_JOB_ID $smartSlurmLogDir/allJobs.txt | awk '{print $1,  $2,  $3}' >> $smartSlurmLogDir/$name.out

            if [ -f $smartSlurmLogDir/$name.adjust ]; then
                echo Already adjusted? Directly release and run.
                scontrol release $id
                continue
            fi

            IFS=' ' read -r inputSize mem min extraMem <<< `estimateResource.sh $program ${ref//\//-} $inputs $name 0 0 adjust`

            #inputSize=${output%% *}; mem=${output% *}; mem=${mem#* }; min=${output##* }

            [[ "$mem" == 0 ]] && echo Fail to estimate new resource. Directly release job. && scontrol release $id && continue

            #[ "$mem" -lt 20 ] && mem=20 # at least 20M

            #echo Got estimation inputsize: $inputSize mem: $mem  time: $min

            #echo Got estimation inputsize: $inputSize mem: $mem  time: $min  >> $smartSlurmLogDir/$name.out
            hours=$((($min + 59) / 60))

            echo looking partition for hour: $hours

            adjustPartition $hours partition

            seconds=$(($min * 60))

            time=`eval "echo $(date -ud "@$seconds" +'$((%s/3600/24))-%H:%M:%S')"`

            #set -x 

            #scontrol show job $id

            echo running: scontrol update jobid=$id timelimit=$time partition=$partition MinMemoryNode=${mem}

            echo running: scontrol update jobid=$id timelimit=$time partition=$partition MinMemoryNode=${mem} >> $smartSlurmLogDir/$name.out

            scontrol update JobId=$id TimeLimit=$time Partition=$partition  MinMemoryNode=${mem}
            #scontrol show job $id

            scontrol release $id

            #scontrol show job $id

            #set +x 
        else
            echo Need wait for other jobs to finish before we can ajust mem and runtime...
        fi
    done
fi
#rm -r $smartSlurmLogDir/downsteamjob.adjusting 2>/dev/null

#summarizeRunReleaseHoldings.sh    

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
            #check if statics is available for new ten holding jobs
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

if [ $running -eq 0 ] && [ $pending -eq 0 ]; then 
    
    echo This should be the last job for this run... >> $smartSlurmLogDir/summary.$SLURMJOB_ID
    
else 
    rm $smartSlurmLogDir/summary.run 2>/dev/null
fi 



