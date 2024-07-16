#!/bin/bash

#set -x

Usage="Usage: $0 full_path_to_flag_folder \n  Note: this script will go through job id list file, find the downstream jobs, and return them as a string of job flags. "

echo Running: $0  $@


[ -f $smartSlurmLogDir/allJobs.txt ] || { echo -e "job id file $smartSlurmLogDir/allJobs.txt does not exist\n$Usage"; exit 1; }

# jobid, deps, flag, software, ref, input, inputSize
text=`cat $smartSlurmLogDir/allJobs.txt`

IFS=$' ';

# directly get id, deps, software, ref, and input here, if input is none, directly skip this job
output=`echo $text | awk '{if ($2 ~ /'"$SLURM_JOBID/"') print $1, $2, $3, $4, $5, $6;}'`

echo -e "Jobs on the same dependency level with current job:\n$output"
[ -z "$output" ] && { echo -e "Downstream job ids not found for $SLURM_JOBID"; summarizeRunReleaseHoldings.sh; exit; }

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
        job=`echo $text | awk '{if ($1 ~ /'"$j/"') print $3;}'`
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

        IFS=' ' read -r inputSize mem min extraMem <<< `estimateResource.sh $program ${ref//\//-} $inputs $flag 0 0 adjust`

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

#rm -r $smartSlurmLogDir/downsteamjob.adjusting 2>/dev/null

summarizeRunReleaseHoldings.sh 
