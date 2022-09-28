#!/bin/sh

#set -x 

Usage="Usage: $0 full_path_to_flag_folder  flag(job name) \n  Note: this script will go through job id list file, find the downstream jobs, and return them as a string of job flags. "
echo Running: $0  $@

path=$1

[ -z "$2" ] && { echo -e "$Usage"; exit; }

[ -f $1/allJobs.txt ] || { echo -e "job id file $path/allJobs.txt does not exist\n$Usage"; exit 1; }

# jobid, deps, flag, software, ref, input, inputSize
text=`cat $path/allJobs.txt`
 
job=$2

IFS=$' ';  

# check the third column for the job name, then find the the job id in column 1
id=`echo $text | awk '{if ($3 ~ /'"$job/"') print $1;}' | tail -n 1`

echo -e "Find current job id (flag: $2):\n$id"
[ -z "$id" ] && { echo -e "job id for job name $job not found!\n$Usage"; exit; }

echo 

echo Find all downstream jobs which depend on current job
idNames=`echo $text | awk '{if ($2 ~ /'"$id/"') print $2, $3;}'`

echo -e "job idNames:\n$idNames"
[ -z "$idNames" ] && { echo -e "Downstream job ids not found for $id"; exit; }

# sleep for a random number of seconds to avoid race condition
sleep $(( ( RANDOM % 120)  + 1 ))

IFS=$'\n';
for i in $idNames; do
    echo 1working on $i
    id1=${i% *}; name=${i#* };
    allDone=""
    IFS=$' '; 
    for j in ${id1//\./ }; do 
        echo 2working on $j
        [[ "$j" == "$id" ]] && echo Ignore. It is the current job. It should adjust the mem and time for the downsteam job.  && continue
        
        echo look for the job flag for $j
        job=`echo $text | awk '{if ($1 ~ /'"$j/"') print $3;}'`
        #[ -z "$job" ] && { echo -e "job name not found!"; exit; }
        [ -f "$path/$job.success" ] && echo -e This job was done! || { echo Dependent job is not done yet: $job; allDone=no; }         
    done
    if [ -z "$allDone" ]; then
        echo Dependants for $name are all done except for the current job. Ready to adjust mem/runtime
        
        # look for the downstream job info:                                     jobID,software, ref, inputs
        output=`cat $path/allJobs.txt| awk '{if ($3 ~ /'"$name/"') print $1, $4, $5, $6;}'`
    
        id2=${output%% *}; inputs=${output##* }; 
        softwareRef=${output#* };  softwareRef=${softwareRef% *}
        justRunStats=""
        mem="" 
        if [ ! -f ~/.smartSlurm/${softwareRef/ /.}.mem.stat.final ]; then    
            echo Do not have a formula. Let us build one...
            jobStatistics.sh $softwareRef 4
            justRunStats=yes
        fi
        if [  -f ~/.smartSlurm/${softwareRef/ /.}.mem.stat.final ]; then
            echo estimating here
            inputSize=`{ du --apparent-size -c -L ${inputs//,/ } 2>/dev/null || echo notExist; } | tail -n 1 | cut -f 1`

            if [[ "$inputSize" == "notExist" ]]; then
                echo Some or all input files not exist: $inputs  && exit 1
            else
                $iputSize=$(($inputSize/1024)) # convert to M
                echo inputSize: $inputSize
                echo running: estimateMemTime.sh  $softwareRef $inputSize
                output=`estimateMemTime.sh $softwareRef $inputSize`    
                
                if [[ "$output" == "outOfRange" ]]; then 
                    echo Input size is too big for the curve to estimate!
                    if [ -z "$justRunStats" ]; then 
                        echo Delete the curve and try to re-run the statistics.
                        rm ~/.smartSlurm/$software.$ref.mem.stat.final ~/.smartSlurm/$software.$ref.time.stat.final ~/.smartSlurm/jobRecord.txt
                        jobStatistics.sh $softwareRef 4 
                        if [ -f ~/.smartSlurm/${softwareRef/ /.}.mem.stat.final ]; then    
                            output=`estimateMemTime.sh $softwareRef $inputSize`
                            if [[ "$output" == "outOfRange" ]]; then 
                                echo Input size is too big for the curve to estimate! Use default mem and runtime to submit job.
                            else
                                mem=${output% *}M; time=${output#* };     
                                echo Got estimation inputsize: $inputSize mem: $mem  time: $time
                            fi 
                        else 
                            echo Building formula failed. Use default mem and runtime to submit job.
                        fi
                    else 
                        echo Use default mem and runtime to submit job.    
                    fi    
                else
                    mem=${output% *}M; time=${output#* };     
                    echo Got estimation inputsize: $inputSize mem: $mem  time: $time
                fi 
                [ -z "$mem" ] && continue
              
                echo Got estimation inputsize: $inputSize mem: $mem  time: $time
                hours=$((($time + 59) / 60))
    
                echo looking partition for hour: $hours
                x=`realpath $0`
                . ${x%\/bin\/adjustDownStreamJobs.sh}/config/partitions.txt || { echo Partition list file not found: partition.txt; exit 1; }
                setPartition $hours partition
                
                seconds=$(($time * 60))
                time=`eval "echo $(date -ud "@$seconds" +'$((%s/3600/24))-%H:%M:%S')"`
                scontrol show job $id2                
                echo running: scontrol update jobid=$id2 timelimit=$time partition=$partition MinMemoryNode=${mem}
                scontrol update JobId=$id2 TimeLimit=$time Partition=$partition  MinMemoryNode=${mem}
                scontrol show job $id2

                echo "#scontrol update JobId=$id2 TimeLimit=$time Partition=$partition  MinMemoryNode=${mem}" >> $path/$name.sh

                memN=$(($mem-5))M; timeN=$(($seconds-300))
                timeN=`eval "echo $(date -ud "@$timeN" +'$((%s/3600/24))-%H:%M:%S')"`

                #old srun command:
                #srun -n 1 -t 0-03:59:55 --mem 40955M bash  
                echo Original slurm script: 
                cat $path/$name.sh 
                sed -i "s/srun -n 1 -t [0-9]*-[0-9]*:[0-9]*:[0-9]* --mem [0-9]*M/srun -n 1 -t $timeN --mem $memN/g" $path/$name.sh 

                # todo: also need modify the noExist input parameter for emmai and record here 
                sed -i "s/notExist/$inputSize/" $path/$name.sh 
                
                echo New slurm script:
                cat $path/$name.sh 
                               
            fi
        else 
            echo Try to build fomular, but it was not successful
        fi
    else 
        echo Need wait for other jobs to finish before we can ajust
    fi   
        
done
    
