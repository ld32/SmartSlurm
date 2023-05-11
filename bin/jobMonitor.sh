#!/bin/bash

#set -x 

sleep 2 

# Set the cgroup path
CGROUP_PATH=/sys/fs/cgroup/memory/slurm/uid_$UID/job_$SLURM_JOB_ID

# Initialize the maximum memory usage to zero
MAX_MEMORY_USAGE=0

counter=1
counter1=0; 

jobName=$1 
originalMem=$2
originalTime=$3

defaultMem=$4

[ -f log/$jobName.adjust ] && originalMem=`cat log/$jobName.adjust | cut -d' ' -f1`
[ -f log/$jobName.adjust ] && originalTime=`cat log/$jobName.adjust | cut -d' ' -f2`

START=`date +%s`

cancelMailSent=""

# Loop indefinitely, checking the memory usage every 10 seconds
while true; do
    #[ -f /tm/job_$SLURM_JOB_ID.done ] && exit 
    # Read the current memory usage from the cgroup's memory.usage_in_bytes file
    #CURRENT_MEMORY_USAGE=$(cat $CGROUP_PATH/memory.usage_in_bytes)
    CURRENT_MEMORY_USAGE=$(grep 'total_rss ' $CGROUP_PATH/memory.stat | cut -d' ' -f2)
    #echo current: $(grep total_rss $CGROUP_PATH/memory.stat)
    # If the current memory usage is greater than the maximum seen so far, update the maximum
    if [ $CURRENT_MEMORY_USAGE -gt $MAX_MEMORY_USAGE ]; then
        MAX_MEMORY_USAGE=$CURRENT_MEMORY_USAGE
        #echo $MAX_MEMORY_USAGE > /tmp/job_$SLURM_JOB_ID.maxMem.txt
    fi
    
    ps -f -o pid,ppid,rss,vsz,args >> log/job_$SLURM_JOB_ID.ps.txt
    echo $CURRENT_MEMORY_USAG >> log/job_$SLURM_JOB_ID.ps.txt

    #echo $MAX_MEMORY_USAGE
    # Wait for 10 seconds before checking again
    sleep 2
    counter=$((counter - 1))

    if [ "$counter" -eq 0 ]; then 
        MAX_MEMORY_USAGE=$(( MAX_MEMORY_USAGE / 1024 / 1024 ))
        echo "$counter1 $MAX_MEMORY_USAGE $(($originalMem - $MAX_MEMORY_USAGE)) $((defaultMem - originalMem))" >> log/job_$SLURM_JOB_ID.mem.txt

        #if [ $originalTime -lt 120 ]; then 
            continue
        #else     
            CURRENT=`date +%s`
            min=$(( $originalTime - ($CURRENT - $START + 59)/60))

            if [[ -z "$cancelMailSent" ]] && [ $min -le 15 ]; then 
                echo "$SLURM_JOB_ID is running out of time. Please contact admin to rextend." | mail -s "$SLURM_JOB_ID is running out of time" $USER
                cancelMailSent=yes
            fi
        #fi

        counter1=$(($counter1 + 1))    
        counter=15
        MAX_MEMORY_USAGE=0 
    fi
 
done

