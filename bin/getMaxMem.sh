#!/bin/bash

# Set the cgroup path
CGROUP_PATH=/sys/fs/cgroup/memory/slurm/uid_$UID/job_$SLURM_JOB_ID

# Initialize the maximum memory usage to zero
MAX_MEMORY_USAGE=0

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
        echo $MAX_MEMORY_USAGE > /tmp/job_$SLURM_JOB_ID.maxMem.txt
    fi
    
    #echo $MAX_MEMORY_USAGE
    # Wait for 10 seconds before checking again
    sleep 2
done

