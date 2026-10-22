#!/bin/bash

#set -x

output="Running: $0"

counter=0

jobName="$SLURM_JOB_NAME"
[ -z "$jobName" ] && jobName=jobWithoutName

# # if jobs has --mem
reservedMem=$SLURM_MEM_PER_NODE

# # if job has --mem-per-cpu and -c
[ -z "$reservedMem" ] && reservedMem=$((SLURM_MEM_PER_CPU * SLURM_JOB_CPUS_PER_NODE))

# # if job has --mem-per-cpu and -n
[ -z "$reservedMem" ] &&  reservedMem=$((SLURM_MEM_PER_CPU * SLURM_CPUS_PER_TASK))

core="$SLURM_CPUS_ON_NODE"

[ -z "$core" ] && core="$SLURM_CPUS_PER_TASK"

function calculate_resource_usage {
    local pid=$1
    local total_memory=0
    local total_cpu=0

    # Get memory and CPU usage of the current process
    local process_info=$(ps -o rss=,%cpu=,pid=,ppid=,cmd= -p $pid)
    #echo $process_info >&2
    [[ "$process_info" == *memCpuMonitor* ]] && echo 5000 2.0 && return

    local memory=$(echo "$process_info" | awk '{print $1}')
    local cpu=$(echo "$process_info" | awk '{print $2}')
    [ -z "$cpu" ] && cpu=0

    # Calculate memory and CPU usage of children
    local children=$(ps --ppid $pid -o pid=)
    for child_pid in $children; do
        local output=$(calculate_resource_usage $child_pid)
        local child_memory=${output% *}
        local child_cpu=${output#* }
        total_memory=$((total_memory + child_memory))
        total_cpu=$(echo "scale=4; $total_cpu + $child_cpu"|bc)
    done

    # Calculate total memory and CPU usage
    total_memory=$((total_memory + memory))
    total_cpu=$(echo "scale=4; $total_cpu + $cpu"|bc)
    
    echo "$total_memory $total_cpu"
}


START=`date +%s`
echo "count memUsed memReserve cpuUsed cpuReserved $START jobID:$SLURM_JOBID time:`date`" >> job_$SLURM_JOB_ID.memCPU.txt

# Loop indefinitely, checking the memory usage every 10 seconds
while true; do
    total_memory_usage=0; total_cpu_usage=0
    job_pids=`ps -AF | grep $SLURM_JOBID | grep slurmstepd|awk '{print $2}'`
    #echo job pids: $job_pids  >&2
    if [ -n "$job_pids" ]; then
        for job_pid in $job_pids; do
            #echo working on: $job_pid
            output=$(calculate_resource_usage $job_pid)
            mem=${output% *}
            cpu=${output#* }
            total_memory_usage=$((total_memory_usage + mem))
            total_cpu_usage=$(echo "scale=4; $total_cpu_usage + $cpu"|bc)
            #echo Total: $total_memory_usage, $total_cpu_usage
        done

        counter=$((counter + 1))
        saved=$((defaultMem - reservedMem))

        total_memory_usage=$((total_memory_usage/1024))

	    echo "$counter $total_memory_usage $$reservedMem ${total_cpu_usage%.*} $core" >> job_$SLURM_JOB_ID.memCPU.txt

    else
        exit
    fi
    

    #todo: need a way to auto quit this background script
    sleep 60
done

