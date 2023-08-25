#!/bin/bash

set -x

echo Running: $0 $@
echo pwd: `pwd`
cd ${1%/log/*}
echo pwd `pwd`

# requeue failed jobs
if ls log/*.requeueCMD && mkdir log/requeue.start; then
    for requeue in log/*.requeueCMD; do
        if grep -q $SLURM_JOBID $requeue; then
            rm $requeue
            rm ${requeue%.requeueCMD}.failed
        else
           sh $requeue && rm $requeue && rm ${requeue%.requeueCMD}.failed
        fi
        sleep 1
    done
    rm -r log/requeue.start
fi
set +x

sleep 2

counter=0

jobName=${1#*/log/}
originalMem=$2
originalTime=$3

defaultMem=$4

[ -f log/$jobName.adjust ] && originalMem=`cat log/$jobName.adjust | cut -d' ' -f1`
[ -f log/$jobName.adjust ] && originalTime=`cat log/$jobName.adjust | cut -d' ' -f2`

cancelMailSent=""

function calculate_resource_usage {
    local pid=$1
    local total_memory=0
    local total_cpu=0

    # Get memory and CPU usage of the current process
    local process_info=$(ps -o rss=,%cpu= -p $pid)
    local memory=$(echo "$process_info" | awk '{print $1}')
    local cpu=$(echo "$process_info" | awk '{print $2}')

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

# Loop indefinitely, checking the memory usage every 10 seconds
while true; do
    job_pid=`ps -AF|grep $SLURM_JOBID|grep slurmstepd|awk '{print $2}'|tail -1`
    if [ -n "$job_pid" ]; then
        output=$(calculate_resource_usage $job_pid)
        total_memory_usage=${output% *}
        total_cpu_usage=${output#* }
        #echo $(date +"%Y-%m-%d %H:%M:%S") $reservedMem $(echo "scale=4; $total_memory_usage/1024" | bc) $(echo "scale=4; $total_memory_usage/10.24/ $reservedMem" | bc) $reservedCpu $total_cpu_usage >> $SLURM_JOBID.memCpuLog
        counter=$((counter + 1))
        saved=$((defaultMem - originalMem))
        [ "$saved" -lt 0 ] && saved=0

        total_memory_usage=$((total_memory_usage/1024))

	    echo "p$counter $total_memory_usage $(($originalMem - $total_memory_usage)) $saved ${total_cpu_usage%.*}" >> log/job_$SLURM_JOB_ID.memCPU.txt

    else
        exit
    fi

   	if [ -z "$cancelMailSent" ] && [ "$originalTime" -ge 120 ]; then
        CURRENT=`date +%s`
        min=$(( $originalTime - ($CURRENT - $START + 59)/60))
        if [ $min -le 15 ]; then
            echo "$SLURM_JOB_ID is running out of time. Please contact admin to rextend." | mail -s "$SLURM_JOB_ID is running out of time" $USER
            cancelMailSent=yes
        fi
    fi
    sleep 5
done
