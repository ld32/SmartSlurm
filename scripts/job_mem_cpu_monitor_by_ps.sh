#!/bin/bash

#set -x

# if jobs has --mem
reservedMem=$SLURM_MEM_PER_NODE

# if job has --mem-per-cpu and -c
[ -z "$reservedMem" ] && reservedMem=$((SLURM_MEM_PER_CPU * SLURM_JOB_CPUS_PER_NODE))

# if job has --mem-per-cpu and -n
[ -z "$reservedMem" ] &&  reservedMem=$((SLURM_MEM_PER_CPU * SLURM_CPUS_PER_TASK))

reservedCpu=$SLURM_JOB_CPUS_PER_NODE

[ -z "$reservedCpu" ] && reservedCpu=SLURM_CPUS_PER_TASK

echo "Timestamp ReservedRam(M) Ram(M) Wasted(M) Ram_Utilization(%) ReservedCpu CPU_Utilization(%)" > $SLURM_JOBID.memCpuLog

function calculate_resource_usage {
    local pid=$1
    local total_memory=0
    local total_cpu=0

    # Get memory and CPU usage of the current process
    local process_info=$(ps -o rss=,%cpu=,cmd= -p $pid)
    [[ "$process_info" == *job_mem_cpu_monitor* ]] && echo 5000 2.0 && return

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

while true; do
    total_memory_usage=0; total_cpu_usage=0
    job_pids=`ps -AF | grep $SLURM_JOBID | grep slurmstepd|awk '{print $2}'`
    if [ -n "$job_pids" ]; then
        for job_pid in $job_pids; do
            output=$(calculate_resource_usage $job_pid)
            mem=${output% *}
            cpu=${output#* }
            total_memory_usage=$((total_memory_usage + mem))
            total_cpu_usage=$(echo "scale=4; $total_cpu_usage + $cpu"|bc)
        done

        total_memory_usage=$((total_memory_usage/1024))

        echo $(date +"%Y-%m-%d %H:%M:%S") $reservedMem $(echo "scale=4; $total_memory_usage/1024" |bc) $(echo "scale=4; $reservedMem - $total_memory_usage/1024" | bc) $(echo "scale=4; $total_memory_usage/10.24/ $reservedMem" | bc) $reservedCpu $total_cpu_usage >> $SLURM_JOBID.memCpuLog
    else

        # time vs. memory for current job
        gnuplot -e "set key outside; set key reverse; set key invert; set datafile separator ' '; set style data histogram; set style histogram rowstacked gap 2; set style fill solid border rgb 'black'; set xtics rotate by -45; set terminal png size 800,600; set output 'job_$SLURM_JOBID.mem.png'; set title 'Time vs. Mem for job $SLURM_JOBID'; set xlabel 'Time (Mins)'; set ylabel 'Mem (M)'; plot 'job_$SLURM_JOBID.memCPU.txt' using 4:xtic(1) title 'Used' lc rgb 'green', '' using 5:xtic(1) title 'Wasted' lc rgb 'red'"

        # time vs. CPU usage for current job
        gnuplot -e "set key outside; set key reverse; set key invert; set datafile separator ' '; set style data histogram; set style histogram rowstacked gap 2; set style fill solid border rgb 'black'; set xtics rotate by -45; set terminal png size 800,600; set output 'job_$SLURM_JOBID.cpu.png'; set title 'Time vs. CPU Usage for job $SLURM_JOBID'; set xlabel 'Time (Mins)'; set ylabel 'CPU Usage (%)'; plot 'job_$SLURM_JOBID.memCPU.txt' using 8:xtic(1) title 'Used' lc rgb 'green'"
        exit
    fi
    sleep 20
done

