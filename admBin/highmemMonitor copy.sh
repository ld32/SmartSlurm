#!/bin/bash

set -x 
set -e

# hihhmemMonitor.sh - A script to monitor high memory usage processes
# Usage: ./hihhmemMonitor.sh [threshold]


function toS ()  {
    set +x 
    time=$1;day=0;hour=0;min=0;sec=0
    if [[ "$time" == *-* ]]; then
        day=${time%-*};tem=${time#*-};hour=${tem%%:*};min=${tem#*:};min=${min%%:*};sec=${tem#$hour:$min};sec=${sec#:}
    else
        if [[ "$time" == *\.* ]]; then
            min=${time%:*};tem=${time#*:};sec=${tem%.*}
        else
            if [[ "$time" =~ ^[0-9]+$ ]]; then
                min=$time;
            else
                sec=${time##*:};min=${time%:*};min=${min##*:};hour=${time%$min:$sec};hour=${hour%:}
            fi
        fi
    fi
    seconds=$((10#$day * 24 * 60 * 60 + 10#$hour * 60 * 60  + 10#$min * 60 + 10#$sec ))
    echo $seconds
}
date=$(date +%Y-%m-%dT)

#users=$(getent group hpc_highmem_users | cut -d: -f4 | tr ',' ' ')
users="yt145 mdp817 ye12 ers288 npp10 mez150 doa456"

echo "Users in hihhmem_user group: $users"
for user in $users; do
  echo "Processing user: $user"
  
  if [ ! -f "highmemJobs.$user.$date.txt" ]; then
    
      # find highmem partition job with starting time less than 1 month ago
      echo "highmemJobs.$user.$date.txt not found, creating it now..."
      # check if sacct command is available

      sacct --format=user,jobid,jobidraw,JobName%22,state%14,NodeList%20,Start%20,Timelimit%14,Elapsed%14,CPUTime,TotalCPU,AllocTres%25,MaxRSS,ExitCode --partition=highmem --starttime=$(date -d '1 month ago' +%Y-%m-%dT%H:%M:%S) -u $user --units=M --parsable2 > highmemJobs.$user.$date.txt

  fi 

  if [ ! -s "highmemJobs.$user.$date.txt" ]; then
    echo "No high memory jobs found for user $user in the last month. We will remove your access to highmem partition soon."\
     | mail -s "High Memory Job Alert" "lingsheng_dong@hms.harvard.edu"
  else
    echo "High memory jobs for user $user have been logged to highmemJobs.$user.$date.txt"

    total_jobs=0
    
    # Memory reservation categories
    mem_res_very_low=0    # < 25GB
    mem_res_low=0         # 25-50GB
    mem_res_medium=0      # 50-100GB
    mem_res_high=0        # > 100GB
    
    # Memory efficiency categories
    mem_eff_very_low=0    # < 25%
    mem_eff_low=0         # 25-50%
    mem_eff_medium=0      # 50-75%
    mem_eff_high=0        # > 75%
    
    # CPU efficiency categories
    cpu_eff_very_low=0    # < 25%
    cpu_eff_low=0         # 25-50%
    cpu_eff_medium=0      # 50-75%
    cpu_eff_high=0        # > 75%
    
    declare -A job_max_rss job_total_cpu job_alloc_mem job_alloc_cpu job_state job_name
    
    while IFS="|" read u jobid jobidraw JobName state nodelist start timelimit elapsed cputime totalcpu alloctres maxrss xcode; do

      total_jobs=$((total_jobs +1))
      [ $total_jobs -le 2 ] && continue # first two rows are headers 
      
      # Skip if no job ID (header or empty line)
      [[ -z "$jobid" ]] && continue
      
      # Extract base job ID (remove step suffix like .batch, .0, etc.)
      base_jobid=${jobid%%.*}
      
      # Skip if job ran less than 5 minutes
      elapsed_sec=$(toS "$elapsed")
      [ $elapsed_sec -lt 300 ] && continue
      
      # Parse AllocTres to get allocated resources
      cpu=$(echo $alloctres | grep -oP 'cpu=\K[0-9]+' || echo "0")
      mem=$(echo $alloctres | grep -oP 'mem=\K[0-9.]+' || echo "0")
      
      # Convert MaxRSS to MB (remove M suffix if present)
      maxrss_mb=$(echo $maxrss | sed 's/M$//' | grep -oP '^[0-9.]+' || echo "0")
      
      # Convert TotalCPU to seconds  
      totalcpu_sec=$(toS "$totalcpu")
      
      # Track maximum memory usage and total CPU for each base job ID
      if (( $(echo "$maxrss_mb > ${job_max_rss[$base_jobid]:-0}" | bc -l) )); then
        job_max_rss[$base_jobid]=$maxrss_mb
      fi
      
      if (( $(echo "$totalcpu_sec > ${job_total_cpu[$base_jobid]:-0}" | bc -l) )); then
        job_total_cpu[$base_jobid]=$totalcpu_sec
      fi
      
      # Store job metadata (only once per job)
      if [[ -z "${job_alloc_mem[$base_jobid]}" ]]; then
        job_alloc_mem[$base_jobid]=$mem
        job_alloc_cpu[$base_jobid]=$cpu
        job_state[$base_jobid]=$state
        job_name[$base_jobid]=$JobName
      fi
      
    done < highmemJobs.$user.$date.txt
    
    # Now process each unique job
    for base_jobid in "${!job_max_rss[@]}"; do
      req_mem_mb=${job_alloc_mem[$base_jobid]}
      used_mem_mb=${job_max_rss[$base_jobid]}
      totalcpu_sec=${job_total_cpu[$base_jobid]}
      cpu=${job_alloc_cpu[$base_jobid]}
      
      # Convert to GB for easier reading
      req_mem_gb=$(echo "scale=2; $req_mem_mb / 1024" | bc)
      used_mem_gb=$(echo "scale=2; $used_mem_mb / 1024" | bc)
      
      # Calculate memory efficiency
      if (( $(echo "$req_mem_gb > 0" | bc -l) )); then
        mem_efficiency=$(echo "scale=2; $used_mem_gb * 100 / $req_mem_gb" | bc)
      else
        mem_efficiency=0
      fi
      
      # Cap at 100%
      if (( $(echo "$mem_efficiency > 100" | bc -l) )); then
        mem_efficiency=100
      fi
      
      # Calculate CPU efficiency (this is tricky - need walltime * cores)
      # For now, let's use a simpler approach based on allocated vs used CPU time
      if [[ $cpu -gt 0 ]] && [[ $totalcpu_sec -gt 0 ]]; then
        # This is a simplified CPU efficiency - you might need to adjust based on actual walltime
        cpueffic=$(echo "scale=2; $totalcpu_sec / ($cpu * $elapsed_sec) * 100" | bc 2>/dev/null || echo "0")
      else
        cpueffic=0
      fi
      
      echo "Job $base_jobid: ReqMem=${req_mem_gb}GB, UsedMem=${used_mem_gb}GB, MemEff=${mem_efficiency}%, CPUEff=${cpueffic}%"
      
      # Count jobs by memory reservation categories
      if (( $(echo "$req_mem_gb < 25" | bc -l) )); then
        mem_res_very_low=$((mem_res_very_low + 1))
      elif (( $(echo "$req_mem_gb < 50" | bc -l) )); then
        mem_res_low=$((mem_res_low + 1))
      elif (( $(echo "$req_mem_gb < 100" | bc -l) )); then
        mem_res_medium=$((mem_res_medium + 1))
      else
        mem_res_high=$((mem_res_high + 1))
      fi
      
      # Memory efficiency categories
      if (( $(echo "$mem_efficiency < 25" | bc -l) )); then
        mem_eff_very_low=$((mem_eff_very_low + 1))
      elif (( $(echo "$mem_efficiency < 50" | bc -l) )); then
        mem_eff_low=$((mem_eff_low + 1))
      elif (( $(echo "$mem_efficiency < 75" | bc -l) )); then
        mem_eff_medium=$((mem_eff_medium + 1))
      else
        mem_eff_high=$((mem_eff_high + 1))
      fi
      
      # CPU efficiency categories
      if (( $(echo "$cpueffic < 25" | bc -l) )); then
        cpu_eff_very_low=$((cpu_eff_very_low + 1))
      elif (( $(echo "$cpueffic < 50" | bc -l) )); then
        cpu_eff_low=$((cpu_eff_low + 1))
      elif (( $(echo "$cpueffic < 75" | bc -l) )); then
        cpu_eff_medium=$((cpu_eff_medium + 1))
      else
        cpu_eff_high=$((cpu_eff_high + 1))
      fi
    done
    
    # Update total_jobs to reflect unique jobs, not job steps
    total_jobs=${#job_max_rss[@]}
  fi 

  if [[ $total_jobs -gt 0 ]]; then

    # Calculate percentages for memory reservation categories
    mem_res_very_low_pct=$(echo "scale=1; $mem_res_very_low * 100 / $total_jobs" | bc)
    mem_res_low_pct=$(echo "scale=1; $mem_res_low * 100 / $total_jobs" | bc)
    mem_res_medium_pct=$(echo "scale=1; $mem_res_medium * 100 / $total_jobs" | bc)
    mem_res_high_pct=$(echo "scale=1; $mem_res_high * 100 / $total_jobs" | bc)
    
    # Calculate percentages for memory efficiency categories
    mem_very_low_pct=$(echo "scale=1; $mem_eff_very_low * 100 / $total_jobs" | bc)
    mem_low_pct=$(echo "scale=1; $mem_eff_low * 100 / $total_jobs" | bc)
    mem_medium_pct=$(echo "scale=1; $mem_eff_medium * 100 / $total_jobs" | bc)
    mem_high_pct=$(echo "scale=1; $mem_eff_high * 100 / $total_jobs" | bc)
    
    # Calculate percentages for CPU efficiency categories
    cpu_very_low_pct=$(echo "scale=1; $cpu_eff_very_low * 100 / $total_jobs" | bc)
    cpu_low_pct=$(echo "scale=1; $cpu_eff_low * 100 / $total_jobs" | bc)
    cpu_medium_pct=$(echo "scale=1; $cpu_eff_medium * 100 / $total_jobs" | bc)
    cpu_high_pct=$(echo "scale=1; $cpu_eff_high * 100 / $total_jobs" | bc)

    # Create detailed report

    # make this report one row and save to highememReport.all.users.csv
    report="User: $user
Total highmem jobs: $total_jobs

Memory Reservation Breakdown:
- Very Low (<25GB): $mem_res_very_low jobs (${mem_res_very_low_pct}%)
- Low (25-50GB): $mem_res_low jobs (${mem_res_low_pct}%)
- Medium (50-100GB): $mem_res_medium jobs (${mem_res_medium_pct}%)
- High (>100GB): $mem_res_high jobs (${mem_res_high_pct}%)

Memory Efficiency Breakdown:
- Very Low (<25%): $mem_eff_very_low jobs (${mem_very_low_pct}%)
- Low (25-50%): $mem_eff_low jobs (${mem_low_pct}%)
- Medium (50-75%): $mem_eff_medium jobs (${mem_medium_pct}%)
- High (>75%): $mem_eff_high jobs (${mem_high_pct}%)

CPU Efficiency Breakdown:
- Very Low (<25%): $cpu_eff_very_low jobs (${cpu_very_low_pct}%)
- Low (25-50%): $cpu_eff_low jobs (${cpu_low_pct}%)
- Medium (50-75%): $cpu_eff_medium jobs (${cpu_medium_pct}%)
- High (>75%): $cpu_eff_high jobs (${cpu_high_pct}%)
"
# Detailed Job Data:
# $(cat highmemJobs.$user.$date.txt)"

    echo "$report" > highmemReport.$user.$date.txt
    
    # Calculate combined low efficiency percentages (very low + low)
    mem_res_combined_low=$(echo "scale=1; $mem_res_very_low_pct + $mem_res_low_pct" | bc)
    cpu_combined_low=$(echo "scale=1; $cpu_very_low_pct + $cpu_low_pct" | bc)
    mem_eff_combined_low=$(echo "scale=1; $mem_very_low_pct + $mem_low_pct" | bc)
    
    # Convert to integers for comparison (remove decimal part)
    mem_res_combined_int=${mem_res_combined_low%.*}
    cpu_combined_int=${cpu_combined_low%.*}
    mem_eff_combined_int=${mem_eff_combined_low%.*}
    
    if [ $mem_res_combined_int -gt 50 ] || [ $cpu_combined_int -gt 50 ] || [ $mem_eff_combined_int -gt 50 ]; then
      echo "Low Resource Efficiency Alert for user $user" >> highmemReport.$user.$date.txt
      echo "$report" | mail -s "High Memory Job Alert - $user" "lingsheng_dong@hms.harvard.edu"
    fi
  fi
  echo "Report generated for user $user: highmemReport.$user.$date.txt"
  
  #break 
done

# please order all users and make sure low efficiency users are no the top of the list

echo "High memory monitoring completed."
