#!/bin/sh

#set -x

usage() {
    echo "Usage: debug.sh [Optional log folder name, such log, or smartSlurmLog]
          This script can guild user to debug smartSlurm runs." 
    exit 1; 
}

if [ -z "$1" ]; then 

    if [ -f ~/.smartSlurm/config/config.txt ]; then
        source ~/.smartSlurm/config/config.txt
    else
        source $(dirname $0)/../config/config.txt || { echoerr Config list file not found: config.txt; exit 1; }
    fi
    # [[ "$smartSlurmLogDir" == /* ]] || export smartSlurmLogDir=$PWD/$smartSlurmLogDir
elif [[ $1 == "-h" ||  $i == "--help" ]]; then
    usage; 
else 
    smartSlurmLogDir=$1
fi

logFolders=`ls -d $smartSlurmLogDir $smartSlurmLogDir.20* 2>/dev/null`

[ -z "$logFolders" ] && echo Log folder name does not exist: $smartSlurmLogDir. Maybe you are not in the project folder? && usage 

while true; do
    echo Availale logFolders:
    count=1
    for i in $logFolders; do
        echo "$count $i"
        count=$((count + 1))
    done
    if [ $count -eq 2 ]; then
        [[ "$xx" == "q" ]] && xx="" && break
        echo Only one logFolder. No need to select.
        x=1
    else
        echo -e "Please select the logFolder you want to check (type q to quit):"
        read -p "" x </dev/tty

        [[ "$x" == q ]] && break;

        [[ "$x" =~ ^[0-9]+$ && "$x" -lt $count && "$x" -ne 0 ]] || { echo "Out of range. Should be between > 0 and < $count"; continue; }
    fi

    logFolder=`echo $logFolders | cut -d' ' -f$x`

    failedJobs=`ls -l $logFolder/*.*.*.*.failed` 

    while true; do

        echo Available failedJobs:
        count=1
        IFS=$'\n'
        for i in $failedJobs; do
            echo $count ${i#*\/}
            count=$((count + 1))
        done

        if [ $count -eq 2 ]; then
            [[ "$xxx" == q ]] && xxx="" && break
            echo Only one failed step. No need to select.
            xx=1
        else
            echo -e "Please select a failed step (type q to quit):"; read -p "" xx </dev/tty;
            [[ "$xx" == q ]] && break;

            [[ "$xx" =~ ^[0-9]+$ && "$xx" -lt $count && "$xx" -ne 0 ]] || { echo "Out of range. Should be > 0 and < $count";  continue; }
        fi
        failedJob=`echo -e "$failedJobs" | head -n $xx | tail -n1 | awk '{print $NF}'`
        
        echo less ${failedJob/.failed/.out} >&2
        less ${failedJob/.failed/.out}
        echo less ${failedJob/.failed/.sh} >&2
    done
done
