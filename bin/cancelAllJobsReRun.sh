#!/bin/sh

#set -x

IFS=$'\n'

out=`squeue -u $USER -t PD,R -o "%.18i"`

[ -z "$out" ] && exit 0;

declare -A nmap

lines=`tail -n +2 log/allJobs.txt | awk 'NF>2{print $1, $2, $3}'`
for line in $lines; do
    if [ ! -z "${line/ /}" ]; then
        id=${line%% *} #`echo $line | cut -d' ' -f1`
        if [[ "$out" == *$id* ]]; then
            ids="$ids $id"
            nmap[$id]=${line##* } #`echo $line | cut -d' ' -f3`
            notDone="$notDone$line\n" #flags="$flags $flag"
        fi
    fi
done

if [ ! -z "${ids/ /}" ]; then
    IFS=' '
    echo the following jobs are not finished yet:
    echo -e "$notDone"
    echo -e "Do you want to cancell all running jobs? " >&2
    echo -e "y:        Cancel all pending/running jobs and re-run the whole pipiline (will ask if you want to re-run successfuly finished jobs if there are any)." >&2
    echo -e "n:        Not cancel the pending/running jobs, only submit other jobs (will ask if you want to re-run successfuly finished jobs if there are any)." >&2
    echo -e "Enter:    Quit, don't anything." >&2
    read -p "" x </dev/tty

    if [[ "$x" == "y" ]]; then
        for id in $ids; do #`echo -e "$ids" | cut -d' ' -f1`; do
          scancel $id
          echo ${nmap[$id]} cancelled.
          touch log/${nmap[$id]}.user.killed
        done
    elif [[ "$x" == "n" ]]; then
        #cat log/allJobs.txt.old >> log/allJobs.txt
        echo -e  "$notDone" > log/keepRunningJobs.txt
        #echo This is still under development!
        #exit 1
    else
        exit 1
    fi
else
    echo Could not find any jobs to cancel.
fi
