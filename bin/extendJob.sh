#!/bin/sh

#set -x
set -eE

usage() { echo -e "Usage: \n${0##*/} <job ID> [Time to extend: day-hour:minute:second. Optional. Default: orignal time]"; exit 1; }

jid=$1
ext=$2

[[ "$ext" == -* ]] && echo Time can not be negtive && usage 

[[ "$jid" =~ ^[0-9]+$ ]] || { echo Job ID looks strange: $jid; usage; }

# from time string to seconds
toS () { 
    time=$1; day=0; hour=0; min=0; sec=0
    #echo Converting $time to seconds >&2
    

    if [[ "$time" == *-* ]]; then 
        day=${time%-*}; tem=${time#*-}; hour=${tem%%:*}; min=${tem#*:}; min=${min%%:*}; sec=${tem#$hour:$min}; sec=${sec#:}
    else 
        if [[ "$time" =~ ^[0-9]+$ ]]; then 
            min=$time 
        else 
            sec=${time##*:}; min=${time%:*}; min=${min##*:}; hour=${time%$min:$sec}; hour=${hour%:}
        fi
    fi    
    seconds=$((10#$day * 24 * 60 * 60 + 10#$hour * 60 * 60  + 10#$min * 60 + 10#$sec ))
    #echo day: $day hour: $hour min: $min  second: $sec seconds: $seconds >&2
    echo $seconds
}

if [ ! -z "$ext" ]; then 
    extS=`toS $ext`
    [[ "$extS" =~ ^[0-9]+$ ]] || { echo Does not look like an time; usage; }
fi 

#squeue --states RUNNING --noheader --Format=jobid:20,username:10,partition:20,timelimit:22,timeleft:22 

sqOut=`squeue  --states RUNNING --noheader --Format=jobid:20,username:10,partition:20,timelimit:22,timeleft:22 -j $jid`

echo $sqOut

original=`echo $sqOut | cut -d' ' -f 4`

[ -z "$original" ] && echo "Job not found" && exit 

originalS=`toS $original`

[ -z "$extS" ] && extS=$originalS 

newTimeS=$((originalS + extS))

thirtyDaysS=`toS 30-0:0:0`

[ "$newTimeS" -gt "$thirtyDaysS" ] && newTimeS=$thirtyDaysS

newTime=`eval "echo $(date -ud "@$newTimeS" +'$((%s/3600/24))-%H:%M:%S')"`  

echo Orignal time: $original New time: $newTime 

echo sudo scontrol update JobId=$jid TimeLimit=$newTime 

echo "Continue? (y/)"
read -p "" xx </dev/tty;

[[ "$xx" == y ]] || exit 

sudo scontrol update JobId=$jid TimeLimit=$newTime

echo -e "\nThanks for contacting us.\nI just extended the job as you requested. Let me close this ticket for now. Please contact us again at rchelp@hms.harvard.edu.\n\nBest"
