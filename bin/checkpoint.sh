#!/bin/sh

#set -x

usage(){            #           1                    2             3              4             5
    echo "checkpoint.sh <one or more commands> <unique name> <total memory> <total time> <extra memory>"
    exit 1
}

CGROUP_PATH=/sys/fs/cgroup/memory/slurm/uid_$UID/job_$SLURM_JOB_ID

export program="$1"

flag=$2

[ -z "$5" ] && usage

totalM=$3
totalT=$4
extraM=$5

# adjusted by upsteam job or re-queue due to OOT or OOM
if [ -f log/$flag.adjust ]; then
   tText=`cat log/$flag.adjust`
   totalM=`echo $tText | cut -d' ' -f1`
   totalT=`echo $tText | cut -d' ' -f2`
   extraM=`echo $tText | cut -d' ' -f3`
fi

checkpointDir="log/$flag"

mkdir -p $checkpointDir/back

portfile=port.$SLURM_JOBID

module load gcc/6.2.0 dmtcp

dmtcp_coordinator --daemon --exit-on-last -p 0 --port-file $portfile # 1>/dev/null 2>&1

hostname=`hostname`

port=`cat $portfile`

export DMTCP_COORD_HOST=$hostname

export DMTCP_COORD_PORT=$port


export DMTCP_TMPDIR=$checkpointDir/tmp

# not sure if this file can be deleted?
rm $portfile

date

START=`date +%s`

if ls $checkpointDir/ckpt_*.dmtcp >/dev/null 2>&1; then
    echo Starting from $checkpointDir checkpoint dir size:
    du -hs $checkpointDir
    dmtcp_restart -h $DMTCP_COORD_HOST -p $DMTCP_COORD_PORT  $checkpointDir/ckpt_*.dmtcp  &
elif ls $checkpointDir/back/ckpt_*.dmtcp >/dev/null 2>&1; then
    mv $checkpointDir/back/ckpt_*.dmtcp $checkpointDir 2>/dev/null
    echo Starting from $checkpointDir/back. checkpoint dir size:
    du -hs $checkpointDir
    dmtcp_restart -h $DMTCP_COORD_HOST -p $DMTCP_COORD_PORT  $checkpointDir/ckpt_*.dmtcp  &
else
    dmtcp_launch --ckptdir $checkpointDir -h $DMTCP_COORD_HOST -p $DMTCP_COORD_PORT  $program &
fi

# need check if resume or sart successful here. If not, requeeu job with more memory
sleep 10
status=`dmtcp_command -h $DMTCP_COORD_HOST -p $DMTCP_COORD_PORT -s`
if [ $? = 0 ] && [[ "$status" == *"RUNNING=yes"* ]]; then
    echo Still running ...
    rm $checkpointDir/back/ckpt_*.dmtcp 2>/dev/null
    mv $checkpointDir/ckpt_*.dmtcp $checkpointDir/back 2>/dev/null
else
    echo Something is wrong here. likely out of memory
    touch log/$flag.likelyCheckpointOOM
    exit 1
fi


checkpointed=""
while true; do
    timeR=$(($(date +%s) - START))
    echo $timeR
    if [ -z "$checkpointed" ]; then
        CURRENT_MEMORY_USAGE=$(grep 'total_rss ' $CGROUP_PATH/memory.stat | cut -d' ' -f2)
        CURRENT_MEMORY_USAGE=$(( CURRENT_MEMORY_USAGE /1024 /1024 ))
        if [ $timeR -gt 60 ]; then # && [ $(( totalM - CURRENT_MEMORY_USAGE )) -gt $checkpointM ]; then
#            echo already run for two minutes and have plenty of memory to do checkpoint

            needCheckpointM=""; needCheckpointT=""
            if [ $(( totalM - CURRENT_MEMORY_USAGE )) -lt $extraM ]; then
                echo mem low
                date
                needCheckpoint="y"
                touch log/$flag.likelyCheckpointOOM
            fi

            if [ $(( totalT * 60 - $timeR )) -lt 60 ]; then
                echo Run out of time...
                date
                needCheckpoint="y"
            fi

            if [ ! -z "$needCheckpoint" ]; then

                # make a mark on memmory plot to see memory usage when checking point
                echo 0 0 0 0 >> log/job_$SLURM_JOB_ID.mem.txt

                echo Doing checkpoint
                echo Running time in minutes: $((timeR/1000)) # this is to refer the memory usage plot to see how much memory is used by dmtcp
                status=`dmtcp_command -h $DMTCP_COORD_HOST -p $DMTCP_COORD_PORT -s`
                if [[ $? == 0 ]] && [[ "$status" == *"RUNNING=yes"* ]]; then

                    echo checkpoint folder size:
                    du -hs $checkpointDir

                    status=`dmtcp_command --ckptdir $checkpointDir -h $DMTCP_COORD_HOST -p $DMTCP_COORD_PORT --ckpt-open-files --bcheckpoint`
                    if [[ $? != 0 ]]; then
                        #[ -f log/$flag.failed ] && exit 1
                        touch log/$flag.likelyCheckpointOOM
                        exit 1
                    else
                        checkpointed=y
                        echo after checkpoint folder size:
                        du -hs $checkpointDir
                        rm dmtcp_restart_script*
                    fi
                else
                    [ -f log/$flag.success ] && { rm log/$flag.adjust 2>/dev/null || : ; exit; } || exit 1
                fi
            fi
        fi
    fi
    sleep 10
    status=`dmtcp_command -h $DMTCP_COORD_HOST -p $DMTCP_COORD_PORT -s`
    if [ $? -ne 0 ] || [[ "$status" != *"RUNNING=yes"* ]]; then
        [ -f log/$flag.success ] && { rm log/$flag.adjust 2>/dev/null || : ; exit; } || exit 1
    fi
done
