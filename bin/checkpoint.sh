#!/bin/sh

set -x 

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
elif [ -f log/$flag/reRun.adjust ]; then 
    tText=`cat log/$flag/reRun.adjust`
    totalM=`echo $tText | cut -d' ' -f1`
    totalT=`echo $tText | cut -d' ' -f2`
    extraM=`echo $tText | cut -d' ' -f3`    
fi



#totalM=$((totalM * 1024 * 1024))

#extraM=$((extraM * 1024 * 1024 * 2))

checkpointM=1 #$((1 * 1024 * 1024))

checkpointDir="log/$flag"

mkdir -p $checkpointDir

portfile=port.$SLURM_JOBID

dmtcp_coordinator --daemon --exit-on-last -p 0 --port-file $portfile # 1>/dev/null 2>&1 

hostname=`hostname`

port=`cat $portfile`

export DMTCP_COORD_HOST=$hostname

export DMTCP_COORD_PORT=$port

# not sure if this file can be deleted?
#rm $portfile

START=`date +%s`
    
if ls $checkpointDir/ckpt_*.dmtcp >/dev/null 2>&1; then 
    dmtcp_restart -h $DMTCP_COORD_HOST -p $DMTCP_COORD_PORT  $checkpointDir/ckpt_*.dmtcp  &

else
    dmtcp_launch --ckptdir $checkpointDir -h $DMTCP_COORD_HOST -p $DMTCP_COORD_PORT  $program &
fi 

# need check if resume or sart successful here. If not, requeeu job with more memory
sleep 5
status=`dmtcp_command -h $DMTCP_COORD_HOST -p $DMTCP_COORD_PORT -s`
if [ $? = 0 ] && [[ "$status" == *"RUNNING=yes"* ]]; then      
    echo Still running ...
else 
    echo Something is wrong here. likely out of memory
    totalM1=$totalM 
    if [ $totalM -lt 200 ]; then 
        totalM1=200
        newFactor=5
    elif [ $totalM -lt 512 ]; then 
        totalM1=512
        newFactor=4
    elif [ $totalM -lt 1024 ]; then 
        totalM1=1024
        newFactor=3
    elif [ $totalM -lt 10240 ]; then
        newFactor=2
    elif [ $totalM -lt 51200 ]; then 
        newFactor=1.5
    else 
        newFactor=1.2
    fi            
    totalM1=`echo "($totalM1*$newFactor+$extraM)/1" | bc`
    echo $totalM1 $totalT $extraM > log/$flag/reRun.adjust
    exit 1
fi 

checkpointed=""
while true; do
    if [ -z "$checkpointed" ]; then
        CURRENT_MEMORY_USAGE=$(grep 'total_rss ' $CGROUP_PATH/memory.stat | cut -d' ' -f2)
        CURRENT_MEMORY_USAGE=$(( CURRENT_MEMORY_USAGE /1024 /1024 ))
        timeR=$(($(date +%s) - START))

        if [ -z "$checkpointed" ] && [ $timeR -gt 120 ] && [ $(( totalM - CURRENT_MEMORY_USAGE )) -gt $checkpointM ]; then   
            echo already run for two minutes and have plenty of memory to do checkpoint
            
            needCheckpoint="" 
            if [ $(( totalM - CURRENT_MEMORY_USAGE )) -lt $extraM ]; then # || 
                echo mem low or runing out of time

                needCheckpoint="y"
                # in case current run fails, rerun will use this mem, time, extra mem to re-submit job
                totalM1=$totalM 
                if [ $totalM -lt 200 ]; then 
                    totalM1=200
                    newFactor=5
                elif [ $totalM -lt 512 ]; then 
                    totalM1=512
                    newFactor=4
                elif [ $totalM -lt 1024 ]; then 
                    totalM1=1024
                    newFactor=3
                elif [ $totalM -lt 10240 ]; then
                    newFactor=2
                elif [ $totalM -lt 51200 ]; then 
                    newFactor=1.5
                else 
                    newFactor=1.2
                fi            

                #newFactor=2
                totalM1=`echo "($totalM1*$newFactor+$extraM)/1" | bc`
                echo $totalM1 $totalT $extraM > log/$flag/reRun.adjust
            elif [ $(( totalT * 60 - $timeR )) -lt 300 ]; then
                needCheckpoint="y"
                echo $totalM $totalT $extraM > log/$flag/reRun.adjust
            fi
            
            if [ ! -z "$needCheckpoint" ]; then 

                status=`dmtcp_command -h $DMTCP_COORD_HOST -p $DMTCP_COORD_PORT -s`
                if [[ $? == 0 ]] && [[ "$status" == *"RUNNING=yes"* ]]; then      
                    #rm -r $checkpointDir/ckpt_*.dmtcp dmtcp_restart_script*
                    status=`dmtcp_command --ckptdir $checkpointDir -h $DMTCP_COORD_HOST -p $DMTCP_COORD_PORT --ckpt-open-files --bcheckpoint`
                    if [[ $? != 0 ]]; then 
                        totalM1=$totalM 
                        if [ $totalM -lt 200 ]; then 
                            totalM1=200
                            newFactor=5
                        elif [ $totalM -lt 512 ]; then 
                            totalM1=512
                            newFactor=4
                        elif [ $totalM -lt 1024 ]; then 
                            totalM1=1024
                            newFactor=3
                        elif [ $totalM -lt 10240 ]; then
                            newFactor=2
                        elif [ $totalM -lt 51200 ]; then 
                            newFactor=1.5
                        else 
                            newFactor=1.2
                        fi            
                        totalM1=`echo "($totalM1*$newFactor+$extraM)/1" | bc`
                        echo $totalM1 $totalT $extraM > log/$flag/reRun.adjust

                        exit 1
                    else  
                        checkpointed=y
                    fi 
                    # todo: check if it is successfule here, if not resume from earlier or rerun? 
                    
                else 
                    [ -f log/$flag/taskProcessStatus.txt ] && exit || exit 1
                fi    
            fi    
        fi    
    fi
    sleep 10
    status=`dmtcp_command -h $DMTCP_COORD_HOST -p $DMTCP_COORD_PORT -s`
    if [ $? = 0 ] && [[ "$status" == *"RUNNING=yes"* ]]; then      
        echo Still running ...
    else 
        [ -f log/$flag/taskProcessStatus.txt ] && exit || exit 1
    fi 
done

