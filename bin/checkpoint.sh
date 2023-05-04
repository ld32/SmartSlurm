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
    totalT=`echo $tText | cut -d' ' -f4`
    extraM=`echo $tText | cut -d' ' -f3`
elif [ -f log/$flag/reRun.adjust ]; then 
    tText=`cat log/$flag/reRun.adjust`
    totalM=`echo $tText | cut -d' ' -f1`
    totalT=`echo $tText | cut -d' ' -f4`
    extraM=`echo $tText | cut -d' ' -f3`    
fi

totalM=$((totalM * 1024 * 1024))

extraM=$((extraM * 1024 * 1024 * 2))

checkpointM=$((10 * 1024 * 1024))

checkpointDir="log/$flag"

mkdir -p $checkpointDir

portfile=port.$SLURM_JOBID

dmtcp_coordinator --daemon --exit-on-last -p 0 --port-file $portfile # 1>/dev/null 2>&1 

hostname=`hostname`

port=`cat $portfile`

export DMTCP_COORD_HOST=$hostname

export DMTCP_COORD_PORT=$port

rm $portfile

START=`date +%s`
    
if ls $checkpointDir/ckpt_*.dmtcp >/dev/null 2>&1; then 
    dmtcp_restart -h $DMTCP_COORD_HOST -p $DMTCP_COORD_PORT  $checkpointDir/ckpt_*.dmtcp  &
else
    dmtcp_launch --ckptdir $checkpointDir -h $DMTCP_COORD_HOST -p $DMTCP_COORD_PORT  $program &
fi 

taskPID="none"

checkpointed=""

sleep 5

[ -f log/$flag/taskProcessID.txt ] && taskPID=`cat log/$flag/taskProcessID.txt` || taskPID=""
while true; do
    # todo: should move this to /tmp for less tranffic 
    #[ -f log/$flag.success ] && sleep 10 && exit
    #[ -f log/$flag.failed ] && sleep 10 && exit 1
    
    #if dmtcp_command -h $DMTCP_COORD_HOST -p $DMTCP_COORD_PORT -s; then
    #    echo Still running
    #else 
    #    exit 0 
    #fi     

    if [ -z "$checkpointed" ]; then
        CURRENT_MEMORY_USAGE=$(grep 'total_rss ' $CGROUP_PATH/memory.stat | cut -d' ' -f2)
    
        timeR=$(($(date +%s) - START))
    
        if [ $(( totalM - CURRENT_MEMORY_USAGE )) -lt $extraM ] || [ $(( totalT * 60 - $timeR )) -lt 300 ]; then
            echo mem low or runing out of time
            if [ $timeR -gt 120 ] && [ $(( totalM - CURRENT_MEMORY_USAGE )) -gt $checkpointM ]; then   
                echo already run for two minutes and have plenty of memory to do checkpoint
                
                if dmtcp_command -h $DMTCP_COORD_HOST -p $DMTCP_COORD_PORT -s; then 
                    #rm -r $checkpointDir/ckpt_*.dmtcp dmtcp_restart_script*
                    dmtcp_command --ckptdir $checkpointDir -h $DMTCP_COORD_HOST -p $DMTCP_COORD_PORT --ckpt-open-files --bcheckpoint
                    checkpointed=y

                    # in case current run fails, rerun will use this mem, time, extra mem to re-submit job
                    if [ $totalM -lt 200 ]; then 
                        totalM=200
                        newFactor=5
                    elif [ $totalM -lt 512 ]; then 
                        totalM=512
                        newFactor=4
                    elif [ $totalM -lt 1024 ]; then 
                        totalM=1024
                        newFactor=3
                    elif [ $totalM -lt 10240 ]; then 
                        newFactor=2
                    elif [ $totalM -lt 51200 ]; then 
                        newFactor=1.5
                    else 
                        newFactor=1.2
                    fi            
                
                    #newFactor=2
                    mem=`echo "($totalM*$newFactor+$extraMem)/1" | bc`
                    echo $mem $totalT $extraM > log/$flag/reRun.adjust
                else 

                    exit 0    
                fi    
                #break
            #else 
                #sacct -j $SLURM_JOBID 
                #[ -f log/$flag.fail ] && exit 1
                #[ -f log/$flag.success ] && exit || { touch log/$flag.failed && exit 1; }
            fi
        fi    
    #else
    #    sacct -j $SLURM_JOBID 
    #    [ -f log/$flag.success ] && exit || { touch log/$flag.failed && exit 1; }
    fi
    
    sleep 10
    # oot and oot happens, command under dmtcp directly killed, trap does not work. So we need manually check the
    # the task process status and tell the job step finished
    
    if [ ! -z "$taskPID" ]; then 
        if ps -p $taskPID > /dev/null; then 
            echo "$taskPID is running"
        else 
            sleep 2
            # killed by cgroup
            echo "$taskPID is finished"
            if [ -f log/$flag/taskProcessStatus.txt ]; then
                status=`cat log/$flag/taskProcessStatus.txt`
                echo process status: $status;
             #   if [[ "$status" == 0 ]]; then
                    exit 0
              #  else
              #      exit 1
               # fi        
            else 
                exit 1 # task was killed due to oom or oot
            fi     
        fi
    fi 
done

# while  [ -z "$checkpointed" ]; do         
#         CURRENT_MEMORY_USAGE=$(grep 'total_rss ' $CGROUP_PATH/memory.stat | cut -d' ' -f2)

#         timeR=$(($(date +%s) - START))

#         if [ $(( totalM - CURRENT_MEMORY_USAGE )) -lt $extraM ] || [ $(( totalT * 60 - $timeR )) -lt 300 ] && [ $timeR -gt 120 ] && [ $(( totalM - CURRENT_MEMORY_USAGE )) -gt $checkpointM ] ; then   
#             #rm -r $checkpointDir/ckpt_*.dmtcp dmtcp_restart_script*
            
#             #break  
#         fi
#      sleep 10
# done

# wait 

#     #else 
#         [ -f log/$flag.success ] && exit || { touch log/$flag.failed && exit 1; }
#     #fi        
#     #sleep 10    
# #done
