#!/bin/sh

#set -x

usage(){            #           1                    2                3             4              5   
    echo "checkpoint.sh    <one or more commands> <unique name> <total memory> <total time> <extra memory>"
    exit 1
}

echo Running $0 $@
echo pwd: `pwd`

#smartSlurmLogDir=$1

export program="$1"

flag=$2

[ -z "$5" ] && usage

totalM=$SLURM_MEM_PER_NODE #$3
totalT=$4
extraM=$5

# adjusted by upsteam job or re-queue due to OOT or OOM
if [ -f $smartSlurmLogDir/$flag.adjust ]; then
   tText=`cat $smartSlurmLogDir/$flag.adjust`
   #totalM=`echo $tText | cut -d' ' -f1`
   totalT=`echo $tText | cut -d' ' -f2`
   extraM=`echo $tText | cut -d' ' -f3`
fi

checkpointDir="$smartSlurmLogDir/$flag"

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
# elif ls $checkpointDir/back/ckpt_*.dmtcp >/dev/null 2>&1; then
#     mv $checkpointDir/back/ckpt_*.dmtcp $checkpointDir 2>/dev/null
#     echo Starting from $checkpointDir/back. checkpoint dir size:
#     du -hs $checkpointDir
#     dmtcp_restart -h $DMTCP_COORD_HOST -p $DMTCP_COORD_PORT  $checkpointDir/ckpt_*.dmtcp  &
else
    dmtcp_launch --ckptdir $checkpointDir -h $DMTCP_COORD_HOST -p $DMTCP_COORD_PORT  $program &
fi

# need check if resume or sart successful here. If not, requeeu job with more memory
sleep 20
status=`dmtcp_command -h $DMTCP_COORD_HOST -p $DMTCP_COORD_PORT -s`
if [ $? = 0 ] && [[ "$status" == *"RUNNING=yes"* ]]; then
    #set -x 
    echo Still running ...
    lastFile=$(ls -r $checkpointDir/job*memCPU.txt 2>/dev/null | tail -n 1 )
    if [ ! -z "$lastFile" ]; then 
        mv $smartSlurmLogDir/job_$SLURM_JOB_ID.memCPU.txt $smartSlurmLogDir/job_$SLURM_JOB_ID.memCPU.txt.overwrittenByPreviousCheckpoint
        cat $lastFile > $smartSlurmLogDir/job_$SLURM_JOB_ID.memCPU.txt
    fi 

    if ls $checkpointDir/ckpt_*.dmtcp 2>/dev/null; then 
        echo -e "Read checkpoint data sucessfully by job $SLURM_JOBID" | mail -s "CPRead:$SLURM_JOBID:$flag" $USER
        touch $smartSlurmLogDir/$flag.startFromCheckpoint
        echo 0 0 0 0 0  jobID:$SLURM_JOBID time:`date` resume from  checkpoint >> $smartSlurmLogDir/job_$SLURM_JOB_ID.memCPU.txt
        rm $checkpointDir/back/ckpt_*.dmtcp 2>/dev/null
        mv $checkpointDir/ckpt_*.dmtcp $checkpointDir/back 2>/dev/null
    # esif ls $checkpointDir/back/ckpt_*.dmtcp; then 
    #      echo -e "Read checkpoint data sucessfully by job $SLURM_JOBID" | mail -s "CPRead:$SLURM_JOBID:$flag" $USER
    #     touch $smartSlurmLogDir/$flag.startFromCheckpoint && echo 0 0 0 0 0  jobID:$SLURM_JOBID time:`date` resume from  checkpoint >> $smartSlurmLogDir/job_$SLURM_JOB_ID.memCPU.txt

    fi 

else
    
    
    if ls $checkpointDir/ckpt_*.dmtcp 2>/dev/null; then
        touch $smartSlurmLogDir/$flag.likelyCheckpointOOM
    else 
        echo Something is wrong here. likely out of memory
        echo -e "Checkpoint failed by job :$SLURM_JOBID" | mail -s "CPFail:$SLURM_JOBID:$flag" $USER
        touch $smartSlurmLogDir/$flag.likelyCheckpointDoNotWork
    fi    
    exit 1
fi

checkpointed=""
while true; do
    timeR=$(($(date +%s) - START))
    #echo $timeR
    if [ -z "$checkpointed" ]; then
        
        #CURRENT_MEMORY_USAGE=$(grep 'total_rss ' $CGROUP_PATH/memory.stat | cut -d' ' -f2)
        CURRENT_MEMORY_USAGE=`tail -1 $smartSlurmLogDir/job_$SLURM_JOB_ID.memCPU.txt | awk '{print $2}'`
        #CURRENT_MEMORY_USAGE=$(( CURRENT_MEMORY_USAGE /1024 ))
        if [ $timeR -gt 120 ]; then # && [ $(( totalM - CURRENT_MEMORY_USAGE )) -gt $checkpointM ]; then
#            echo already run for two minutes and have plenty of memory to do checkpoint

            if [ $(( totalM - CURRENT_MEMORY_USAGE )) -lt 100 ]; then
                echo mem low ...
                date
                needCheckpoint="Mem"
                touch $smartSlurmLogDir/$flag.likelyCheckpointOOM
            fi

            if [ $(( totalT * 60 - $timeR )) -lt 120 ]; then
                echo time low ...
                date
                needCheckpoint="Time"
            fi

            if [ ! -z "$needCheckpoint" ]; then

                # make a mark on memmory plot to see memory usage when checking point
                echo 0 0 0 0 0  jobID:$SLURM_JOBID time:`date` start checkpoint >> $smartSlurmLogDir/job_$SLURM_JOB_ID.memCPU.txt

                echo Doing checkpoint
                echo Running time in minutes: $((timeR/1000)) # this is to refer the memory usage plot to see how much memory is used by dmtcp
                status=`dmtcp_command -h $DMTCP_COORD_HOST -p $DMTCP_COORD_PORT -s`
                if [[ $? == 0 ]] && [[ "$status" == *"RUNNING=yes"* ]]; then

                    echo checkpoint folder size:
                    du -hs $checkpointDir
                    
                    #status=`dmtcp_command --ckptdir $checkpointDir -h $DMTCP_COORD_HOST -p $DMTCP_COORD_PORT --ckpt-open-files --bcheckpoint`
                    status=`dmtcp_command --bcheckpoint`
                    if [[ $? != 0 ]]; then
                        #[ -f $smartSlurmLogDir/$flag.failed ] && exit 1
                        echo -e "Wrote checkpoint data failed by job $SLURM_JOBID" | mail -s "CPFail:$SLURM_JOBID:$flag" $USER
                        touch $smartSlurmLogDir/$flag.likelyCheckpointOOM
                        #exit 1
                    else
                        echo -e "Wrote checkpoint data sucessfully by job $SLURM_JOBID" | mail -s "CPWrite.$needCheckpoint:$SLURM_JOBID:$flag" $USER
                        checkpointed=y
                        echo after checkpoint folder size:
                        du -hs $checkpointDir
                        rm dmtcp_restart_script*

                        echo 0 0 0 0 0  jobID:$SLURM_JOBID time:`date` end checkpoint >> $smartSlurmLogDir/job_$SLURM_JOB_ID.memCPU.txt

                        # copy mem and cpu usage to checkpoint folder
                        cp $smartSlurmLogDir/job_$SLURM_JOB_ID.memCPU.txt $checkpointDir

                        #touch $smartSlurmLogDir/$flag/checkpointDone$needCheckpoint

                        #cat /var/spool/slurmd/job$SLURM_JOBID/slurm_script
                        #exit 1
                    fi
                else

                    echo -e "Checkpoint failed by job $SLURM_JOBID" | mail -s "CPFail:$SLURM_JOBID:$flag" $USER
                    #[ -f $smartSlurmLogDir/$flag.success ] && { rm $smartSlurmLogDir/$flag.adjust 2>/dev/null || : ; exit; } || exit 1
                    touch $smartSlurmLogDir/$flag.likelyCheckpointOOM
                    #exit 1
                    #rm $smartSlurmLogDir/$flag.adjust 2>/dev/null || :;  exit 1
                fi
            fi
        fi
    fi
    sleep 10
    
    [ -f $smartSlurmLogDir/$flag.success ] && exit
    # status=`dmtcp_command -h $DMTCP_COORD_HOST -p $DMTCP_COORD_PORT -s 2>/dev/null`
    # if [ $? -ne 0 ] || [[ "$status" != *"RUNNING=yes"* ]]; then
    #     echo All done!
    #     #[ -f $smartSlurmLogDir/$flag.success ] && { rm $smartSlurmLogDir/$flag.adjust 2>/dev/null || : ; exit; } || exit 1
    #     #rm $smartSlurmLogDir/$flag.adjust 2>/dev/null || : ;
    #     exit
    # fi
done
