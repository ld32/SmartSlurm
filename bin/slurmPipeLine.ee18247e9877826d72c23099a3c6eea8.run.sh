#!/bin/sh
#from cmd: /home/ld32/SmartSlurm/bin/runAsPipeline /home/ld32/SmartSlurm/bin/../scripts/bashScriptV2x.sh 12345 sbatch -A rccg -p short -c 1 --mem 4G -t 50:0 noTmp run
#run from: /home/ld32/scratch/smartSlurmTest
#set -x
echo Running $0 $@
module list
xsubd="sbatch -A rccg -p short -c 1 --mem 4G -t 50:0"
rm $smartSlurmLogDir/keepRunningJobs.txt 2>/dev/null
if [ -f $smartSlurmLogDir/allJobs.txt ]; then
    cancelAllJobsReRun.sh
    [ $? == 1 ] && exit 0;
fi
cwd=`realpath ./`
[ -f $smartSlurmLogDir/allJobs.txt ] && cp $smartSlurmLogDir/allJobs.txt $smartSlurmLogDir/allJobs.txt.old
rm -r $smartSlurmLogDir/downsteamjob.adjusting $smartSlurmLogDir/job.adjusting.lock $smartSlurmLogDir/requeue.start $smartSlurmLogDir/skipAllSuccessJobs*.txt $smartSlurmLogDir/reRunAllSuccessJobs*.txt $smartSlurmLogDir/*requeued* 2>/dev/null
[ ! -f $smartSlurmLogDir/keepRunningJobs.txt ] && printf "%-10s   %-20s   %-10s   %-10s  %-10s %-10s %-10s\n" job_id depend_on job_flag program reference inputs comment > $smartSlurmLogDir/allJobs.txt || echo -e "\n`date`" >> $smartSlurmLogDir/allJobs.txt
echo ---------------------------------------------------------
#!/bin/sh
number=$1
[ -z "$number" ] && echo -e "Error: number is missing.\nUsage: bashScriptV1.sh <numbert>" && exit 1
for i in {1..5}; do
    
    input=numbers$i.txt
    
    xsub="sbatch -p short -c 1 --mem 4G -t 50:0  -A rccg"
    #@1,0,findNumber,,,sbatch -p short -c 1 --mem 4G -t 50:0 
    echo; echo step: 1, depends on: 0, job name: findNumber, flag: 1.0.findNumber.$i   
    flag=1.0.findNumber.$i
    flag=${flag//\//_}
    deps=""
    if [ -z "$deps" ]; then
            nJobsInSameStep=$(echo "${jobIDs[1]}" | tr -cd ':' | wc -c)
            if [[ "$nJobsInSameStep" -ge 1 ]]; then
                deps="-H"
            fi
    else
        deps="-H -d afterok$deps"
    fi
    id=$(ssbatch -P findNumber -R none -F $flag -I none $deps -A rccg ${xsub#sbatch }  --wrap "findNumber.sh $number $input > $number.$i.txt;" run)
    echo -e "Got output from ssbatch: $id"
    [[ "$id" == *"missingInputFile"* ]] && echo -e "$id\n" && exit 
    [[ "$id" == *"Submitted batch job "* ]] && id=${id##*Submitted batch job } && id=${id%% *}
    [[ "$id" == *"thisJobStillRunning"* ]] && id=${id##*thisJobStillRunning } && id=${id%% *}
    if [[ "$id" == *skipThisJob* ]]; then
        jobID[1]=""
    elif [[ "$id" == *"thisJobStillRunning"* ]]; then
        id=${id##*thisJobStillRunning } && jobID[1]=${id%% *}
    elif [[ "$id" =~ ^[0-9]+$ ]]; then
        alljobs="$alljobs $id"
        startNewLoop[0]="no"
        [ -z ${startNewLoop[1]} ] && jobIDs[1]="" && startNewLoop[1]="no" 
        jobID[1]=$id
        [ -z "${jobIDs[1]}" ] && jobIDs[1]=$id || jobIDs[1]=${jobIDs[1]}:$id
    else
        echo  job $flag is not submitted
        exit 1
    fi
    xsub="sbatch -p short -c 1 --mem 4G -t 50:0  -A rccg"
    #@2,0,findNumber,,,sbatch -p short -c 1 --mem 4G -t 50:0 
    echo; echo step: 2, depends on: 0, job name: findNumber, flag: 2.0.findNumber.$i   
    flag=2.0.findNumber.$i
    flag=${flag//\//_}
    deps=""
    if [ -z "$deps" ]; then
            nJobsInSameStep=$(echo "${jobIDs[2]}" | tr -cd ':' | wc -c)
            if [[ "$nJobsInSameStep" -ge 1 ]]; then
                deps="-H"
            fi
    else
        deps="-H -d afterok$deps"
    fi
    id=$(ssbatch -P findNumber -R none -F $flag -I none $deps -A rccg ${xsub#sbatch }  --wrap "findNumber.sh $number $input > $number.$i.txt;" run)
    echo -e "Got output from ssbatch: $id"
    [[ "$id" == *"missingInputFile"* ]] && echo -e "$id\n" && exit 
    [[ "$id" == *"Submitted batch job "* ]] && id=${id##*Submitted batch job } && id=${id%% *}
    [[ "$id" == *"thisJobStillRunning"* ]] && id=${id##*thisJobStillRunning } && id=${id%% *}
    if [[ "$id" == *skipThisJob* ]]; then
        jobID[2]=""
    elif [[ "$id" == *"thisJobStillRunning"* ]]; then
        id=${id##*thisJobStillRunning } && jobID[2]=${id%% *}
    elif [[ "$id" =~ ^[0-9]+$ ]]; then
        alljobs="$alljobs $id"
        startNewLoop[0]="no"
        [ -z ${startNewLoop[2]} ] && jobIDs[2]="" && startNewLoop[2]="no" 
        jobID[2]=$id
        [ -z "${jobIDs[2]}" ] && jobIDs[2]=$id || jobIDs[2]=${jobIDs[2]}:$id
    else
        echo  job $flag is not submitted
        exit 1
    fi
done
##@2,1,mergeNumber,,,sbatch -p short -c 1 --mem 4G -t 50:0 
#cat $number.*.txt > all$number.txt
echo; echo All submitted jobs: 
awk '{print substr($0, 1, 155)}' $smartSlurmLogDir/allJobs.txt
echo ---------------------------------------------------------
