#!/bin/bash

## notes: to estimate resource, need at least three jobs if no input size is given, and 5 jobs if there is input size given

## features for regular sbatch (when calling without -L)
# Auto adjust partition according to run-time request if they does not match up
# Auto check if slurm script exists
# Auto create output and erro folders if not exist
# Auto adjust memory and run time based on earlier mem time usage
# Support re-run checking (using output as checkpoint)
# Auto adjust memory and run time based on earlier mem time usage
# Auto requeue job with double time when out of time
# Auto requeue job with double memory when out of memory
# Keep good log and send informative email
# allow user level config file in ~/.smartSlurm/config/config.txt

# Limitations:
# does not take logs in logDir, does not take care of re-run because there is not unique flag for the job, does not take care of downstream jobs

# todo:
# Forjob with breaking-point. Need more work to adjust how job record the run time and memory to jobRecords.txt.

# the logic:
# If job does not depends on other job, and has input, estimate mem/timed and submit
# If a job has no inputs, directly submit
# Adjust mem/time from upstream job for jobs with input. todo list: adjust jobs without input as well?

#set -e
#set -x
#set -u

parentCmd=`ps -h --pid $(ps -h --pid $$ -o ppid) -o command`

# nextflow
if [ -f .command.sh ] && [ -f .command.run ]; then 
    echoerr() { echo "$@" >> ../../../.nextflow.log; }
# snakemake
else #if [[ "$parentCmd" == */bin/snakemake* ]]; then
    echoerr() { echo "$@" >> .smartSlurm.log; }
fi 

echoerr1() { echo "$@" >&2; }

#echoerr parent: $parentCmd 

usage() { echo -e "Short Usage: \n${0##*/} [-P progam] [-R reference] [-F uniqueJobFlag (for smartPieline)s] [-I inputList] [-d Dependencies (for smartPipeline)] [sbatch options] [run]\nDetail Detail Usage:\n${0##*/} [-P progam, optional. Such as: bowtie2-4core. If empty, use script name as program name] [-R reference, optional. Such as: hg19] [-F uniqueJobFlag, optional. Such as 1.bowtie.s1.fa] [-I inputFileOrFolderList, optional. Such as: read1.fq,read2.fq] [ -d dependent jobs] <regular sbatch options, optional. Such as: job.sh or -p short -c 1 -t 2:0:0 --mem 2G --wrap \"my_application para1 para2\"> [run, optional: run will submit job, empty will do a dry run without submitting a job.]"; exit 1; }

[ -z "$1" ] || [[ "-h" == "$1" ]] || [[ "--help" == "$1" ]] && usage

if [ -z "$smartSlurmJobRecordDir" ]; then
    if [ -f ~/.smartSlurm/config/config.txt ]; then
        source ~/.smartSlurm/config/config.txt
    else
        source $(dirname $0)/../config/config.txt || { echo Config list file not found: config.txt; exit 1; }
    fi
fi

[[ "$smartSlurmLogDir" == /* ]] || export smartSlurmLogDir=$PWD/$smartSlurmLogDir

mkdir -p $smartSlurmJobRecordDir/stats

[ -f  $smartSlurmJobRecordDir/jobRecord.txt ] || echo 1jobID,2inputSize,3memDefault,4timeDefaut,5memAllocated,6timeAllocated,7memUsed,8timeUsed,9JobStatus,10useID,11saccMem,12program,13reference,14flag,15core,16extraMem,17extraTime,18date > $smartSlurmJobRecordDir/jobRecord.txt

current_date=$(date '+%Y-%m-%d_%H-%M-%S_%4N') 
echoerr ssbatch run date: $current_date
echoerr Running: $0 "$@"
echoerr pwd: `pwd`

cmd="ssbatch"
whitespace="[[:space:]]"
for i in "$@"; do
    if [[ $i =~ $whitespace ]]; then
        i="\"$i\""
    fi
    cmd="$cmd  $i"
done
echoerr $cmd

testRun=${@: -1}

#unset P R I d C

echoerr
array=( "$@" )
arrayC=("${array[@]}")  

export alwaysRequeueIfFail=""
caller=`ps -f $PPID` && [[ "$caller" == *slurmPipeLine* ]] && runningSingleJob="" || runningSingleJob=y

[ ! -z "$runningSingleJob" ] && set -- "$@" run && testRun=run

# nextflow pipeline, log dir is on the top folder
if [ -f .command.sh ] && [ -f .command.run ]; then 
    export smartSlurmLogDir="../../../$smartSlurmLogDir"
fi 
mkdir -p $smartSlurmLogDir/

#slurmScriptPosition=$(($(($#))-2))

# get the first few parameters for ssbatch
for (( i=0; i<$(($#)); i++ )); do
    [ -z "${array[$i]}" ] && continue
  	#echoerr $i " / " $(($#)) " : " ${array[$i]}
  	case "${array[$i]}" in
        "-H"            )   onhold=-H;;
        "-D"            )   chdir="${array[$i+1]}" && array[$i+1]="";;
        --chdir=*       )   chdir="${array[$i]}" && chdir="${chdir/--chdir=/}";;
        --comment=*     )   comment="${array[$i]}";; # && comment=${comment/--comment=/};; # && for c in "${arr[@]}"; do eval $c; done;;
  		"-P" 			)   program="${array[$i+1]}" && array[$i+1]="";;
  		--program=* 	)   program="${array[$i]}" && program="${program/--program=/}";;
  		"-R" 			)   ref="${array[$i+1]}" && array[$i+1]="";;
  		--Ref=* 		)   ref="${array[$i]}" && ref="${ref/--Ref=/}";;
  		"-F" 			)   flag="${array[$i+1]}" && array[$i+1]="";;
  		--Flag=* 	    )   flag="${array[$i]}" && flag="${flag/--Flag=/}";;
        "-I" 			)   inputs="${array[$i+1]}" && array[$i+1]="";;
  		--Inputs=* 	    )   inputs="${array[$i]}" && inputs="${inputs/--Inputs=/}";;
        "-d"            )   deps="${array[$i+1]}" && array[$i+1]="";;
        "--Dependency=" )   deps="${array[$i]}" && deps="${deps/--Dependency=/}";;
        "--alwaysRequeueIfFail") export alwaysRequeueIfFail=true;;
  	    "-p"            )   partition="${array[$i+1]}" && array[$i+1]="";; # will set partition later
  		--partition=*   )   partition="${array[$i]}" && partition=${partition/--partition=/};;
  		--mem-per-cpu=* ) 	mem1="${array[$i]}" && mem1="${mem1/--mem-per-cpu=/}";;
  		"-c" 	 		) 	core="${array[$i+1]}" && array[$i+1]="";;
  		--cpus-per-task=*)  core="${array[$i]}" && core="${core/--cpus-per-task=/}";;
  		"-n" 	 		) 	task="${array[$i+1]}" && array[$i+1]="";;
  		--ntasks=*      )   task="${array[$i]}" && task="${task/--ntasks=/}";;
  		"-N" 	 		) 	node="${array[$i+1]}" && array[$i+1]="";;
  		--nodes=*       )   node="${array[$i]}" && node="${node/--nodes=/}";;
  		"--mem"  		) 	mem="${array[$i+1]}" && array[$i+1]="";; # will set memory later
  		--mem=* 		)  	[ -z "$mem" ] && mem="${array[$i]}" && mem=${mem/--mem=/};;
  		"-t" 			)   time="${array[$i+1]}" && array[$i+1]="";;
  		--time=* 		)   time="${array[$i]}" && time="${time/--time=/}";;
  		"-J" 			)   name="${array[$i+1]}" && array[$i+1]="";;
  		--job-name=* 	)   name="${array[$i]}" && name="${time/--job-name=/}";;
  		"-o" 			)   out="${array[$i+1]}" && array[$i+1]="";;
  		--output=* 		)   out="${array[$i]}" && out="${out/--oupput=/}";;
  		"-e" 			)   err="${array[$i+1]}" && array[$i+1]="";;
  		--error=* 		)   err="${array[$i]}" && err="${err/--error=/}";;
        "-A"            )   slurmAcc="-A ${array[$i+1]}" && array[$i+1]="";;
        "-d" 			)   deps="${array[$i+1]}" && array[$i+1]="";;
  	    --dependency=* 	)   deps="${array[$i]}" && deps="${deps/--dependency=/}";;
        "test"          )   [ $i -eq $(($#)) ] && continue;;
        "run"           )   [ $i -eq $(($#)) ] && continue;;
        "--wrap"        )   wrapCMD="${array[$i+1]}" && array[$i+1]="";;
        --wrap=*        )   wrapCMD="${array[$i]}" && wrapCMD="${wrapCMD/--wrap=}";;
        *               )   { [ -z "$slurmScript" ] && [ -f "${array[$i]}" ] && [[ " -a -A -b -c -d -D -e  -i -J -L -M -m -n -N -o -p -q -P -t -F -w -x -B -G --nice --export " != *" ${arrayC[$i-1]} "* ]] && slurmScript="${array[$i]}" && continue; [ ! -z "$slurmScript" ] && slurmScriptParas="$slurmScriptParas ${array[$i]}" || additionalPara="$additionalPara ${array[$i]}"; };;
    esac

    # [ -z "$wrapCMD" ] && [[ ${array[$i]} == "--wrap" ]] && echoerr Found --wrap ${array[$i+1]} && wrapCMD="${array[$i+1]}" && array[$i]="" && array[$i+1]=""
  	# [ -z "$wrapCMD" ] && [[ ${array[$i]} ==	--wrap=* ]] && echoerr Found --wrap= && wrapCMD="${array[$i]}" && wrapCMD="${wrapCMD/--wrap=}" && array[$i]="" || CMDWithoutWrap="$CMDWithoutWrap ${array[$i]}"
    # [ -z "$slurmScript$wrapCMD" ] && [ -f "${array[$i]}" ] && echoerr Found slurmScript ${array[$i]} && slurmScript="${array[$i]}" && array[$i]=""
    # [ -z "$slurmScript" ] && CMDWithoutSlurmCMD="$CMDWithoutSlurmCMD ${array[$i]}" || slurmScriptParas="$slurmScriptParas ${array[$i]}"
done

echoerr Parsing result from sbatch commandline:
echoerr program: $program ref: $ref jobFlag: $flag inputs: $inputs deps: $deps
echoerr partition: $partition time: $time mem: $mem mem-per-cpu: $mem1 task: $task core: $core node: $node out: $out err: $err additionalPara: $additionalPara

if [ -z "$slurmScript" ]; then
    echoerr wrapCMD: $wrapCMD
else
    echoerr slurmScript: $slurmScript
    echoerr slurmScriptParas: $slurmScriptParas
fi
#echoerr test or run: $testRun
#echoerr

if [ ! -z "$slurmScript" ]; then
  echoerr Validating slurmScript:
  firstRow=`head -n 1 $slurmScript`
  echoerr FirstRow of the script: $firstRow
  [[ "$firstRow" =~ ^#\!/bin/bash ]] || [[ "$firstRow" =~ ^#\!/usr/bin/bash ]] || [[ "$firstRow" =~ ^#\!/bin/sh ]] || [[ "$firstRow" =~ ^#\!/usr/bin/sh ]]|| { echo "Error: first row of the slurm script ($slurmScript) does not start with #!/bin/bash, #!/usr/bin/bash, #!/bin/sh, #!/usr/bin/sh";  exit 1;}


  # get program (rule), input, threads from
  if [[ "$parentCmd" == */bin/snakemake* ]]; then
    #properties='{"local": false, "params": {}, "input": ["data/genome.fa", "data/samples/A.fastq"], "rule": "bwa_map", "resources": {}, "output": ["mapped_reads/A.bam"], "cluster": {}, "threads": 2}'
    properties=`grep "^# properties =" $slurmScript`
    echoerr Got snakemake parameter: $properties
    properties=${properties#\# properties =}
    
    [ -z "$program" ] && program=$(echo $properties | jq -r '.rule')
    [ -z "$core" ] && core=$(echo $properties | jq -r '.threads')
    if [ -z "$inputs" ]; then 
        inputs=$(echo $properties | jq -r '.input | join(" ")') 
    
        # need to ignore non exist file
        tm=""
        for fx in $inputs; do   
            [ -e "$fx" ] && tm="$tm $fx"
        done    
        inputs=$tm
    fi 

    
    flag=$program.${inputs##*/}
    
    echoerr Got snakemake parameter: rule: $program, threads: $core input: $inputs 

  # todo: need to check if there two inputs   
  elif  [ -f .command.sh ] && [ -f .command.run ] ; then 

    for i in `awk '/nxf_stage\(\) {/, /\}/{if ($0 ~ /ln -s /) print $0; }' .command.run | cut -d' ' -f7 `; do
        inputs="$inputs $i"
    done 

    #inputs=`grep -A 2 "# stage input files" .command.run | tail -n 1 | tr -s " " | cut -d' ' -f4`; 
    
    program=`grep "# NEXTFLOW TASK:" .command.run`; program=${program#*NEXTFLOW TASK: }; 

    flag=$program.${inputs##*/}

    echoerr Got parameter for nextflow: program: $program flag: $flag inputs: $inputs
  fi  
fi

[ -z "$wrapCMD" ] && [ -z "$slurmScript" ] && echo Error: Did not find --wrap, did not find slurmScript either. && exit 1

if [[ ! -z "$slurmScript" ]]; then
    echoerr
    echoerr Parsing slurm script ...
    while IFS=$'\n' read line; do
        #echoerr Got line: $line
        arr=($line)
        if [[ ${arr[0]} == "#SBATCH" ]]; then
            case "${arr[1]}" in
                --comment=*    )   echoerr found --comment && comment="${arr[1]}";; # && comment=${comment/--comment=/};; # && for c in "${ar[@]}"; do v=${c%=*}; [ -z "${!v}" ] && eval $c; done;;

                "-A"   )   echoerr Found -A  && [ -z "$slurmAcc" ] && slurmAcc="-A ${arr[2]}";;
                "--mem" ) [ -z "$mem" ] && echoerr Found --mem && mem="${arr[2]}";;
                --mem=* ) [ -z "$mem" ] && echoerr Found --mem= && mem="${arr[1]}" && mem="${mem/--mem=/}";;
                "--mem-per-cpu") [ -z "$mem1" ] && echoerr Found --mem-per-cpu && mem1="${arr[2]}";;
                --mem-per-cpu=*) [ -z "$mem1" ] && echoerr Found --mem-per-cpu= && mem1="${arr[1]}" && mem1="${mem1/--mem-per-cpu=/}";;
                "-c"    ) [ -z "$core" ] && echoerr Found -c && core="${arr[2]}";;
                --cpus-per-task=* ) [ -z "$core" ] && echoerr Found --cpus-per-task= && core="${arr[1]}" && core="${core/--cpus-per-task=/}";;
                "-n"    ) [ -z "$task" ] && echoerr Found -n && task="${arr[2]}";;
                --ntasks=* ) [ -z "$task" ] && echoerr Found --ntasks= && task="${arr[1]}" && task="${task/--ntasks=/}";;
                "-N"    ) [ -z "$node" ] && echoerr Found -N && node="${arr[2]}";;
                --nodes=*   ) [ -z "$node" ] && echoerr Found --nodes= && node="${arr[1]}" && node="${node/--nodes=/}";;
                "-t"    ) [ -z "$time" ] && echoerr Found -t && time="${arr[2]}";;
                --time=* ) [ -z "$time" ] && echoerr Found --time= && time="${arr[1]}" && time="${time/--time=/}";;
                "-J"    ) [ -z "$name" ] && echoerr Found -J && name="${arr[2]}";;
                --job-name=* ) [ -z "$name" ] && echoerr Found --job-name= && name="${arr[1]}" && name="${time/--job-name=/}";;
                "-o"    ) [ -z "$out" ] && echoerr Found -o && out="${arr[2]}";;
                --output=* ) [ -z "$out" ] && echoerr Found --output= && out="${arr[1]}" && out="${out/--output=/}";;
                "-e" ) [ -z "$err" ] && echoerr Found -e && err="${arr[2]}";;
                --error=* ) [ -z "$err" ] && echoerr Found --error= && err="${arr[1]}" && err="${err/--error=/}";;

                "-d"    ) [ -z "$deps" ] && echoerr Found -d && deps="${arr[2]}";;
                --dependency=* ) [ -z "$deps" ] && echoerr Found --dependency= && deps="${arr[1]}" && deps="${deps/--dependency=/}";;
            esac
        fi
    done < $slurmScript

    echoerr
    echoerr Parsing result from slurm script:
    echoerr time: $time mem: $mem mem-per-cpu: $mem1 task: $task core: $core node: $node out: $out err: $err deps: $deps

    echoerr slurmScript: $slurmScript
    echoerr slurmScriptParas: $slurmScriptParas
fi

[ -z "$core" ] && { [ ! -z "$task" ] && core="$task" || core=1; }

#[ -z "$program" ] && [ ! -z "$P" ] && program="$P"
#[ -z "$ref" ] && [ ! -z "$R" ] && ref="$R"
#[ -z "$input" ] && [ ! -z "$I" ] && inputs="$I"
#[ -z "$deps" ] && [ ! -z "$d" ] && deps="$d"
#[[ "$program" == *Checkpoint ]] || program=$program$C
[ -z "$out" ] || out=`realpath $out`
[ -z "$err" ] || err=`realpath $err`
#todo: err is ignored. will work on it later

if [ -z "$flag" ]; then   
    tm=`mktemp XXXXXXXX --dry-run`
    if [ ! -z "$program" ]; then 
        [ ! -z "$inputs" ] && flag=$program.${inputs##*/} || flag=$program.$tm  
    elif [ ! -z "$name" ]; then 
        flag=$name.$tm
    else 
        # todo: figure out how to find program name and input automatically
        if [ ! -z "$wrapCMD" ]; then
            arr=($wrapCMD)
            # need remove the parameters, to get unique name for job types
            if [[ "${arr[0]}" == sh || "${arr[0]}" == bash || "${arr[0]}" == python || ""${arr[0]}" == python3  || "${arr[0]}" == matlab || ""${arr[0]}" == Rscript ]]; then
                program="${arr[1]}"
            else
                program="${arr[0]}"
            fi
            program=${program##*/}
            #flag=${wrapCMD// /.}; flag=${flag%;}.$tm
        else
            program=${slurmScript##*/}
            #flag=$program.${slurmScriptParas// /.}
            #flag=${flag%.}.$tm
        fi
        [ ! -z "$inputs" ] && flag=$program.${inputs##*/} || flag=$program.$tm 

    fi
fi

flag=${flag//\//-}
flag=${flag//../.}

[ -z "$name" ] && name=$flag

mkdir -p $smartSlurmLogDir/
job=$smartSlurmLogDir/$flag.sh
succFlag=$smartSlurmLogDir/$flag.success
failFlag=$smartSlurmLogDir/$flag.failed
#startFlag=$smartSlurmLogDir/$flag.start
killFlag=$smartSlurmLogDir/$flag.user.killed
outFlag=$smartSlurmLogDir/$flag.out
errFlag=$smartSlurmLogDir/$flag.err

#deps=${deps#.}

[ -z "$ref" ] && ref=none
[ -z "$inputs" ] && inputs=none
if [ -z "$deps" ]; then 
    depsID=null
else 
    depsID=${deps#afterok:}; 
    deps="-d $deps"
fi
[ -z "$chdir" ] || chdir="-D $chdir"

if [ -z "$mem" ]; then
    if [ -z "$mem1" ]; then
        mem=$defaultMem
    else
        [[ "$mem1" == *G ]] && mem=$(( ${mem1%G} * $core * 1024)) || mem=$(( ${mem1%M} * $core ))
    fi
else
    [[ "$mem" == *G ]] && mem=$(( ${mem%G} * 1024 )) || mem=${mem%M}
fi

if [ -z "$time" ]; then
    min=$defaultTime
else
    [[ "$time" == *-* ]] && { day=${time%-*}; tem=${time#*-}; hour=${tem%%:*}; min=${tem#*:}; min=${min%%:*}; sec=${tem#$hour:$min}; sec=${sec#:}; } || { [[ "$time" =~ ^[0-9]+$ ]] && min=$time || { sec=${time##*:}; min=${time%:*}; min=${min##*:}; hour=${time%$min:$sec}; hour=${hour%:}; day=0;} }

    [ -z "$day" ] && day=0; [ -z "$hour" ] && hour=0; [ -z "$min" ] && min=0;[ -z "$sec" ] && sec=0

    min=$(($day * 24 * 60 + $hour * 60  + $min))
fi

resAjust="#Depend on: $depsID\n"

# if [[ $deps == "null" ]]; then
#     deps=""
#     resAjust="#Depend on no job\n"

# elif [[ $deps == ${deps/\./} ]]; then
#     resAjust="#Depend on single job\n"
#     deps="--dependency=afterok:${deps/\./}"
# else
#     resAjust="#Depend on multiple jobs\n"
#     tmp=""
#     for t in ${deps//\./ }; do
#         #echoerr working on $t  test
#         tmp="$tmp:$t"
#     done
#     [ ! -z "$tmp" ] && deps="--dependency=afterok$tmp"
# fi

resAjust="#Original mem $mem M, Original time: $min mins\n"

# check job done before or not
if [[ "$testRun" == "run" ]]; then
    
    #echoerr ps command output 
    
    if [ -f $succFlag ] && [[ "$parentCmd" != */bin/snakemake* ]] && [ ! -f .command.run ] && [ ! -f command.sh ]; then
        stepID=${flag%%.*}
        if ([ -f $smartSlurmLogDir/skipAllSuccessJobs.txt ] || [ -f $smartSlurmLogDir/skipAllSuccessJobs$stepID.txt ]) && [ -z "$deps" ]; then
            echo skipThisJob
            exit
        elif [ ! -f $smartSlurmLogDir/reRunAllSuccessJobs.txt ] && [ ! -f $smartSlurmLogDir/reRunAllSuccessJobs$stepID.txt ] && [ -z "$deps" ]; then
            stepName=${flag#*.}; stepName=${stepName#*.}; stepName=${stepName%%.*}
            echoerr1 $flag was done before, do you want to re-run it?
            echoerr1 -e "y:        To re-run this job, press y, then enter key."
            [ -z "$runningSingleJob" ] && echoerr1 -e "ystep:    To re-run all jobs for step $stepID: $stepName, type ystep, then press enter key."
            [ -z "$runningSingleJob" ] && echoerr1 -e "yall:     To re-run all jobs, type yall, then press enter key."
            echoerr -e "enter:    To not re-run this job, directly press enter key."
            [ -z "$runningSingleJob" ] && echoerr1 -e "nstep:    To not re-run all successful jobs for step $stepID: $stepName, type nstep, then press enter key."
            [ -z "$runningSingleJob" ] && echoerr1 -e "nall:     To not re-run all successful jobs, type nall, then press enter key."
            read -p "" x </dev/tty

            echoerr1 You typed: \"$x\"

            if [[ "$x" == "y" ]]; then
                [ -z "$runningSingleJob" ] && echoerr "Will re-run the down stream steps even if they are done before (because they have deps - see code in row 70)."
                [ -z "$runningSingleJob" ] && echoerr1 "Will re-run the down stream steps even if they are done before (because they have deps - see code in row 70)."
            elif [[ "$x" == "ystep" ]]; then
                touch $smartSlurmLogDir/reRunAllSuccessJobs$stepID.txt
            elif [[ "$x" == "nstep" ]]; then
                touch $smartSlurmLogDir/skipAllSuccessJobs$stepID.txt
                echo skipThisJob
                exit
            elif [[ "$x" == "yall" ]]; then
                touch $smartSlurmLogDir/reRunAllSuccessJobs.txt
            elif [[ "$x" == "nall" ]]; then
                touch $smartSlurmLogDir/skipAllSuccessJobs.txt
                echo skipThisJob
                exit
            else
                echo skipThisJob
                exit
            fi
        fi
        #echoerr
        #rm $succFlag
        #ls -l $smartSlurmLogDir/ 1>&2
    fi

    if [ -f $smartSlurmLogDir/keepRunningJobs.txt ]; then
        #echoerr $flag
        # check the third column for the job name, then find the the job id in column 1
        out=`cat $smartSlurmLogDir/allJobs.txt | awk '{if ($3 ~ /'"$flag/"') print $1, $2;}' | tail -n 1`
        id=${out%% *}; de=${out##* }

        # this job was submitted earlier
        if [ ! -z "$id" ]; then

            # this job still running or pending
            # todo: need troubleshoot this part. If dependency condition changed, need resubmit the job
            if grep ^$id $smartSlurmLogDir/keepRunningJobs.txt 1>/dev/null; then

                # check if upsteam job is resubmitted
                upstreamResubmitted=""
                if [[ $de != null ]]; then
                    if [[ $de == ${de/:/} ]]; then
                        [ -f $smartSlurmLogDir/${de/:/}.resubmitted ] && upstreamResubmitted=y
                    else
                        for t in ${de//:/ }; do
                            [ -f $smartSlurmLogDir/$t.resubmitted ] && upstreamResubmitted=y
                        done
                    fi
                fi

                # the job still running, but the dependency contidtion is change, need cancel it.
                if [ ! -z "$upstreamResubmitted" ] || [[ "$de" != "$depsID" ]]; then
                    echoerr "Running job job is cancelled: $id" && scancel $id
                else

                    echoerr thisJobStillRunning $id here
                    echo $id
                    exit
                fi
            else
                if [ -z "$deps" ]; then
                    depsR="" #dependency=\"\""
                    #resAjust="#Depend on

                #elif [[ $depsID == ${depsID/:/} ]]; then
                #    #resAjust="#Depend on single job\n"
                #    depsR="dependency=afterok:${depsID/:/}"
                else
                    #resAjust="#Depend on multiple jobs\n"
                    # tmp=""
                    # for t in ${depsID//:/ }; do
                    #     #echoerr working on $t  test
                    #     tmp="$tmp:$t"
                    # done
                    # [ ! -z "$tmp" ] && depsR="dependency=afterok$tmp"

                    epsR="dependency=$deps"
                fi

                scontrol requeue $id

                touch $smartSlurmLogDir/$id.resubmitted

                [ ! -z "$depsR" ] && scontrol update job=$id $depsR && scontrol hold $id 
                echo "Resubmit#$id-#$de-#$flag" >> $smartSlurmLogDir/allJobs.txt
                echo $id
                #rm $smartSlurmLogDir/$flag.success 2>/dev/null
                rm -r $smartSlurmLogDir/$flag $smartSlurmLogDir/$flag.* 2>/dev/null || : # in case user has set -e and file not exist,  we will not exit
                exit
            fi
        fi
        #set +x
    fi
    #ls $flagDir/$flag.* 1>&2
    rm -r $smartSlurmLogDir/$flag $smartSlurmLogDir/$flag.* 2>/dev/null || : # in case user has set -e and file not exist,  we will not exit
fi

#[ -f $smartSlurmJobRecordDir/stats/extraMem.$program.$ref ] && maxExtra=`sort -n $smartSlurmJobRecordDir/stats/extraMem.$program.$ref | tail -n1 | cut -d' ' -f1` && oomCount=`wc -l $smartSlurmJobRecordDir/stats/extraMem.$program.$ref | cut -d' ' -f1` && extraMem=$(( $maxExtra * $oomCount ))

[ -f $smartSlurmJobRecordDir/stats/extraMem.$program.$ref ] && maxExtra=`sort -n $smartSlurmJobRecordDir/stats/extraMem.$program.$ref | tail -n1 | cut -d' ' -f1` && extraMem=$(( $maxExtra * 2 )) || extraMem=$(($defaultExtraMem * 2))

 #extraMem=`sort $smartSlurmJobRecordDir/stats/extraMem.$program.$ref | tail -n1`

# original mem and time
memO=$mem; minO=$min
inputSize=0

echoerr
echoerr Check if there input file list and this job does not depend on other jobs

# has checkpoint ready to resume from
if [ -z "$deps" ] && [[ "$program" == *.Checkpoint ]] && ls $smartSlurmLogDir/$flag/ckpt_*.dmtcp >/dev/null 2>&1 && [ -f $smartSlurmLogDir/$flag.adjust ]; then
    tText=`cat $smartSlurmLogDir/$flag.adjust`
    mem=`echo $tText | cut -d' ' -f1`
    min=`echo $tText | cut -d' ' -f2`
    extraMem=`echo $tText | cut -d' ' -f3`
    echo Got mem/time from $smartSlurmLogDir/$flag.adjust: > $outFlag
    echo $tText >> $outFlag
else
    rm -r $smartSlurmLogDir/$flag $smartSlurmLogDir/$flag.* 2>/dev/null || :

    # do not have input file list    #todo: we can try to figure out the input from the program parameters here
    if [[ "$inputs" == "none" ]]; then
        echoerr No inputs
        resAjust="$resAjust#This job does not have input.\n"
        ref=${ref//\//-}

        rows=`( wc -l $smartSlurmJobRecordDir/stats/$program.$ref.mem.noInput 2>/dev/null || echo 0 ) | awk '{print $1}'`
        #echoerr rows  $rows
        # empty or more than 60 minutes but less than 4 records
        if test `find $smartSlurmJobRecordDir/stats/$program.$ref.mem.noInput -mmin +20 2>/dev/null` && [ $rows -lt 200 ] || [ $rows -eq 0 ]; then
            mkdir -p $smartSlurmJobRecordDir/stats/
            #cat /home/*/smartSlurm/stats/myJobRecord.txt > $smartSlurmJobRecordDir/jobRecord.txt
            #test `find $smartSlurmJobRecordDir/jobRecord.txt -mmin +20` && echoerr jobRecord.txt synced within 21 hour. No need to re-sync || cat $HOME/smartSlurm/myJobRecord.txt > $smartSlurmJobRecordDir/jobRecord.txt

            # todo: could use single file here
            grep COMPLETED $smartSlurmJobRecordDir/jobRecord.txt 2>/dev/null | awk -F, -v a=$program -v b=$ref '{ if($12 == a && $13 == b) {print $7 }}' | uniq > $smartSlurmJobRecordDir/stats/$program.$ref.mem.noInput

            grep COMPLETED $smartSlurmJobRecordDir/jobRecord.txt 2>/dev/null | awk -F, -v a=$program -v b=$ref '{ if($12 == a && $13 == b) {print $8 }}' | uniq > $smartSlurmJobRecordDir/stats/$program.$ref.time.noInput

            if [ -s $smartSlurmJobRecordDir/stats/$program.$ref.mem.noInput ] && [ -s $smartSlurmJobRecordDir/stats/$program.$ref.time.noInput ]; then

                #OUT="$(mktemp -d)"
                paste $smartSlurmJobRecordDir/stats/$program.$ref.time.noInput $smartSlurmJobRecordDir/stats/$program.$ref.mem.noInput | column -s $'\t' -t | sed '$ d' > $smartSlurmJobRecordDir/stats/$program.$ref.timeMem.noInput

                #cp $smartSlurmJobRecordDir/stats/$program.$ref.timeMem.noInput $OUT/timeMem.txt
                #cd $OUT
                # make plot and calculate statistics
                gnuplot -e 'set key outside; set key reverse; set key invert; set term png; set output "'"$smartSlurmJobRecordDir/stats/$program.$ref.stat.noInput.png"'"; set title "Time vs. Memory Usage"; set xlabel "Time(Min)"; set ylabel "Memory(M)"; f(x)=a*x+b; fit f(x) "'"$smartSlurmJobRecordDir/stats/$program.$ref.timeMem.noInput"'" u 1:2 via a, b; t(a,b)=sprintf("f(x) = %.2fx + %.2f", a, b); plot "'"$smartSlurmJobRecordDir/stats/$program.$ref.timeMem.noInput"'" u 1:2,f(x) t t(a,b); print "Finala=", a; print "Finalb=",b; stats "'"$smartSlurmJobRecordDir/stats/$program.$ref.timeMem.noInput"'" u 1 ' 2>&1 | grep 'Final\| M' | awk 'NF<5{print $1, $2}' | sed 's/:/=/' | sed 's/ //g' > $smartSlurmJobRecordDir/stats/$program.$ref.timeMem.stat.noInput

                # cd -
                # convert $OUT/timeMem.pdf -background White -flatten $smartSlurmJobRecordDir/stats/$program.$ref.stat.noInput.pdf 2>/dev/null
                # pdftoppm $smartSlurmJobRecordDir/stats/$program.$ref.stat.noInput.pdf  -png > $smartSlurmJobRecordDir/stats/$program.$ref.stat.noInput.png 2>/dev/null

                rows=`{ wc -l $smartSlurmJobRecordDir/stats/$program.$ref.mem.noInput 2>/dev/null || echo 0; } | cut -f 1 -d " "`

                #rm -r $OUT 2>/dev/null
            fi
        fi

        # at least 3 records
        if [ $rows -ge 3 ]; then
            cutoffRow=$(( ($row - 1)  / 10  + 1)) # top 10

            mem=`cat $smartSlurmJobRecordDir/stats/$program.$ref.mem.noInput | sort -nr | tr '\n' ' ' | cut -f $cutoffRow -d " "`
            mem=$((${mem/\.*/} + extraMem))
            # mem=${txt[$cutoffRow]}; mem=$((${mem/\.*/} + extraMem))M

            min=`cat $smartSlurmJobRecordDir/stats/$program.$ref.time.noInput | sort -nr | tr '\n' ' ' | cut -f $cutoffRow -d " "`
            min=$((${min/\.*/} + defaultExtraTime))

            resAjust="$resAjust#Got estimation based on program.reference: $program.$ref.\n"
            resAjust="$resAjust#Give ${extraMem} M extra memory and $defaultExtraTime more minutes. \n#So use this to submit the job: $mem M  ${min} min"

        else
            resAjust="$resAjust#There are less than 3 job records. Use default mem and time."
        fi


    # not depends on other job, and there is input file list
    elif [ -z "$deps" ]; then
        echoerr Does not depend on other jobs and have inputs.
        inputSize=`{ du --apparent-size -c -L ${inputs//,/ } 2>/dev/null || echo notExist; } | tail -n 1 | cut -f 1`

        if [[ "$inputSize" == "notExist" ]]; then
            resAjust="$resAjust#Some or all input files not exist: $inputs\n"
            echoerr Error! missingInputFile: ${inputs//,/ }
            
            [[ "$testRun" == "run" ]] && exit
        else
            #inputSize=$(($inputSize/1024)); # convert to M
            resAjust="$resAjust#InputSize: $inputSize\n"

            #rm ~/.rcbio/$program.$ref.mem.stat # for testing

            if [ -f $smartSlurmJobRecordDir/stats/$program.${ref//\//-}.mem.stat ]; then
                if [ -f .command.sh ] && [ -f .command.run ]; then 
                    output=`estimateMemTime.sh $program ${ref//\//-} $inputSize 2>> ../../../.nextflow.log`
                else #if [[ "$parentCmd" == */bin/snakemake* ]]; then
                    output=`estimateMemTime.sh $program ${ref//\//-} $inputSize 2>> .smartSlurm.log`
                fi 

                
                resAjust="$resAjust\n`cat $smartSlurmJobRecordDir/stats/$program.${ref//\//-}.mem.stat`\n"
                resAjust="$resAjust\n#Output from estimateMemTime.sh: $output \n"
                
                if [[ "$output" == "outOfRange" ]]; then
                    resAjust="$resAjust#Input size is too big for the curve to estimate! Use default mem and runtime to submit job.\n"
                    # not deleting mem.stat, so other jobs will not re-build it within 60 minutes
                elif [ ! -z "$output" ]; then
                    output=${output% *}
                    [[ ${output% *} != 0 ]] && mem=$((${output% *}+extraMem)) && resAjust="$resAjust\n#Give ${extraMem} M extra memory. "
                    [[ ${output#* } != 0 ]] && min=$((${output#* }+defaultExtraTime)) && resAjust="$resAjust\n#Give $defaultExtraTime mins more time."
                    resAjust="$resAjust\n#So use this to submit the job: $mem M ${min} mins"
                fi

            fi

            if [[ "$output" == "outOfRange" ]] && test `find $smartSlurmJobRecordDir/stats/$program.${ref//\//-}.mem.stat -mmin +20` || [ ! -f $smartSlurmJobRecordDir/stats/$program.${ref//\//-}.mem.stat ]; then
                resAjust="$resAjust#Do not have a formula, or it is old and out of range. Let us build one...\n"

                #[ test `find $smartSlurmJobRecordDir/jobRecord.txt -mmin -20` ] && echoerr jobRecord.txt synced within 20 hour. No need to re-sync || cat $HOME/smartSlurm/myJobRecord.txt > $smartSlurmJobRecordDir/jobRecord.txt

                #jobStatistics.sh $program ${ref//\//-} 4 1>&2


                #filter by program and reference
                # todo: maybe able to replace / in ref at begaining of the script?
                ref=${ref//\//-}
                grep COMPLETED $smartSlurmJobRecordDir/jobRecord.txt 2>/dev/null | awk -F"," -v a=$program -v b=$ref '{ if($12 == a && $13 == b && $2!=0) {print $2, $7 }}' | sort -r  -k1,1 -k2,2 | sort -u -k1,1 > $smartSlurmJobRecordDir/stats/$program.$ref.mem

                grep COMPLETED $smartSlurmJobRecordDir/jobRecord.txt 2>/dev/null | awk -F"," -v a=$program -v b=$ref '{ if($12 == a && $13 == b && $2!=0) {print $2, $8 }}' | sort -r  -k1,1 -k2,2 | sort -u -k1,1 > $smartSlurmJobRecordDir/stats/$program.$ref.time

                echoerr "Got mem data from jobRecord.txt (content of $smartSlurmJobRecordDir/stats/$program.$ref.mem):"
                echoerr `cat $smartSlurmJobRecordDir/stats/$program.$ref.mem`

                echoerr "Got time data from jobRecord.txt (content of $smartSlurmJobRecordDir/stats/$program.$ref.time):"
                echoerr `cat $smartSlurmJobRecordDir/stats/$program.$ref.time`

                if [[ $(wc -l <$smartSlurmJobRecordDir/stats/$program.$ref.mem) -lt 3 ]]; then
                    echoerr There are less than 3 records. No way to fit a curve. User defaut values...
                    resAjust="$resAjust#There are less than 3 records. No way to fit a curve.\n"

                else



                    # make plot and calculate statistics
                    # gnuplot -e 'set term pdf; set output "mem.pdf"; set title "Input Size vs. Memory Usage" font "Helvetica Bold,18"; set xlabel "Input Size(K)"; set ylabel "Memory Usage(M)"; f(x)=a*x+b; fit f(x) "mem.txt" u 1:2 via a, b; t(a,b)=sprintf("f(x) = %.2fx + %.2f", a, b); plot "mem.txt" u 1:2,f(x) t t(a,b); print "Finala=", a; print "Finalb=",b; stats "mem.txt" u 1 ' 2>&1 | grep 'Final\| M' | awk 'NF<4{print $1, $2}' |sed 's/:/=/' | sed 's/ //g' > mem.stat.txt; echo STDFIT=`cat fit.log | grep FIT_STDFIT | tail -n 1 | awk '{print $8}'` >> mem.stat.txt
                    # echo RSquare="$(gnuplot -e 'stats "mem.txt" using 1:2;' 2>&1| grep Correlation | cut -d' ' -f7 | awk '{print $1 * $1 }')" >> mem.stat.txt

                    gnuplot -e 'set key outside; set key reverse; set key invert; set term png; set output "'"$smartSlurmJobRecordDir/stats/$program.$ref.mem.png"'"; set title "Input Size vs. Memory Usage"; set xlabel "Input Size(K)"; set ylabel "Memory Usage(M)"; f(x)=a*x+b; fit f(x) "'"$smartSlurmJobRecordDir/stats/$program.$ref.mem"'" u 1:2 via a, b; t(a,b)=sprintf("f(x) = %.2fx + %.2f", a, b); plot "'"$smartSlurmJobRecordDir/stats/$program.$ref.mem"'" u 1:2,f(x) t t(a,b); print "Finala=", a; print "Finalb=",b; stats "'"$smartSlurmJobRecordDir/stats/$program.$ref.mem"'" u 1 ' 2>&1 | grep 'Final\| M' | awk 'NF<5{print $1, $2}' |sed 's/:/=/' | sed 's/ //g' > $smartSlurmJobRecordDir/stats/$program.$ref.mem.stat ; echo STDFIT=`cat fit.log | grep FIT_STDFIT | tail -n 1 | awk '{print $8}'` >> $smartSlurmJobRecordDir/stats/$program.$ref.mem.stat

                    echo RSquare="$(gnuplot -e 'stats "'"$smartSlurmJobRecordDir/stats/$program.$ref.mem"'" using 1:2;' 2>&1| grep Correlation | cut -d' ' -f7 | awk '{print $1 * $1 }')" >> $smartSlurmJobRecordDir/stats/$program.$ref.mem.stat
                  
                    echo SCount=$(wc -l $smartSlurmJobRecordDir/stats/$program.$ref.mem | cut -d' ' -f1) >> $smartSlurmJobRecordDir/stats/$program.$ref.mem.stat

                    sed -i 's/\x0//g' $smartSlurmJobRecordDir/stats/$program.$ref.mem.stat

                    gnuplot -e 'set key outside; set key reverse; set key invert; set term png; set output "'"$smartSlurmJobRecordDir/stats/$program.$ref.time.png"'"; set title "Input Size vs. Time Usage"; set xlabel "Input Size(K)"; set ylabel "Time Usage(Min)"; f(x)=a*x+b; fit f(x) "'"$smartSlurmJobRecordDir/stats/$program.$ref.time"'" u 1:2 via a, b; t(a,b)=sprintf("f(x) = %.2fx + %.2f", a, b); plot "'"$smartSlurmJobRecordDir/stats/$program.$ref.time"'" u 1:2,f(x) t t(a,b); print "Finala=", a; print "Finalb=",b; stats "'"$smartSlurmJobRecordDir/stats/$program.$ref.time"'" u 1 ' 2>&1 | grep 'Final\| M' | awk 'NF<5{print $1, $2}' |sed 's/:/=/' | sed 's/ //g' > $smartSlurmJobRecordDir/stats/$program.$ref.time.stat ; echo STDFIT=`cat fit.log | grep FIT_STDFIT | tail -n 1 | awk '{print $8}'` >> $smartSlurmJobRecordDir/stats/$program.$ref.time.stat

                    echo RSquare="$(gnuplot -e 'stats "'"$smartSlurmJobRecordDir/stats/$program.$ref.time"'" using 1:2;' 2>&1| grep Correlation | cut -d' ' -f7 | awk '{print $1 * $1 }')" >> $smartSlurmJobRecordDir/stats/$program.$ref.time.stat

                    sed -i 's/\x0//g' $smartSlurmJobRecordDir/stats/$program.$ref.time.stat



                    # make plot and calculate statistics
                    # gnuplot -e 'set term pdf; set output "time.pdf"; set title "Input Size vs. Time Usage" font "Helvetica Bold,18"; set xlabel "Input Size(K)"; set ylabel "Time(Min)"; f(x)=a*x+b; fit f(x) "time.txt" u 1:2 via a, b; t(a,b)=sprintf("f(x) = %.2fx + %.2f", a, b); plot "time.txt" u 1:2,f(x) t t(a,b); print "Finala=", a; print "Finalb=",b; stats "time.txt" u 1 ' 2>&1 | grep 'Final\| M' | awk 'NF<4{print $1, $2}' |sed 's/:/=/' | sed 's/ //g' > time.stat.txt; echo STDFIT=`cat fit.log | grep FIT_STDFIT | tail -n 1 | awk '{print $8}'` >> time.stat.txt
                    # echo RSquare="$(gnuplot -e 'stats "time.txt" using -1:2;' 2>&1| grep Correlation | cut -d' ' -f7 | awk '{print $1 * $1 }')" >> time.stat.txt

                    echoerr There are more than 3 $program $ref jobs already run for this program, statics is ready for current job:
                    # echoerr Memeory statisics:
                    # echoerr "inputsize mem(M)"
                    # cat $smartSlurmJobRecordDir/stats/$program.$ref.mem.stat
                    # echoerr
                    # echoerr Time statistics:
                    # echoerr "inputsize time(minute)"
                    # cat $smartSlurmJobRecordDir/stats/$program.$ref.time.stat

                    # mv $OUT/mem.txt $smartSlurmJobRecordDir/stats/$program.$ref.mem
                    # mv $OUT/time.txt $smartSlurmJobRecordDir/stats/$program.$ref.time

                    # convert $OUT/mem.pdf -background White -flatten $smartSlurmJobRecordDir/stats/$program.$ref.mem.pdf
                    # convert $OUT/time.pdf -background White -flatten $smartSlurmJobRecordDir/stats/$program.$ref.time.pdf
                    # pdftoppm $smartSlurmJobRecordDir/stats/$program.$ref.mem.pdf  -png > $smartSlurmJobRecordDir/stats/$program.$ref.mem.png
                    # pdftoppm $smartSlurmJobRecordDir/stats/$program.$ref.time.pdf  -png > $smartSlurmJobRecordDir/stats/$program.$ref.time.png

                    echoerr
                    echoerr You can see the plot using commands:
                    echoerr display $smartSlurmJobRecordDir/stats/$program.$ref.mem.png
                    echoerr display $smartSlurmJobRecordDir/stats/$program.$ref.time.png

                    # cd -



                    #echoerr got files in $smartSlurmJobRecordDir/stats:
                    #ls -lrt $smartSlurmJobRecordDir/stats
                    if [ -f $smartSlurmJobRecordDir/stats/$program.${ref//\//-}.mem.stat ]; then
                        if [ -f .command.sh ] && [ -f .command.run ]; then 
                            output=`estimateMemTime.sh $program ${ref//\//-} $inputSize 2>> ../../../.nextflow.log`
                        else #if [[ "$parentCmd" == */bin/snakemake* ]]; then
                            output=`estimateMemTime.sh $program ${ref//\//-} $inputSize 2>> .smartSlurm.log`
                        fi 

                        resAjust="$resAjust\n`cat $smartSlurmJobRecordDir/stats/$program.${ref//\//-}.mem.stat`\n"
                        resAjust="$resAjust\n#Output from estimateMemTime.sh: $output \n"
                        if [[ "$output" == "outOfRange" ]]; then
                            resAjust="$resAjust#Input size is too big for the curve to estimate! Use default mem and runtime to submit job.\n"
                            # not deleting mem.stat, so other jobs will not re-build it within 60 minutes
                        elif [ ! -z "$output" ]; then
                            output=${output% *}
                            [[ ${output% *} != 0 ]] && mem=$((${output% *}+extraMem)) && resAjust="$resAjust\n#Give ${extraMem} M extra memory. "
                            [[ ${output#* } != 0 ]] && min=$((${output#* }+defaultExtraTime)) && resAjust="$resAjust\n#Give $defaultExtraTime more minutes."
                            resAjust="$resAjust\n#So use this to submit the job: $mem M ${min} mins"

                        fi
                        #echoerr got estimation $output
                    fi

                fi
                #rm -r $OUT 2>/dev/null
            fi

        fi
    else
        echoerr Has input, but depends on other jobs
        resAjust="$resAjust#Use default mem and time. Has input, but this job depends on $deps."

    fi
    
fi

#echoerr

# for OOM memtesting
#mem=4000

[ -z "$min" ] && { echo did not find time limit; exit 1; }

[ -z "$mem" ] && { echo did not find mem limit; exit 1; }


#time=$min

#[[ "$time" == *-* ]] && { day=${time%-*}; tem=${time#*-}; hour=${tem%%:*}; min=${tem#*:}; min=${min%%:*}; sec=${tem#$hour:$min}; sec=${sec#:}; } || { [[ "$time" =~ ^[0-9]+$ ]] && min=$time || { sec=${time##*:}; min=${time%:*}; min=${min##*:}; hour=${time%$min:$sec}; hour=${hour%:}; day=0;} }

#[ -z "$day" ] && day=0; [ -z "$hour" ] && hour=0; [ -z "$min" ] && min=0;[ -z "$sec" ] && sec=0

# resAjust="$resAjust#day $day,  hour $hour,  min $min,  sec $sec"

#echoerr looking partition for hour: $hours ##

[ "$mem" -lt 100 ] && mem=100 && resAjust="$resAjust\n#Mem is reset to 100M. "

[ "$min" -lt 10 ] && min=10 && resAjust="$resAjust\n#Time is reset to 10min. "

echo -e "$resAjust\n" > $outFlag

# testtime
#min=10   # min
# testmem
# mem=20   # M

adjustPartition $(( ( $min + 59 ) / 60 )) $partition

seconds=$(( $min * 60 ))

#echoerr srun seconds: $seconds

time=`eval "echo $(date -ud "@$seconds" +'$((%s/3600/24))-%H:%M:%S')"`

#echoerr New time for srun: $timeN

#echoerr

# mem for srun is 10M less than sbatch command
#[[ "$mem" == *G ]] && memN=$(( 1024 * ${mem%G} -10 ))M || memN=$(( ${mem%M} -10 ))M

#[ ${memN%M} -le 1 ] && { echoerr Error: --mem for sbatch command should bigger than 11M; usage; }
#[ ${memN%M} -le 1 ] && mem=20M && memN=10M

#echoerr Mew mem for srun: $memN

#[ -z "$out" ] && out="slurm-\$SLURM_JOBID.out"
#[ -z "$err" ] && err="slurm-\$SLURM_JOBID.err"
if [ -f .command.sh ] && [ -f .command.run ]; then 
    cp .command.run .command.run.back 
    sed -i "s/exit \$exit_status/cleanUp.sh $flag $program \\\"${ref//\//-}\\\" $inputSize $core $memO $minO $mem $min $partition  \\\"${slurmAcc#*-A }\\\" \\\"${inputs##*/}\\\" $extraMem $defaultExtraTime\n    exit \$exit_status/" .command.run >> ../../../.nextflow.log
    
    sed -i "s/\$NXF_ENTRY/memCpuMonitor.sh $flag $mem $min $memO $core \&\n\$NXF_ENTRY/" .command.run >> ../../../.nextflow.log

    diff .command.run .command.run.back >> ../../../.nextflow.log
    cmd="/usr/bin/sbatch --mail-type=FAIL --requeue --parsable -p $partition --mem $mem -t $time .command.run"
    
else 

    echoerr
    echoerr Building new sbatch command ...
    #echo "#!/bin/bash" > $job

    #echo -e "#!/bin/bash\ndate\n\ntrap \"{ cleanUp.sh $flag $program ${ref//\//-} $inputSize $core $memO $minO $mem $min $partition  \\\"${slurmAcc#*-A }\\\" \\\"$inputs\\\" $extraMem $defaultExtraTime; }\" EXIT\nmemCpuMonitor.sh $flag $mem $min $memO $core &\n" > $job

    # echo -e "#!/bin/bash\ndate\n\ntrap \"{ cleanUp.sh $flag $program ${ref//\//-} $inputSize $core $memO $minO $mem $min $partition  \\\"${slurmAcc#*-A }\\\" \\\"$inputs\\\" $extraMem $defaultExtraTime; }\" EXIT\nmemCpuMonitor.sh $flag $program ${ref//\//-} $inputSize $core $memO $minO $mem $min $partition \\\"${slurmAcc#*-A }\\\" \\\"$inputs\\\" $extraMem $defaultExtraTime &\n" > $job


    # [[ "$mem" == *G ]] && totalM=$(( 1024 * ${mem%G})) || totalM=${mem%M}

    # #echo "set -x" >> $job
    # echo "if [ ! -f $smartSlurmLogDir/$flag.adjust ]; then" >> $job
    # echo "   totalM=$totalM" >> $job
    # echo "else " >> $job
    # echo "   totalM=\`cat $smartSlurmLogDir/$flag.adjust\`" >> $job
    # ## this might work if cluster is bussy, may need sleep more or get the data from a text file
    # #echo "   sleep 10; sacct=\`sacct --format=JobID,Submit,Start,End,MaxRSS,State,NodeList%30,Partition,ReqTRES%30,TotalCPU,Elapsed%14,Timelimit%14 --units=M -j \$SLURM_JOBID\` " #>> $job
    # #echo "   totalM=\${sacct#*,mem=}; totalM=\${totalM%%M,n*}" >> $job
    # echo "fi" >> $job

    # when submit jobs from interactive job, need this. Otherwise get error:
    # echo -e "unset SLURM_MEM_PER_CPU SLURM_MEM_PER_GPU SLURM_MEM_PER_NODE SLURM_CPU_BIND" >> $job


    #echo -e "unset SLURM_MEM_PER_CPU SLURM_MEM_PER_GPU SLURM_MEM_PER_NODE SLURM_CPU_BIND" >> $job


    if [ -z "$slurmScript" ]; then
        echo -e "#!/bin/bash\ndate\n\ntrap \"{ cleanUp.sh $flag $program ${ref//\//-} $inputSize $core $memO $minO $mem $min $partition  \\\"${slurmAcc#*-A }\\\" \\\"$inputs\\\" $extraMem $defaultExtraTime; }\" EXIT\nmemCpuMonitor.sh $flag $mem $min $memO $core &\n" > $job
        echo -e "unset SLURM_CPU_BIND" >> $job

        #echo "touch $startFlag" >> $job
        [[ "$wrapCMD" == *" " ]] && wrapCMD=${wrapCMD% *}; wrapCMD=${wrapCMD%;} # remove ending space and ; from command
        # `echo -e "$wrapCMD" | xargs echo -ne`; wrapCMD=${wrapCMD%;} # remove ending space and ; from command
        if [[ "$program" == *.Checkpoint ]]; then
            #echo -e "set -x\ntrap \"[ -f $succFlag ] || touch $failFlag\" EXIT SIGSTOP SIGHUP\nsh -e -c \"$wrapCMD\" && touch $succFlag" > ${job%.sh}.cmd

            # echo -e "set -x\nps -u $USER\nsh -eo -c \"$wrapCMD;\" &" > ${job%.sh}.cmd
            # echo -e "pid=\$$" >> ${job%.sh}.cmd
            # echo -e "echo \$pid > ${job%.sh}/taskProcessID.txt" >> ${job%.sh}.cmd
            # echo -e "ps -u $USER\nwait \$pid" >> ${job%.sh}.cmd
            # echo -e "echo \$? > ${job%.sh}/taskProcessStatus.txt" >> ${job%.sh}.cmd

            # note: if job OOT or OOM, .success or .failed will not be created
            echo -e "for i in \`seq 10000\`; do echo Program running time: \${i}0 s; sleep 10; done &\nsh -e -o pipefail  -c \"$wrapCMD; \" && touch $succFlag" > ${job%.sh}.cmd
            ##echo -e "x=\"\$wrapCMD\" >> ${job%.sh}.cmd
            ##echo -e "ps -fu $USER | grep -v awk | grep -v srun | grep -v grep | grep \"$flag\"" >> ${job%.sh}.cmd
            #echo -e "ps -fu $USER | grep \"$smartSlurmLogDir/$flag.cmd\"" >> ${job%.sh}.cmd

            #echo -e "pss=\`ps -fu $USER | grep \"$smartSlurmLogDir/$flag.cmd\"\`" >> ${job%.sh}.cmd
            ##echo -e "pid=\`ps -fu $USER | grep -v awk | grep -v srun | grep -v grep | grep \" $flag \" | tail -n 2 | head -n 1 | tr -s ' ' | cut -d' ' -f2\`" >> ${job%.sh}.cmd

            #echo -e "pid=\`echo -e \"\$pss\" | head -n 4 | tail -n 1 | tr -s ' ' | cut -d' ' -f2\`" >> ${job%.sh}.cmd
            #echo -e "echo \$pid > ${job%.sh}/taskProcessID.txt" >> ${job%.sh}.cmd
            #echo -e "ps -fu $USER\nwait" >> ${job%.sh}.cmd
            #echo -e "echo done > ${job%.sh}.success" >> ${job%.sh}.cmd

            # x="bash -e -c \"set -e; sh sleeping1.sh; \""
            # pid=`ps -f -u $USER | grep -v awk | grep -v srun | awk -v pat="$x"  '$0 ~ pat { print $2 }' | tail -n 1`
            # echo $pid > /n/scratch3/users/l/ld32/sethTest/log/2.0.useSomeMemTimeAccordingInputSize.sh.1/taskProcessID.txt

            #pid=$!
            # #echo $pid >
            # ps -f -u ld32
            # wait

            # echo done > /n/scratch3/users/l/ld32/sethTest/log/2.0.useSomeMemTimeAccordingInputSize.sh.1/taskProcessStatus.txt



            echo "srun -n 1 $slurmAcc sh -e -c \"checkpoint.sh \\\"sh ${job%.sh}.cmd\\\" $flag $mem $min $extraMem \"" >> $job
        else
            
            if [ ! -z "$runningSingleJob" ]; then

                echo "srun -n 1 $slurmAcc sh -e -o pipefail -c '$wrapCMD; ' && touch $succFlag " >> $job
            else 

                # !!! runAsPipeline jobs will land here !!!
                echo "srun -n 1 $slurmAcc sh -e -o pipefail -c '$wrapCMD; ' && touch $succFlag && adjustDownStreamJobs.sh" >> $job
                # echo "srun -n 1 $slurmAcc sh -e -o pipefail -c '$wrapCMD; ' && adjustDownStreamJobs.sh $smartSlurmLogDir && touch $succFlag " >> $job
            fi 
        fi
    else
        echo -e "#!/bin/bash" > $job

        grep "^#SBATCH" $slurmScript >> $job || true >> $job
        echo -e "\ndate\n\ntrap \"{ cleanUp.sh $flag $program ${ref//\//-} $inputSize $core $memO $minO $mem $min $partition  \\\"${slurmAcc#*-A }\\\" \\\"$inputs\\\" $extraMem $defaultExtraTime; }\" EXIT\nmemCpuMonitor.sh $flag $mem $min $memO $core &\n" >> $job

        #echo "touch $startFlag" >> $job
        echo -e "unset SLURM_CPU_BIND" >> $job
        if [[ "$program" == *.Checkpoint ]]; then
            echo "srun -n 1 $slurmAcc sh -e -c \"checkpoint.sh \\\"sh $slurmScript $slurmScriptParas && touch $succFlag\\\" $flag $mem $min $extraMem && adjustDownStreamJobs.sh $smartSlurmLogDir; \"" >> $job
        else
            
            if [ ! -z "$runningSingleJob" ]; then

                # snakemake should land here 
                echo "srun -n 1 $slurmAcc sh -e -o pipefail -c \"sh $slurmScript $slurmScriptParas; \" && touch $succFlag" >> $job  
            else 

                # this is never called? Because runAsPipeline usa --wrap ? 
                echo "srun -n 1 $slurmAcc sh -e -o pipefail -c \"sh $slurmScript $slurmScriptParas; \" && adjustDownStreamJobs.sh $smartSlurmLogDir && touch $succFlag" >> $job
            fi     
        fi
    fi
    #echo "kill -9 \$mypid" >> $job
    #echo "sleep 15 # wait slurm get the job status into its database" >> $job

    #echo "echo Job done. Summary:" >> $job

    # echo "sacct --format=JobID,Submit,Start,End,MaxRSS,State,NodeList%30,Partition,ReqTRES%30,TotalCPU,Elapsed%14,Timelimit%14 --units=M -j \$SLURM_JOBID" >> $job
    # #echo SLURM_JOBID=\$SLURM_JOBID >> $job

    # echo "emailAndRecord.sh \"$program\" \"${ref//\//-}\" \"$flag\" \"$inputSize\" \"$core\" \"$memO\" \"$timeO\" \"$mem\" \"$time\"" >> $job #  >/dev/null" >> $job

    # echo "adjustDownStreamJobs.sh $smartSlurmLogDir/" >> $job

    # echo "[ -f $succFlag ] ||  { touch $failFlag; exit 1; }" >> $job

    #[ -z "$dir" ] && dir="./" || mkdir -p $dir

    ## todo: need fix this when path not starting with /
    [ -z "$out" ] || ( [[ "$out" == /* ]] && mkdir -p $(dirname $out) || mkdir -p $dir/$(dirname $out) )
    [ -z "$err" ] || ( [[ "$err" == /* ]] && mkdir -p $(dirname $err) || mkdir -p $dir/$(dirname $err) )

    # [ -f $outFlag ] && [ ! -z "$out" ] && [[ ! $outFlag == "$out" ]] && echo ln -s $outFlag $out >> $job
    # [ -f $errFlag ] && [ ! -z "$err" ] && [[ ! $errFlag == "$err" ]] && echo ln -s $errFlag $err >> $job


    #if [ -z "$smartSlurmLogDir/" ]; then
    #    cmd="/usr/bin/sbatch --mail-type=FAIL --requeue --parsable -p $partition --mem $mem -t $time --open-mode=append $slurmAcc"
    #else
    [ -z "$err" ] && errFlag=$outFlag
    #[ -z "$deps" ] || deps="-H $deps" # if depends on other jobs, hold the job
    
        # if [ -z "$deps" ]; then #errFlag=$smartSlurmLogDir/$flag.err
        #     cmd="/usr/bin/sbatch --mail-type=FAIL --requeue --parsable -p $partition -c $core --mem $mem -t $time --open-mode=append -o $outFlag -e $outFlag -J $name $slurmAcc"
        #     echo -ne "\n#Command used to submit the job: /usr/bin/sbatch --mail-type=FAIL --requeue --parsable -p \$myPartition -c $core --mem \$myMem -t \$myTime --open-mode=append -o $outFlag -e $outFlag -J $name $deps $slurmAcc" >> $job
        # else
            cmd="/usr/bin/sbatch --mail-type=FAIL --requeue --parsable -p $partition -c $core --mem $mem -t $time --open-mode=append -o $outFlag -e $errFlag -J $name $deps $chdir $onhold $slurmAcc"
            echo -ne "\n#Command used to submit the job: /usr/bin/sbatch --mail-type=FAIL --requeue --parsable -p \$myPartition -c $core --mem \$myMem -t \$myTime --open-mode=append -o $outFlag -e $errFlag -J $name $deps $chdir $slurmAcc" >> $job
        #fi

        #rm $startFile $failFile $flag.failed $flag.killed # 2>/dev/null  || :
    #fi
    if [ -z "$slurmScript" ]; then
        cmd="$cmd $additionalPara $job"
        echo -e " $additionalPara $job" >> $job
    else
        cmd="$cmd $additionalPara $job $slurmScriptParas"
        echo -e " $additionalPara $job $slurmScriptParas" >> $job
    fi



    #echo -e "\n#Command used to submit the job: $cmd" >> $job

    #echoerr $resAjust
fi 

echoerr New sbatch command to submit job:
echoerr $cmd
if [[ "$testRun" == "run" ]]; then
    echoerr Start submtting job...
    #output=`$cmd`
    #echoerr $output
    #jobID=`$cmd`

     #set -x

    # try to submit the job three time in case it fails
    for i in {1..3}; do
        jobID=$($cmd 2>jobErr)
        jobErr=$(< jobErr); rm jobErr

        #catch jobID jobErr '$cmd'
        #echoerr id: .$jobID.
        #echoerr err: .$jobErr.
        [ -z "$jobErr" ] || echo $jobErr >&2 

        # no error and a number job ID
        if [ ! -z "$jobID" ] && [ "$jobID" -eq "$jobID" ]; then
            touch $outFlag.$jobID
            break
        elif [ ! -z "$jobErr" ]; then
            echoerr "Got error $jobErr"

            # job submitted but the jobID is empty
            if [[ "jobErr" == *"Socket timed out"* ]]; then
                sleep 2
                jobID=`sacct --name $flag --format jobid | head -n 1`
                touch $outFlag.$jobID
                break
            else
                echoerr "not sure what error it is: $jobErr"
                exit
            fi
        fi
        sleep 2
    done

    #echoerr final output is: $jobID

    #[ ! -z "$deps" ] && [[ "$jobID" =~ ^[0-9]+$ ]] && scontrol hold $jobID # todo: move this to sbatch command using -H

    #set +x

    #touch $submitFlag
else
    jobID=$(date +"%4N")
    #echoerr This is a testing, not really running a job...
    echoerr "This is testing, so no job is submitted. In real run it should submit job such as: Submitted batch job $jobID"
    #jobID=
fi

echo -e "\n#myMem=$mem myTime=$time" >> $job

# add this to the job script
# echo -e "\n#Sbatch command output:\n#Submitted batch job $jobID" >> $job

# if [ ! -z "$runningSingleJob" ]; then
#     if [ -z "$out" ]; then
#         touch $outFlag && ln -s $outFlag slurm-$jobID.out
#     else
#         out=${out/\%j/$jobID};
#         touch $outFlag && ln -s $outFlag $out
#     fi
#     if [ ! -z "$err" ]; then
#         err=${err/\%j/$jobID};
#         touch $errFlag && ln -s $errFlag $err
#     fi
# fi

#[ ! -z "$err" ] && ln -s $outFlag $err
#if [ -z "$smartSlurmLogDir/" ]; then
#    out=${out/\%j/$jobID};
#    ln -s $job ${out/.out/}.sh
#else

#if [[ "$parentCmd" != */bin/snakemake* ]]; then
    printf "%-10s  %-20s  %-10s %-10s %-10s %-10s %-15s\n" $jobID $depsID $flag $program $ref $inputs>> $smartSlurmLogDir/allJobs.txt;
#fi

#echoerr New slurmScirpt is ready. The content is:
#cat $job >&2


echo Submitted batch job $jobID
echoerr -e "Submitted batch job $jobID" 

[ -z "$runningSingleJob" ] || echo -e "With mem: $mem and time: $time\nLog is appended .smartSlurm.log \nPlease search $current_date to find it.\n" >&2


