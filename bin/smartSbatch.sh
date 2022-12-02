#!/bin/sh

# features:
# Auto adjust partition according to run-time request if they does not match up
# Auto check if slurm script exists
# Auto create output and erro folders if not exist
# Auto adjust memory and run time based on earlier mem time usage 
# Support re-run checking
# Auto requeue when out of time
# Auto requeue when out of memory
# Keep good log and send informative email

# Auto adjust partition according to run-time request if they does not match up
# Auto adjust memory and run time based on earlier mem time usage
# Support re-run from breaking-point
# Auto requeue job with double time when out of time
# Auto requeue job with double memory when out of memory

# todo:
# 
#

#set -x 

echoerr() { echo "$@" 1>&2; }

usage() { echoerr -e "Short Usage: \n${0##*/} [-P logDir]  <-S software> [-R reference] [-F uniqueJobFlag] [-I inputList] [sbatch options] [run]\nDetail Detail Usage:\n${0##*/} [-P logDir, optional. Such as: ./] ,-S software, required. Such as: bowtie2-4core> [-R reference, optional. Such as: hg19] [-F uniqueJobFlag, optional. Such as 1.bowtie.s1.fa] [-I inputFileOrFolderList, optional. Such as: read1.fq,read2.fq] <regular sbatch options, optional. Such as: job.sh or -p short -c 1 -t 2:0:0 --mem 2G --wrap \"my_application para1 para2\"> [run, optional: run will submit job, empty will do a dry run without submitting a job.]"; exit 1; }

[ -z "$1" ] && usage

[[ "-h" == "$1" ]] && usage

echoerr -e "Running:"

cmd="smartSbatch"
whitespace="[[:space:]]"
for i in "$@"; do
    if [[ $i =~ $whitespace ]]; then
        i="\"$i\""
    fi
    cmd="$cmd  $i"
done
echoerr $cmd

testRun=${@: -1} 

echoerr 
array=( "$@" )

# get the 6 parameters for ssbatch
for (( i=0; i<$(($#)); i++ )); do
    [ -z "${array[$i]}" ] && continue
  	echoerr $i " / " $(($#)) " : " ${array[$i]}
  	case "${array[$i]}" in
  	    "-P"            )   echoerr found -P && logDir="${array[$i+1]}" && array[$i]="" && array[$i+1]="";; 
  		--LogDir=*  )   echoerr found --LogDir && logDir="${array[$i]}" && logDir=${logDir/--LogDir=/} && array[$i]="";;
  		"-S" 			)   echoerr found -S && software="${array[$i+1]}" && array[$i]="" && array[$i+1]="";;
  		--Software=* 	)   echoerr Found --Software= && software="${array[$i]}" && software="${software/--Software=/}" && array[$i]="";;
  		"-R" 			)   echoerr Found -R && ref="${array[$i+1]}" && array[$i]="" && array[$i+1]="";;
  		--Ref=* 		)   echoerr Found --Ref= && ref="${array[$i]}" && ref="${ref/--Ref=/}" && array[$i]="";;
  		"-F" 			)   echoerr Found -F && flag="${array[$i+1]}" && array[$i]="" && array[$i+1]="";;
  		--Flag=* 	    )   echoerr Found --Flag= && flag="${array[$i]}" && flag="${flag/--Flag=/}" && array[$i]="";;
        "-I" 			)   echoerr Found -I && inputs="${array[$i+1]}" && array[$i]="" && array[$i+1]="";;
  		--Inputs=* 	    )   echoerr Found --Inputs= && inputs="${array[$i]}" && inputs="${inputs/--Inputs=/}" && array[$i]="";;
        "-D"            )   echoerr Found -D  && deps="${array[$i+1]}" && array[$i]="" && array[$i+1]="";;
        "--Dependency=" )   echoerr Found -D  && deps="${array[$i]}" && deps="${deps/--Dependency=/}" && array[$i]="";;
        *               )   echoerr Found sbatch options &&  sbatchLocation=$i && break
    esac
done 

# this is regular sbatch without smart options added
[ -z "$logDir$software$ref$flag$inputs$deps" ] && set -- "$@" run && testRun=run

echoerr options before sbatch options: slurmAcc: $slurmAcc logDir: $logDir software: $software ref: $ref jobFlag: $flag inputs: $inputs deps: $deps

#exit

depsO=$deps 

#echoerr -e "Parameters:\n${array[*]} " 
echoerr Parsing sbatch command... 
for (( i=$sbatchLocation; i<`[[ "$testRun" == "run" ]] && echo $(($# -1))  || echo $(($#))`; i++ )); do
#for (( i=$sbatchLocation; i< $(($#)); i++ )); do
    [ -z "${array[$i]}" ] && continue
  	echoerr $i " / " $(($#)) " : " ${array[$i]}
  	case "${array[$i]}" in
  	    "-p"            )   echoerr Found -p && partition="${array[$i+1]}" && array[$i]="" && array[$i+1]="";; # will set partition later
  		--partition=*   )   echoerr Found --partition && partition="${array[$i]}" && partition=${partition/--partition=/} && array[$i]="";;
  		--mem-per-cpu=* ) 	echoerr Found --mem-per-cpu= && mem1="${array[$i]}" && mem1="${mem1/--mem-per-cpu=/}";;
  		"-c" 	 		) 	echoerr Found -c && core="${array[$i+1]}";;
  		--cpus-per-task=*)  echoerr Found --cpus-per-task= && core="${array[$i]}" && core="${core/--cpus-per-task=/}";;
  		"-n" 	 		) 	echoerr Found -n && task="${array[$i+1]}";;
  		--ntasks=*      )   echoerr Found --ntasks= && task="${array[$i]}" && task="${task/--ntasks=/}";;
  		"-N" 	 		) 	echoerr Found -N && node="${array[$i+1]}" && array[$i]="" && array[$i+1]="";;
  		--nodes=*       )   echoerr Found --nodes= && node="${array[$i]}" && node="${node/--nodes=/}" && array[$i]="";;
  		"--mem"  		) 	echoerr Found --mem && mem="${array[$i+1]}"&& array[$i]="" && array[$i+1]="";; # will set memory later
  		--mem=* 		)  	echoerr Found --mem= && [ -z "$mem" ] && mem="${array[$i]}" && mem=${mem/--mem=/} && array[$i]="";;
  		"-t" 			)   echoerr Found time && time="${array[$i+1]}" && array[$i]="" && array[$i+1]="";;
  		--time=* 		)   echoerr Found --time= && time="${array[$i]}" && time="${time/--time=/}" && array[$i]="";;
  		"-o" 			)   echoerr Found -o && out="${array[$i+1]}" && array[$i]="" && array[$i+1]="";;
  		--output=* 		)   echoerr Found --output= && out="${array[$i]}" && out="${out/--oupput=/}" && array[$i]="";;
  		"-e" 			)   echoerr Found -o && err="${array[$i+1]}" && array[$i]="" && array[$i+1]="";;
  		--error=* 		)   echoerr Found --error= && err="${array[$i]}" && err="${err/--error=/}" && array[$i]="";;
        "-A"            )   echoerr Found -A  && [ -z "$slurmAcc" ] && slurmAcc="-A ${array[$i+1]}" && array[$i]="" && array[$i+1]="";;
        "-d" 			)   echoerr Found -d && dep="${array[$i+1]}";;
  		--dependency=* 	)   echoerr Found --dependency= && dep="${array[$i]}" && dep="${dep/--dependency=/}";;
  		
    esac
    #[ ! -z "$out" ] &&  CMDWithoutWrap="$CMDWithoutWrap -o $out" && CMDWithoutSlurmCMD="$CMDWithoutSlurmCMD -o $out"
    #[ ! -z "$err" ] &&  CMDWithoutWrap="$CMDWithoutWrap -e $err" && CMDWithoutSlurmCMD="$CMDWithoutSlurmCMD -e $err"
    [ -z "$wrapCMD" ] && [[ ${array[$i]} == "--wrap" ]] && echoerr Found --wrap ${array[$i+1]} && wrapCMD="${array[$i+1]}" && array[$i]="" && array[$i+1]="" 
  	[ -z "$wrapCMD" ] && [[ ${array[$i]} ==	--wrap=* ]] && echoerr Found --wrap= && wrapCMD="${array[$i]}" && wrapCMD="${wrapCMD/--wrap=}" && array[$i]="" || CMDWithoutWrap="$CMDWithoutWrap ${array[$i]}"
    [ -z "$slurmScript$wrapCMD" ] && [ -f "${array[$i]}" ] && echoerr Found slurmScript ${array[$i]} && slurmScript="${array[$i]}" && array[$i]=""
    [ -z "$slurmScript" ] && CMDWithoutSlurmCMD="$CMDWithoutSlurmCMD ${array[$i]}" || slurmScriptParas="$slurmScriptParas ${array[$i]}"
done 

echoerr 
echoerr Parsing result from sbatch commandline: 
echoerr sbatch options: time: $time mem: $mem mem-per-cpu: $mem1 task: $task core: $core node: $node out: $out err: $err dep: $dep

if [ -z "$slurmScript" ]; then
    echoerr wrapCMD: $wrapCMD 
    echoerr additional sbatch parameter: $CMDWithoutWrap 
else 
    echoerr additional sbatch parameter: $CMDWithoutSlurmCMD
    echoerr slurmScript: $slurmScript 
    echoerr slurmScriptParas: $slurmScriptParas
fi
echoerr test or run: $testRun

#stamp=`date -d "today" +"%Y.%m.%d.%H.%M-%N"`

if [ ! -z "$slurmScript" ]; then 
  echoerr Validating slurmScript: 
  firstRow=`head -n 1 $slurmScript`
  echoerr FirstRow of the script: $firstRow
  [[ "$firstRow" =~ ^#\!/bin/bash ]] || [[ "$firstRow" =~ ^#\!/usr/bin/bash ]] || [[ "$firstRow" =~ ^#\!/bin/sh ]] || [[ "$firstRow" =~ ^#\!/usr/bin/sh ]]|| { echoerr "Error: first row of the slurm script ($slurmScript) does not start with #!/bin/bash, #!/usr/bin/bash, #!/bin/sh, #!/usr/bin/sh";  exit 1;}
fi

[ -z "$wrapCMD" ] && [ -z "$slurmScript" ] && echoerr Error: Did not find --wrap, did not find slurmScript either. && exit 1 

if [[ ! -z "$slurmScript" ]]; then 
    echoerr 
    echoerr Parsing slurm script ... 
    while IFS=$'\n' read line; do
        echoerr Got line: $line
        arr=($line)
        if [[ ${arr[0]} == *SBATCH ]]; then 
            case "${arr[1]}" in 
                "-P"            )   echoerr found -P && logDir=${arr[2]};; 
                --LogDir=*  )   echoerr found --LogDir && logDir="${arr[1]}" && logDir="${logDir/--LogDir=/}";;
                "-S" 			)   echoerr found -S && software=${arr[2]};;
                --Software=* 	)   echoerr Found --Software= && software="${arr[1]}" && software="${software/--Software=/}";;
                "-R" 			)   echoerr Found -R && ref=${arr[2]};;
                --Ref=* 		)   echoerr Found --Ref= && ref="${arr[1]}" && ref="${ref/--Ref=/}";;
                "-F" 			)   echoerr Found -F && flag=${arr[2]};;
                --Flag=* 	    )   echoerr Found --Flag= && flag="${arr[1]}" && flag="${flag/--Flag=/}";;
                "-I" 			)   echoerr Found -I && inputs=${arr[2]};;
                --Inputs=* 	    )   echoerr Found --Inputs= && inputs="${arr[1]}" && inputs="${inputs/--Inputs=/}";;
                "-D"            )   echoerr Found -D  && deps=${arr[2]};;
                "--Dependency"  )   echoerr Found -D  && deps="${arr[1]}" && deps="${deps/--Dependency/}";;
                
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
                # "-o"    ) [ -z "$out" ] && echoerr Found -o && out="${arr[2]}" && out=${out/\%j/\$\{SLURM_JOBID\}} && out=${out/\%A/\$\{SLURM_ARRAY_JOB_ID\}} && out=${out/\%a/\$\{SLURM_ARRAY_TASK_ID\}};;
                # --output=* ) [ -z "$out" ] && echoerr Found --output= && out="${arr[1]}" && out="${out/--output=/}" && out=${out/\%j/\$\{SLURM_JOBID\}} && out=${out/\%A/\$\{SLURM_ARRAY_JOB_ID\}} && out=${out/\%a/\$\{SLURM_ARRAY_TASK_ID\}};;
                # "-e" ) [ -z "$err" ] && echoerr Found -e && err="${arr[2]}" && err=${err/\%j/\$\{SLURM_JOBID\}} && err=${err/\%A/\$\{SLURM_ARRAY_JOB_ID\}} && err=${err/\%a/\$\{SLURM_ARRAY_TASK_ID\}};;
                # --error=* ) [ -z "$err" ] && echoerr Found --error= && err="${arr[1]}" && err="${err/--error=/}" && err=${err/\%j/\$\{SLURM_JOBID\}}  && err=${err/\%A/\$\{SLURM_ARRAY_JOB_ID\}} && err=${err/\%a/\$\{SLURM_ARRAY_TASK_ID\}};;
                
                "-o"    ) [ -z "$out" ] && echoerr Found -o && out="${arr[2]}";;
                --output=* ) [ -z "$out" ] && echoerr Found --output= && out="${arr[1]}" && out="${out/--output=/}";;
                "-e" ) [ -z "$err" ] && echoerr Found -e && err="${arr[2]}";;
                --error=* ) [ -z "$err" ] && echoerr Found --error= && err="${arr[1]}" && err="${err/--error=/}";;
                
                "-d"    ) [ -z "$dep" ] && echoerr Found -d && dep="${arr[2]}";;
                --dependency=* ) [ -z "$dep" ] && echoerr Found --dependency= && dep="${arr[1]}" && dep"${dep/--dependency=/}";;
            esac
        fi     
    done < $slurmScript 


    echoerr
    echoerr Parsing result from slurm script:  
    echoerr time: $time mem: $mem mem-per-cpu: $mem1 task: $task core: $core node: $node out: $out err: $err dep: $dep

    echoerr slurmScript: $slurmScript 
    echoerr additional sbatch parameter: $CMDWithoutSlurmCMD
    echoerr slurmScriptParas: $slurmScriptParas
fi

# this is regular sbatch without smart options added
if [ -z "$logDir$software$ref$flag$inputs$deps" ]; then  
    software=regularSbatch
    tm=`mktemp XXXXXXXX --dry-run`
    if [ ! -z "$wrapCMD" ]; then
        arr=($$wrapCMD)
        # need remove the parameters, to get unique name for job types
        if [[ "$arr[0]" == sh || "$arr[0]" == bash || "$arr[0]" == python || "$arr[0]" == python3  || "$arr[0]" == matlab || "$arr[0]" == Rscript ]]; then
          ref=$arr[1]  
        else
          ref=$arr[0]  
        fi
        ref=${ref##*/}
        [ -z "$flag" ] && flag=${wrapCMD// /}.$tm
        
    else
        ref=${slurmScript##*/}
        [ -z "$flag" ] && flag=$ref.${slurmScriptParas// /}.$tm
    fi
     
fi    

[ -z "$logDir" ] && logDir=~/.smartSlurm
logDir=${logDir%/}
[ -z "$software" ] && { echoerr software is emtpty; usage; }
[ -z "$ref" ] && ref=none
[ -z "$inputs" ] && inputs=none
[ -z "$deps" ] && deps=null

[ -d "$logDir" ] || { echoerr Directory not exist: $logDir; usage; }


flag=${flag//\//-}

flagDir=$logDir/log
[ ! -d $flagDir ] && mkdir $flagDir 
job=$flagDir/$flag.sh
#submitFlag=$flagDir/$flag.submitted
succFlag=$flagDir/$flag.success
failFlag=$flagDir/$flag.failed
#startFlag=$flagDir/$flag.start
killFlag=$flagDir/$flag.user.killed
outFlag=$flagDir/$flag.out
errFlag=$flagDir/$flag.err

if [ -z "$mem" ]; then
    if [ -z "$mem1" ]; then
        mem=1G
    else
        [[ "$mem1" == *G ]] && mem=$(( ${mem1%G} * $core ))G || mem=$(( ${mem1%M} * $core ))M
    fi
fi

deps=${deps#.}

#echo d1 ${d/\./} >&2
if [[ $deps == null ]]; then
    deps=""
    echo depend on no job >&2
elif [[ $deps == ${deps/\./} ]]; then
    echo depend on single job >&2
    deps="--dependency=afterok:${deps/\./}"
else
    echo depend on multiple jobs >&2
    tmp=""
    for t in ${deps//\./ }; do
       echo working on $t >&2
        tmp="$tmp:$t"
    done
    [ ! -z "$tmp" ] && deps="--dependency=afterok$tmp"
fi

# check job done before or not
if [[ "$testRun" == "run" ]]; then 
    if [ -f $succFlag ]; then
        stepID=${flag%%.*}
        if ([ -f $flagDir/skipAllSuccessJobs.txt ] || [ -f $flagDir/skipAllSuccessJobs$stepID.txt ]) && [ -z "$deps" ]; then 
            exit    
        elif [[ "$te"  != test ]] && [ ! -f $flagDir/reRunAllSuccessJobs.txt ] && [ ! -f $flagDir/reRunAllSuccessJobs$stepID.txt ] && [ -z "$deps" ]; then
            stepName=${flag#*.}; stepName=${stepName#*.}; stepName=${stepName%%.*}
            echo $flag was done before, do you want to re-run it? >&2
            echo -e "y:        To re-run this job, press y, then enter key." >&2
            echo -e "ystep:    To re-run all jobs for step $stepID: $stepName, type ystep, then press enter key." >&2
            echo -e "yall:     To re-run all jobs, type yall, then press enter key." >&2
            echo -e "enter:    To not re-run this job, directly press enter key." >&2
            echo -e "nstep:    To not re-run all successful jobs for step $stepID: $stepName, type nstep, then press enter key." >&2
            echo -e "nall:     To not re-run all successful jobs, type nall, then press enter key." >&2
            read -p "" x </dev/tty

            if [[ "$x" == "y" ]]; then
                echo "Will re-run the down stream steps even if they are done before (because they have deps - see code in row 70)." >&2 
            elif [[ "$x" == "ystep" ]]; then 
                touch $flagDir/reRunAllSuccessJobs$stepID.txt
            elif [[ "$x" == "nstep" ]]; then
                touch $flagDir/skipAllSuccessJobs$stepID.txt
                exit
            elif [[ "$x" == "yall" ]]; then 
                touch $flagDir/reRunAllSuccessJobs.txt
            elif [[ "$x" == "nall" ]]; then
                touch $flagDir/skipAllSuccessJobs.txt
                exit    
            else
                exit
            fi  
        fi
        echoerr
    fi    
fi

rm $flagDir/$flag.*  2>/dev/null 

# original mem and time 
memO=$mem; timeO=$time
inputSize=0

echoerr 
echoerr Check if there input file list and this job does not depend on other jobs

# not depends on other job, and there is input file list
if [ -z "$dep" ] &&  [[ "$inputs" != "none" ]]; then
    inputSize=`{ du --apparent-size -c -L ${inputs//,/ } 2>/dev/null || echo notExist; } | tail -n 1 | cut -f 1`

    if [[ "$inputSize" == "notExist" ]]; then 
        echoerr Some or all input files not exist: $inputs  && exit 1
    else 
        inputSize=$(($inputSize/1024)); # convert to M      
        echoerr inputSize: $inputSize 
        #rm ~/.rcbio/$software.$ref.mem.stat.final # for testing 
        # there is a job record curve fit already?    
        
        if [ ! -f ~/.smartSlurm/$software.${ref//\//-}.mem.stat.final ]; then   
            mkdir -p ~/.smartSlurm
            echoerr Do not have a formula. Let us build one...
            jobStatistics.sh $software ${ref//\//-} 4 1>&2
            justRunStats=yes
        fi
        if [ -f ~/.smartSlurm/$software.${ref//\//-}.mem.stat.final ]; then    
            output=`estimateMemTime.sh $software ${ref//\//-} $inputSize`
            if [[ "$output" == "outOfRange" ]]; then 
                echoerr Input size is too big for the curve to estimate!
                if [ -z "$justRunStats" ]; then 
                    echoerr Delete the curve and try to re-run the statistics.
                    rm ~/.smartSlurm/$software.${ref//\//-}.mem.stat.final ~/.smartSlurm/$software.${ref//\//-}.time.stat.final ~/.smartSlurm/jobRecord.txt
                    jobStatistics.sh $software ${ref//\//-} 4 1>&2 
                    if [ -f ~/.smartSlurm/$software.${ref//\//-}.mem.stat.final ]; then    
                        output=`estimateMemTime.sh $software ${ref//\//-} $inputSize`
                        if [[ "$output" == "outOfRange" ]]; then 
                            echoerr Input size is too big for the curve to estimate! Use default mem and runtime to submit job.
                        else
                            mem=${output% *}M; time=${output#* };     
                            echoerr Got estimation inputsize: $inputSize mem: $mem  time: $time
                        fi 
                    else 
                        echoerr Building formula failed. Use default mem and runtime to submit job.
                    fi
                else 
                    echoerr Use default mem and runtime to submit job.    
                fi    
            else
                mem=${output% *}M; time=${output#* };     
                echoerr Got estimation inputsize: $inputSize mem: $mem  time: $time
            fi 
        else 
            echoerr Building formula failed. Use default mem and runtime to submit job.
        fi
    fi
# do not have input file list    #todo: we can try to figure out the input from the software parameters here
elif [[ "$inputs" == "none" ]]; then
    ref=${ref//\//-}
    if [ ! -f ~/.smartSlurm/$software.$ref.mem.stat.noInput ] || [ -z "`find ~/.smartSlurm/$software.$ref.mem.stat.noInput -mmin -60`" ]; then  
    
        #cat /home/*/.smartSlurm/myJobRecord.txt > ~/.smartSlurm/jobRecord.txt 
        cat $HOME/.smartSlurm/myJobRecord.txt > ~/.smartSlurm/jobRecord.txt
        #filter by software and reference
        # grep COMPLETED ~/.smartSlurm/jobRecord.txt | awk -F, -v a=$software -v b=$ref '{ if($3 == a && $4 == b) {print $13 }}' | sort -rn > ~/.smartSlurm/$software.$ref.mem.stat.noInput
        # grep COMPLETED ~/.smartSlurm/jobRecord.txt | awk -F, -v a=$software -v b=$ref '{ if($3 == a && $4 == b) {print $14 }}' | sort -rn > ~/.smartSlurm/$software.$ref.time.stat.noInput
        
        grep COMPLETED ~/.smartSlurm/jobRecord.txt | awk -F, -v a=$software -v b=$ref '{ if($11 == a && $12 == b) {print $7 }}' | sort -rn > ~/.smartSlurm/$software.$ref.mem.stat.noInput
        grep COMPLETED ~/.smartSlurm/jobRecord.txt | awk -F, -v a=$software -v b=$ref '{ if($11 == a && $12 == b) {print $8 }}' | sort -rn > ~/.smartSlurm/$software.$ref.time.stat.noInput
        [ -s ~/.smartSlurm/$software.$ref.mem.stat.noInput ] || rm ~/.smartSlurm/$software.$ref.mem.stat.noInput        
        [ -s ~/.smartSlurm/$software.$ref.time.stat.noInput ] || rm ~/.smartSlurm/$software.$ref.time.stat.noInput
    fi
    if [ -f ~/.smartSlurm/$software.$ref.mem.stat.noInput ]; then 
        totalRow=`wc -l < ~/.smartSlurm/$software.$ref.mem.stat.noInput`
        if [[ ! -z "$totalRow" && $totalRow -ge 5 ]]; then 
            cutoffRow=$(( totalRow*10/100 )) # top 10
            
            [[ $cutoffRow == 0 ]] && cutoffRow=1
            
            mem=`head -n $cutoffRow ~/.smartSlurm/$software.$ref.mem.stat.noInput | tail -n 1`; mem=${mem/\.*/}M
            time=`head -n $cutoffRow ~/.smartSlurm/$software.$ref.time.stat.noInput | tail -n 1`        
        else 
            echoerr There is less than 5 records. No way to fit a curve. Exiting...

        fi
    fi
fi

echoerr

[ -z "$time" ] && { echoerr did not find time limit; exit 1; }
 
[ -z "$mem" ] && { echoerr did not find mem limit; exit 1; } 

[ -z "$core" ] && { [ ! -z "$task" ] && core="$task" || core=1; } 

[[ "$time" == *-* ]] && { day=${time%-*}; tem=${time#*-}; hour=${tem%%:*}; min=${tem#*:}; min=${min%%:*}; sec=${tem#$hour:$min}; sec=${sec#:}; } || { [[ "$time" =~ ^[0-9]+$ ]] && min=$time || { sec=${time##*:}; min=${time%:*}; min=${min##*:}; hour=${time%$min:$sec}; hour=${hour%:}; day=0;} }

[ -z "$day" ] && day=0; [ -z "$hour" ] && hour=0; [ -z "$min" ] && min=0;[ -z "$sec" ] && sec=0

echoerr $day $day,  $hour hour,  $min min,  $sec sec

# how many hours for sbatch command
hours=$(($day * 24 + $hour + ($min + 59 + ($sec + 59) / 60 ) / 60))

#echoerr looking partition for hour: $hours 
x=`realpath $0` 
. ${x%\/bin\/smartSbatch.sh}/config/partitions.txt || { echoerr Partition list file not found: partition.txt; exit 1; }

partition=`adjustPartition $hours $partition`

time=$day-$hour:$min:$sec

# 10 minutes less than the time in sbatch command
#seconds=$(($day * 24 * 60 * 60 + $hour * 60 * 60  + $min * 60 + $sec - 600))

#[ $seconds -le 60 ] && time=11:0 && seconds=60

#echoerr srun seconds: $seconds

#timeN=`eval "echo $(date -ud "@$seconds" +'$((%s/3600/24))-%H:%M:%S')"`

#echoerr New time for srun: $timeN

#echoerr 

# mem for srun is 10M less than sbatch command
#[[ "$mem" == *G ]] && memN=$(( 1024 * ${mem%G} -10 ))M || memN=$(( ${mem%M} -10 ))M

#[ ${memN%M} -le 1 ] && { echoerr Error: --mem for sbatch command should bigger than 11M; usage; }
#[ ${memN%M} -le 1 ] && mem=11M && memN=10M

#echoerr Mew mem for srun: $memN

#[ -z "$out" ] && out="slurm-\$SLURM_JOBID.out"
#[ -z "$err" ] && err="slurm-\$SLURM_JOBID.err"

echoerr 
echoerr Building new sbatch command ...
#echo "#!/bin/bash" > $job 
echo -e "#!/bin/bash\ntrap \"{ cleanUp.sh \"$logDir\" \"$software\" \"${ref//\//-}\" \"$flag\" \"$inputSize\" \"$core\" \"$memO\" \"$timeO\" \"$mem\" \"$time\" \"$partition\"  \"${slurmAcc#*-A }\" \\\"$0 $@\\\"; }\" EXIT" > $job

if [ -z "$slurmScript" ]; then 
    #echo "touch $startFlag" >> $job
    wrapCMD=`echo $wrapCMD | xargs echo -n`; wrapCMD=${wrapCMD%;} # remove ending space and ; from command
    echo "srun -n 1 $slurmAcc bash -e -c \"{ $wrapCMD; } && touch $succFlag\"" >> $job 
else 
    grep "^#SBATCH" $slurmScript >> $job || echo >> $job
    #echo "touch $startFlag" >> $job
    echo "srun -n 1 $slurmAcc bash -e -c \"{ sh $slurmScript $slurmScriptParas; } && touch $succFlag\"" >> $job 
fi
#echo "sleep 15 # wait slurm get the job status into its database" >> $job 

#echo "echo Job done. Summary:" >> $job 

# echo "sacct --format=JobID,Submit,Start,End,MaxRSS,State,NodeList%30,Partition,ReqTRES%30,TotalCPU,Elapsed%14,Timelimit%14 --units=M -j \$SLURM_JOBID" >> $job 
# #echo SLURM_JOBID=\$SLURM_JOBID >> $job

# echo "emailAndRecord.sh \"$software\" \"${ref//\//-}\" \"$flag\" \"$inputSize\" \"$core\" \"$memO\" \"$timeO\" \"$mem\" \"$time\"" >> $job #  >/dev/null" >> $job 

# echo "adjustDownStreamJobs.sh $flagDir $flag" >> $job 

# echo "[ -f $succFlag ] ||  { touch $failFlag; exit 1; }" >> $job

#[ -z "$dir" ] && dir="./" || mkdir -p $dir

## todo: need fix this when path not starting with /
[ -z "$out" ] || ( [[ "$out" == /* ]] && mkdir -p $(dirname $out) || mkdir -p $dir/$(dirname $out) )
[ -z "$err" ] || ( [[ "$err" == /* ]] && mkdir -p $(dirname $err) || mkdir -p $dir/$(dirname $err) )

# [ -f $outFlag ] && [ ! -z "$out" ] && [[ ! $outFlag == "$out" ]] && echo ln -s $outFlag $out >> $job
# [ -f $errFlag ] && [ ! -z "$err" ] && [[ ! $errFlag == "$err" ]] && echo ln -s $errFlag $err >> $job

echoerr New slurmScirpt is ready. The content is:
cat $job 1>&2

cmd="/usr/bin/sbatch --mail-type=FAIL --requeue --parsable -p $partition --mem $mem -t $time --open-mode=append -o $outFlag -e $errFlag -J $flag $deps $slurmAcc" 
if [ -z "$slurmScript" ]; then 
    cmd="$cmd $CMDWithoutWrap $job" 
else 
    cmd="$cmd $CMDWithoutSlurmCMD $job $slurmScriptParas" 
fi

echo -e "\n#Command used to submit the job:" >> $job
echo "#$cmd" >> $job 

rm $startFile $failFile $flag.failed $flag.killed  2>/dev/null  || echoerr 

echoerr New sbatch command to submit job: 
echoerr $cmd 
if [[ "$testRun" == "run" ]]; then
    echoerr Start submtting job...
    jobID=`$cmd 2>&1`
    #touch $submitFlag 
else 
    echoerr This is a testing, not really running a job...
    jobID="123"
fi    

out=${out/\%j/$jobID}
err=${err/\%j/$jobID}

# todo: including job id in flag and use link for flag.out and flag.err? Might need move this to cleanUp.sh
[ -z "$out" ] && out=slurm-$jobID.out  
[ -z "$err" ] && err=slurm-$jobID.err

ln -s $outFlag $out
ln -s $errFlag $err

ln -s $outFlag $outFlag.id.$jobID
ln -s $errFlag $errFlag.id.$jobID

# add this to the job script
echo -e "\n#Sbatch command output:\n#Submitted batch job $jobID" >> $job 

# print out on screen
echoerr -e "Sbatch command output:\nSubmitted batch job $jobID"

# this is for the clalling softwar to get the job id
echo  $jobID

#dep=${dep#afterok}; dep=${dep//:/.}; [ -z "$dep" ] && dep=null 

# for check dependency, rerun, kill downstream jobs, estimate memmory and run time
#[[ "$jobID" =~ ^[0-9]+$ ]]  &&  echo $jobID $dep $flag $software $ref $inputs >> $flagDir/allJobs.txt

printf "%-10s  %-20s  %-10s\n" $jobID $depsO $flag >> $flagDir/allJobs.txt
