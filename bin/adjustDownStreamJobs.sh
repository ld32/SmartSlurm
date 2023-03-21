#!/bin/sh

#set -x 

Usage="Usage: $0 full_path_to_flag_folder  flag(job name) \n  Note: this script will go through job id list file, find the downstream jobs, and return them as a string of job flags. "

echo 

extraMin=20
extraMem=500

x=`realpath $0`
. ${x%\/bin\/adjustDownStreamJobs.sh}/config/partitions.txt || { echo Partition list file not found: partition.txt; exit 1; }

echo Running: $0  $@

path=$1

[ -z "$2" ] && { echo -e "$Usage"; exit; }

[ -f $1/allJobs.txt ] || { echo -e "job id file $path/allJobs.txt does not exist\n$Usage"; exit 1; }

# jobid, deps, flag, software, ref, input, inputSize
text=`cat $path/allJobs.txt`
 
job=$2

IFS=$' ';  

# check the third column for the job name, then find the the job id in column 1
id=`echo $text | awk '{if ($3 ~ /'"$job/"') print $1;}' | tail -n 1`

echo -e "Find current job id (flag: $job):\n$id"
[ -z "$id" ] && { echo -e "job id for job name $job not found in $path/allJobs.txt!"; exit; }

echo 

idNames=`echo $text | awk '{if ($2 ~ /'"$id/"') print $2, $3;}'`

echo -e "Jobs on the save dependency level with current job:\n$idNames"
[ -z "$idNames" ] && { echo -e "Downstream job ids not found for $id"; exit; }

# sleep for a random number of seconds to avoid race condition
#sleep $(( ( RANDOM % 120)  + 1 ))

IFS=$'\n';
for i in $idNames; do
    echo 1working on $i
    id1=${i% *}; name=${i#* };
    allDone=""
    IFS=$' '; 
    for j in ${id1//\./ }; do 
        echo 2working on $j
        [[ "$j" == "$id" ]] && echo Ignore. It is the current job. It should adjust the mem and time for the downsteam job.  && continue
        
        echo look for the job flag for $j
        job=`echo $text | awk '{if ($1 ~ /'"$j/"') print $3;}'`
        #[ -z "$job" ] && { echo -e "job name not found!"; exit; }
        [ -f "$path/$job.success" ] && echo -e This job was done! || { echo This job is not done yet: $job; allDone=no; }         
    done
    if [ -z "$allDone" ]; then
        echo Dependants for $name are all done except for the current job. Ready to adjust mem/runtime
        
        # look for the downstream job info:                            jobID,software, ref, inputs
        output=`cat $path/allJobs.txt | awk '{if ($3 ~ /'"$name/"') print $1, $4, $5, $6;}'`
    
        IFS=' ' read -a arrIN <<< "$output"
    
        id2=${arrIN[0]}
        software=${arrIN[1]}   
        ref=${arrIN[2]}
        inputs=${arrIN[3]}    
        
        mem="" 
        [[ "$inputs" == "none" ]] && exit
        
        inputSize=`{ du --apparent-size -c -L ${inputs//,/ } 2>/dev/null || echo notExist; } | tail -n 1 | cut -f 1`
        if [ -f ~/smartSlurm/stats/$software.${ref//\//-}.mem.stat ]; then    
            output=`estimateMemTime.sh $software ${ref//\//-} $inputSize`
            echo "Output from estimateMemTime.sh: $output"
            if [[ "$output" == "outOfRange" ]]; then 
                echo "Input size is too big for the curve to estimate! Use default mem and runtime to submit job."
                # not deleting mem.stat, so other jobs will not re-build it within 60 minutes
            elif [ ! -z "$output" ]; then
                mem=$((${output% *}+extraMem)); time=$((${output#* }+extraMin));     
                echo "Give ${extraMem}M extra memory and $extraMin more minutes. \nSo use this to submit the job: $mem M  ${time}m"
            fi 
        
        fi
        
        if [[ "$output" == "outOfRange" ]] && test `find ~/smartSlurm/stats/$software.${ref//\//-}.mem.stat -mmin +60` || [ ! -f ~/smartSlurm/stats/$software.${ref//\//-}.mem.stat ]; then  
            echo "Do not have a formula, or it is old and out of range. Let us build one..."   

            #[ test `find ~/smartSlurm/jobRecord.txt -mmin -20` ] && echo jobRecord.txt synced within 20 hour. No need to re-sync || cat $HOME/smartSlurm/myJobRecord.txt > ~/smartSlurm/jobRecord.txt  
            
            #jobStatistics.sh $software ${ref//\//-} 4 1>&2 
            
            OUT="$(mktemp -d)"

            #filter by software and reference
            # todo: maybe able to replace / in ref at begaining of the script?
            ref1=${ref//\//-}
            grep COMPLETED ~/smartSlurm/jobRecord.txt | awk -F"," -v a=$software -v b=$ref1 '{ if($12 == a && $13 == b) {print $2, $7 }}' | sort -r  -k1,1 -k2,2 | sort -u -k1,1 > $OUT/mem.txt
            grep COMPLETED ~/smartSlurm/jobRecord.txt | awk -F"," -v a=$software -v b=$ref1 '{ if($12 == a && $13 == b) {print $2, $8 }}' | sort -r  -k1,1 -k2,2 | sort -u -k1,1 > $OUT/time.txt

            echo "Got mem data from jobRecord.txt (content of mem.txt):"
            cat $OUT/mem.txt 

            echo "Got time data from jobRecord.txt (content of time.txt):"
            cat $OUT/time.txt 

            if [[ $(wc -l <$OUT/mem.txt) -lt 3 ]]; then
                echo There are less than 3 records. No way to fit a curve. Exiting...
                echo There are less than 3 records. No way to fit a curve.
                 
            else

                cd $OUT

                # make plot and calculate statistics
                gnuplot -e 'set term pdf; set output "mem.pdf"; set title "Input Size vs. Memory Usage" font "Helvetica Bold,18"; set xlabel "Input Size(K)"; set ylabel "Memory Usage(M)"; f(x)=a*x+b; fit f(x) "mem.txt" u 1:2 via a, b; t(a,b)=sprintf("f(x) = %.2fx + %.2f", a, b); plot "mem.txt" u 1:2,f(x) t t(a,b); print "Finala=", a; print "Finalb=",b; stats "mem.txt" u 1 ' 2>&1 | grep 'Final\| M' | awk 'NF<4{print $1, $2}' |sed 's/:/=/' | sed 's/ //g' > mem.stat.txt

                # make plot and calculate statistics
                gnuplot -e 'set term pdf; set output "time.pdf"; set title "Input Size vs. Time Usage" font "Helvetica Bold,18"; set xlabel "Input Size(K)"; set ylabel "Time(Min)"; f(x)=a*x+b; fit f(x) "time.txt" u 1:2 via a, b; t(a,b)=sprintf("f(x) = %.2fx + %.2f", a, b); plot "time.txt" u 1:2,f(x) t t(a,b); print "Finala=", a; print "Finalb=",b; stats "time.txt" u 1 ' 2>&1 | grep 'Final\| M' | awk 'NF<4{print $1, $2}' |sed 's/:/=/' | sed 's/ //g' > time.stat.txt
                
                mv $OUT/mem.stat.txt ~/smartSlurm/stats/$software.$ref.mem.stat 
                mv $OUT/time.stat.txt ~/smartSlurm/stats/$software.$ref.time.stat  
                echo There are more than $ref jobs already run for this software, statics is ready for current job: 
                echo Memeory statisics:
                echo "inputsize mem(M)"
                cat ~/smartSlurm/stats/$software.$ref.mem.stat
                echo
                echo Time statistics:
                echo "inputsize time(minute)"
                cat ~/smartSlurm/stats/$software.$ref.time.stat
 
                mv $OUT/mem.txt ~/smartSlurm/stats/$software.$ref.mem.txt
                mv $OUT/time.txt ~/smartSlurm/stats/$software.$ref.time.txt

                convert $OUT/mem.pdf -background White -flatten ~/smartSlurm/stats/$software.$ref.mem.pdf
                convert $OUT/time.pdf -background White -flatten ~/smartSlurm/stats/$software.$ref.time.pdf
                pdftoppm ~/smartSlurm/stats/$software.$ref.mem.pdf  -png > ~/smartSlurm/stats/$software.$ref.mem.png
                pdftoppm ~/smartSlurm/stats/$software.$ref.time.pdf  -png > ~/smartSlurm/stats/$software.$ref.time.png

                echo
                echo You can see the plot using commands:
                echo display ~/smartSlurm/stats/$software.$ref.mem.pdf
                echo display ~/smartSlurm/stats/$software.$ref.time.pdf
                
                cd -

                rm -r $OUT 2>/dev/null

                echo got files in ~/smartSlurm/stats:  
                ls -lrt ~/smartSlurm/stats
                
                if [ -f ~/smartSlurm/stats/$software.${ref//\//-}.mem.stat ]; then    
                    output=`estimateMemTime.sh $software ${ref//\//-} $inputSize`
                    echo "Output from estimateMemTime.sh: $output"
                    if [[ "$output" == "outOfRange" ]]; then 
                        echo Input size is too big for the curve to estimate! Use default mem and runtime to submit job.
                        # not deleting mem.stat, so other jobs will not re-build it within 60 minutes
                    elif [ ! -z "$output" ]; then
                        mem=$((${output% *}+extraMem)); time=$((${output#* }+extraMin));     
                        echo Give ${extraMem}M extra memory and $extraMin more minutes. 
                    fi 
                fi        
            fi
        fi
    
        [ -z "$mem" ] && continue
        
        echo Got estimation inputsize: $inputSize mem: $mem  time: $time
        hours=$((($time + 59) / 60))

        echo looking partition for hour: $hours
        
        setPartition $hours partition
        
        seconds=$(($time * 60))
        time=`eval "echo $(date -ud "@$seconds" +'$((%s/3600/24))-%H:%M:%S')"`
        scontrol show job $id2                
        echo running: scontrol update jobid=$id2 timelimit=$time partition=$partition MinMemoryNode=${mem}
        scontrol update JobId=$id2 TimeLimit=$time Partition=$partition  MinMemoryNode=${mem}
        scontrol show job $id2

        #echo "scontrol update JobId=$id2 TimeLimit=$time Partition=$partition  MinMemoryNode=${mem}" >> $path/$name.sh
    else 
        echo Need wait for other jobs to finish before we can ajust mem and runtime...
    fi           
done
    
