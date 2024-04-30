#!/bin/sh

#set -x

Usage="Usage: $0 full_path_to_flag_folder \n  Note: this script will go through job id list file, find the downstream jobs, and return them as a string of job flags. "

echo
#for i in {1..200}; do sleep 1; echo adjusting $i; done & 

echo Running: $0  $@

smartSlurmLogDir=$1

[ -f $smartSlurmLogDir/allJobs.txt ] || { echo -e "job id file $smartSlurmLogDir/allJobs.txt does not exist\n$Usage"; exit 1; }

# jobid, deps, flag, software, ref, input, inputSize
text=`cat $smartSlurmLogDir/allJobs.txt`

#job=$2

IFS=$' ';

# check the third column for the job name, then find the the job id in column 1
#id=`echo $text | awk '{if ($3 ~ /'"$job/"') print $1;}' | tail -n 1`

#echo -e "Find current job id (flag: $job):\n$id"
#[ -z "$id" ] && { echo -e "job id for job name $job not found in $smartSlurmLogDir/allJobs.txt!"; exit; }

#echo

# while true; do
#     if `mkdir $smartSlurmLogDir/downsteamjob.adjusting  2>/dev/null`; then
#         break
#     fi
#     echo waiting for lock to adjust downsteam jobs 
#     sleep 1
# done 

# directly get id, deps, software, ref, and input here, if input is none, directly skip this job
output=`echo $text | awk '{if ($2 ~ /'"$SLURM_JOBID/"') print $1, $2, $3, $4, $5, $6;}'`

echo -e "Jobs on the same dependency level with current job:\n$output"
[ -z "$output" ] && { echo -e "Downstream job ids not found for $SLURM_JOBID"; exit; }

IFS=$'\n';
for i in $output; do
    echo 1working on $i
    IFS=' ' read -a arrIN <<< "$i"

    id=${arrIN[0]}
    deps=${arrIN[1]}
    name=${arrIN[2]}
    software=${arrIN[3]}
    ref=${arrIN[4]}  ; ref=${ref//\//-}
    inputs=${arrIN[5]}

    #[ -f $smartSlurmJobRecordDir/stats/extraMem.$software.$ref ] && extraMem=`sort $smartSlurmJobRecordDir/stats/extraMem.$software.$ref | tail -n1`
    #[ -f $smartSlurmJobRecordDir/stats/extraMem.$software.$ref ] && maxExtra=`sort -n $smartSlurmJobRecordDir/stats/extraMem.$software.$ref | tail -n1 | cut -d' ' -f1` && oomCount=`wc -l $smartSlurmJobRecordDir/stats/extraMem.$software.$ref | cut -d' ' -f1` && extraMem=$(( $maxExtra * $oomCount ))

    [ -f $smartSlurmJobRecordDir/stats/extraMem.$software.$ref ] && maxExtra=`sort -n $smartSlurmJobRecordDir/stats/extraMem.$software.$ref | tail -n1 | cut -d' ' -f1` && extraMem=$(( $maxExtra * 2 )) || extraMem=$(( $defaultExtraMem * 2 ))

    allDone=""
    IFS=$' ';
    for j in ${deps//:/ }; do
        echo 2working on $j
        [[ "$j" == "$SLURM_JOBID" ]] && continue; # ignore current job
         echo look for the job flag for $j
        job=`echo $text | awk '{if ($1 ~ /'"$j/"') print $3;}'`
        #[ -z "$job" ] && { echo -e "job name not found!"; exit; }
        [ -f "$smartSlurmLogDir/$job.success" ] && echo This job was done! $job || { echo This job is not done yet: $job; allDone=no; break;}
    done

    if [ -z "$allDone" ]; then
        date
        # todo: even this is no input, we may need to modify the runtime becaue we might have new stats from jobs finished after the job is submitted.
        if [[ "$inputs" == "none" ]]; then
            echo No input, do not need to adjust, directly release and run.
            scontrol release $id
            continue
        fi

        if [ -f $smartSlurmLogDir/$name.adjust ]; then
            echo Already adjusted? Directly release and run.
            scontrol release $id
            continue
        fi

        #ls -lrt $smartSlurmLogDir
        echo Dependants for $name are all done. Ready to adjust mem/runtime...

        echo -e "Re-adjust resource by upsteam job job $SLURM_JOB_ID:" >> $smartSlurmLogDir/$name.out
        grep ^$SLURM_JOB_ID $smartSlurmLogDir/allJobs.txt | awk '{print $1,  $2,  $3}' >> $smartSlurmLogDir/$name.out

        inputSize=`{ du --apparent-size -c -L ${inputs//,/ } 2>/dev/null || echo notExist; } | tail -n 1 | cut -f 1`

        if [[ "$inputSize" == "notExist" ]]; then
            scancel $id
            pwd >> $smartSlurmLogDir/$name.out
            echo One or multiple inputs are missing for this job. Cancelling it... >> $smartSlurmLogDir/$name.out
            echo -e "inputs:.${inputs}.-.${inputs//,/}." >> $smartSlurmLogDir/$name.out
            echo ${inputs//,/ } >> $smartSlurmLogDir/$name.out
            touch $smartSlurmLogDir/$name.missingInnput.has.to.cancel
            toSend=`cat $smartSlurmLogDir/$name.out`
            s="Cancel:$id:MissingInput:${inputs//,/ }"
            echo -e "$toSend" | mail -s "$s" $USER && echo Cancel email sent by second try. || \
            { echo Cancel email still not sent!! Try again.; echo -e "Subject: $s\n$toSend" | sendmail `head -n 1 ~/.forward` && echo Cancel email sent by second try. || echo Cancel email still not sent!!; }

            continue
        fi

        if [ -f $smartSlurmJobRecordDir/stats/$software.$ref.mem.stat ]; then
            output=`estimateMemTime.sh $software $ref $inputSize`
            resAjust="$resAjust\nInputSize: $inputSize\nHere is the mem fomular:\n`cat $smartSlurmJobRecordDir/stats/$software.$ref.mem.stat`\n"
            resAjust="$resAjust\nInputSize: $inputSize\nHere is the time fomular:\n`cat $smartSlurmJobRecordDir/stats/$software.$ref.time.stat`\n"
            resAjust="$resAjust\n#Output from estimateMemTime.sh: $output \n"
            echo "Output from estimateMemTime.sh: $output"

            if [[ "$output" == "outOfRange" ]]; then
                echo Input size is too big for the curve to estimate! Use default mem and runtime to adjust job.
                # not deleting mem.stat, so other jobs will not re-build it within 60 minutes
            elif [ ! -z "$output" ]; then
                output=${output% *}
                #m=${output% *}
                echo extra: .$extraMem.$defaultExtraTime.
                [ "${output% *}" -gt 0 ] && mem=$((${output% *} + extraMem)) && resAjust="$resAjust\n#Give ${extraMem}M extra memory. "
                #t=${output#* }
                [ "${output#* }" -gt 0 ] && min=$((${output#* } + defaultExtraTime)) && resAjust="$resAjust\n#Give $defaultExtraTime more minutes."
                resAjust="$resAjust\n#So use this to adjust the job: $mem M ${min} mins"
            fi
        fi

        if [[ "$output" == "outOfRange" ]] && test `find $smartSlurmJobRecordDir/stats/$software.$ref.mem.stat -mmin +20` || [ ! -f $smartSlurmJobRecordDir/stats/$software.$ref.mem.stat ]; then
            echo "Do not have a formula, or it is old and out of range. Let us build one..."

            echo "Do not have a formula, or it is old and out of range. Let us build one..."  >> $smartSlurmLogDir/$name.out

            #[ test `find $smartSlurmJobRecordDir/jobRecord.txt -mmin -20` ] && echo jobRecord.txt synced within 20 hour. No need to re-sync || cat $HOME/smartSlurm/myJobRecord.txt > $smartSlurmJobRecordDir/jobRecord.txt

            #jobStatistics.sh $software ${ref//\//-} 4 1>&2

            #OUT="$(mktemp -d)"

            #filter by software and reference

            grep COMPLETED $smartSlurmJobRecordDir/jobRecord.txt | awk -F"," -v a=$software -v b=$ref '{ if($12 == a && $13 == b) {print $2, $7 }}' | sort -r  -k1,1 -k2,2 | sort -u -k1,1 > $smartSlurmJobRecordDir/stats/$software.$ref.mem.txt
            grep COMPLETED $smartSlurmJobRecordDir/jobRecord.txt | awk -F"," -v a=$software -v b=$ref '{ if($12 == a && $13 == b) {print $2, $8 }}' | sort -r  -k1,1 -k2,2 | sort -u -k1,1 > $smartSlurmJobRecordDir/stats/$software.$ref.time.txt

            echo "Got mem data from jobRecord.txt (content of $smartSlurmJobRecordDir/stats/$software.$ref.mem.txt):"
            cat $smartSlurmJobRecordDir/stats/$software.$ref.mem.txt

            echo "Got time data from jobRecord.txt (content of $smartSlurmJobRecordDir/stats/$software.$ref.time.txt):"
            cat $smartSlurmJobRecordDir/stats/$software.$ref.time.txt

            if [[ $(wc -l <$smartSlurmJobRecordDir/stats/$software.$ref.mem.txt) -lt 3 ]]; then
                echo There are less than 3 records. No way to fit a curve.
                echo There are less than 3 records. No way to fit a curve. >> $smartSlurmLogDir/$name.out

            else
                gnuplot -e 'set key outside; set key reverse; set key invert; set term png; set output "'"$smartSlurmJobRecordDir/stats/$software.$ref.mem.png"'"; set title "Input Size vs. Memory Usage"; set xlabel "Input Size(K)"; set ylabel "Memory Usage(M)"; f(x)=a*x+b; fit f(x) "'"$smartSlurmJobRecordDir/stats/$software.$ref.mem.txt"'" u 1:2 via a, b; t(a,b)=sprintf("f(x) = %.2fx + %.2f", a, b); plot "'"$smartSlurmJobRecordDir/stats/$software.$ref.mem.txt"'" u 1:2,f(x) t t(a,b); print "Finala=", a; print "Finalb=",b; stats "'"$smartSlurmJobRecordDir/stats/$software.$ref.mem.txt"'" u 1 ' 2>&1 | grep 'Final\| M' | awk 'NF<5{print $1, $2}' |sed 's/:/=/' | sed 's/ //g' > $smartSlurmJobRecordDir/stats/$software.$ref.mem.stat ; echo STDFIT=`cat fit.log | grep FIT_STDFIT | tail -n 1 | awk '{print $8}'` >> $smartSlurmJobRecordDir/stats/$software.$ref.mem.stat

                echo RSquare="$(gnuplot -e 'stats "'"$smartSlurmJobRecordDir/stats/$software.$ref.mem.txt"'" using 1:2;' 2>&1| grep Correlation | cut -d' ' -f7 | awk '{print $1 * $1 }')" >> $smartSlurmJobRecordDir/stats/$software.$ref.mem.stat

                sed -i 's/\x0//g' $smartSlurmJobRecordDir/stats/$software.$ref.mem.stat

                gnuplot -e 'set key outside; set key reverse; set key invert; set term png; set output "'"$smartSlurmJobRecordDir/stats/$software.$ref.time.png"'"; set title "Input Size vs. Time Usage"; set xlabel "Input Size(K)"; set ylabel "Time Usage(Min)"; f(x)=a*x+b; fit f(x) "'"$smartSlurmJobRecordDir/stats/$software.$ref.time.txt"'" u 1:2 via a, b; t(a,b)=sprintf("f(x) = %.2fx + %.2f", a, b); plot "'"$smartSlurmJobRecordDir/stats/$software.$ref.time.txt"'" u 1:2,f(x) t t(a,b); print "Finala=", a; print "Finalb=",b; stats "'"$smartSlurmJobRecordDir/stats/$software.$ref.time.txt"'" u 1 ' 2>&1 | grep 'Final\| M' | awk 'NF<5{print $1, $2}' |sed 's/:/=/' | sed 's/ //g' > $smartSlurmJobRecordDir/stats/$software.$ref.time.stat ; echo STDFIT=`cat fit.log | grep FIT_STDFIT | tail -n 1 | awk '{print $8}'` >> $smartSlurmJobRecordDir/stats/$software.$ref.time.stat

                echo RSquare="$(gnuplot -e 'stats "'"$smartSlurmJobRecordDir/stats/$software.$ref.time.txt"'" using 1:2;' 2>&1| grep Correlation | cut -d' ' -f7 | awk '{print $1 * $1 }')" >> $smartSlurmJobRecordDir/stats/$software.$ref.time.stat

                sed -i 's/\x0//g' $smartSlurmJobRecordDir/stats/$software.$ref.time.stat

                # make plot and calculate statistics
                # gnuplot -e 'set term pdf; set output "time.pdf"; set title "Input Size vs. Time Usage" font "Helvetica Bold,18"; set xlabel "Input Size(K)"; set ylabel "Time(Min)"; f(x)=a*x+b; fit f(x) "time.txt" u 1:2 via a, b; t(a,b)=sprintf("f(x) = %.2fx + %.2f", a, b); plot "time.txt" u 1:2,f(x) t t(a,b); print "Finala=", a; print "Finalb=",b; stats "time.txt" u 1 ' 2>&1 | grep 'Final\| M' | awk 'NF<4{print $1, $2}' |sed 's/:/=/' | sed 's/ //g' > time.stat.txt; echo STDFIT=`cat fit.log | grep FIT_STDFIT | tail -n 1 | awk '{print $8}'` >> time.stat.txt
                # echo RSquare="$(gnuplot -e 'stats "time.txt" using 1:2;' 2>&1| grep Correlation | cut -d' ' -f7 | awk '{print $1 * $1 }')" >> time.stat.txt

                echo There are more than 3 $software $ref jobs already run for this software, statics is ready for current job:
               
                echo
                echo You can see the plot using commands:
                echo display $smartSlurmJobRecordDir/stats/$software.$ref.mem.png
                echo display $smartSlurmJobRecordDir/stats/$software.$ref.time.png

                # cd -



                # echo got files in $smartSlurmJobRecordDir/stats:
                # ls -lrt $smartSlurmJobRecordDir/stats

                if [ -f $smartSlurmJobRecordDir/stats/$software.$ref.mem.stat ]; then
                    output=`estimateMemTime.sh $software $ref $inputSize`
                    #resAjust="$resAjust`cat $smartSlurmJobRecordDir/stats/$software.$ref.mem.stat`\n"
                    resAjust="$resAjust\nInputSize: $inputSize\nHere is the mem fomular:\n`cat $smartSlurmJobRecordDir/stats/$software.$ref.mem.stat`\n"

                    resAjust="$resAjust\nInputSize: $inputSize\nHere is the time fomular:\n`cat $smartSlurmJobRecordDir/stats/$software.$ref.time.stat`\n"

                    resAjust="$resAjust\n#Output from estimateMemTime.sh: $output \n"
                    echo "Output from estimateMemTime.sh: $output"

                    if [[ "$output" == "outOfRange" ]]; then
                        echo Input size is too big for the curve to estimate! Use default mem and runtime to adjust job.
                        # not deleting mem.stat, so other jobs will not re-build it within 60 minutes
                    elif [ ! -z "$output" ]; then
                        output=${output% *}
                        [[ ${output% *} != 0 ]] && mem=$((${output% *}+extraMem)) && resAjust="$resAjust\n#Give ${extraMem}M extra memory. "
                        [[ ${output#* } != 0 ]] && min=$((${output#* }+defaultExtraTime)) && resAjust="$resAjust\n#Give $defaultExtraTime more minutes."
                        resAjust="$resAjust\n#So use this to adjust the job: $mem M ${min} mins"
                    fi
                fi
            fi
            #rm -r $OUT 2>/dev/null
        fi

    # test oom
    #mem=1412

        echo -e "$resAjust"

        echo -e "$resAjust\n" >> $smartSlurmLogDir/$name.out

        
        [ -z "$mem" ] && echo Fail to estimate new resource. Directly release job. && scontrol release $id && continue

        #[ "$mem" -lt 20 ] && mem=20 # at least 20M

        #echo Got estimation inputsize: $inputSize mem: $mem  time: $min

        #echo Got estimation inputsize: $inputSize mem: $mem  time: $min  >> $smartSlurmLogDir/$name.out
        hours=$((($min + 59) / 60))

        echo looking partition for hour: $hours

        adjustPartition $hours partition

        seconds=$(($min * 60))

        time=`eval "echo $(date -ud "@$seconds" +'$((%s/3600/24))-%H:%M:%S')"`

        #set -x 

        #scontrol show job $id

        echo running: scontrol update jobid=$id timelimit=$time partition=$partition MinMemoryNode=${mem}

        echo running: scontrol update jobid=$id timelimit=$time partition=$partition MinMemoryNode=${mem} >> $smartSlurmLogDir/$name.out

        scontrol update JobId=$id TimeLimit=$time Partition=$partition  MinMemoryNode=${mem}
        #scontrol show job $id

        scontrol release $id

        #scontrol show job $id

        echo $mem $min $extraMem > $smartSlurmLogDir/$name.adjust

        #echo -e "Adjusted mem: $mem time: $min (including exralMem: $extraMem)\n" >> $smartSlurmLogDir/$name.out

        #echo $mem $min> $smartSlurmLogDir/$name.adjust
        #touch $smartSlurmLogDir/$name.adjusted
        #echo "scontrol update JobId=$id TimeLimit=$time Partition=$partition  MinMemoryNode=${mem}" >> $smartSlurmLogDir/$name.sh

        #set +x 
    else
        echo Need wait for other jobs to finish before we can ajust mem and runtime...
    fi
done

#rm -r $smartSlurmLogDir/downsteamjob.adjusting 2>/dev/null
