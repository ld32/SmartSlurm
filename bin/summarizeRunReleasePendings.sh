#!/bin/sh
Usage="Usage: $0 [ workDir/log, the log folder name. ]  \nThis script will go through job name list in allJobs.txt to see if the jobs finish successfully or not."

#set -x

echo Running: $0 $@

#cd $1 #logDir

if [ -f $smartSlurmLogDir/allJobs.txt ]; then 
    lines=`tail -n +2 $smartSlurmLogDir/allJobs.txt` # | awk 'NF>2{print $1, $2, $3}'`
else 
    exit 1; 
    #lines="$SLURMJOB_ID $2" # tail -n +2 allJobs.txt | awk 'NF>2{print $1, $2, $3}'`
fi 
 

echo pwd: `pwd`

IFS=$'\n'

out=`squeue -u $USER -t PD,R --noheader -o "%.18i-%t"`

current=0; succ=0; fail=0; running=0; pending=0; requeue=0; unknown=0; unholdCounter=0; 
toSend="Summery for jobs in allJobs.txt:"
for line in $lines; do
    if [ ! -z "${line/ /}" ]; then
        #id=${line%% *}; name=${line##* }
        IFS=' ' read -a arrIN <<< "$line"

        id=${arrIN[0]}
        deps=${arrIN[1]}
        name=${arrIN[2]}
        program=${arrIN[3]}
        ref=${arrIN[4]}  ; ref=${ref//\//-}
        inputs="${arrIN[5]}"

        if [ -f $smartSlurmLogDir/$name.success ]; then
            toSend="$toSend\n${line:0:40} Done"
            succ=$((succ + 1))
        elif [ -f $smartSlurmLogDir/$name.failed ]; then
            toSend="$toSend\n${line:0:40} Failed"
            fail=$((fail + 1))
        elif [[ "$out" == *$id-R* ]]; then # && [[ "$id" != "$SLURM_JOBID" ]]; then
            toSend="$toSend\n${line:0:40} Running"
            running=$((running + 1))
        elif [[ "$out" == *$id-P* ]]; then # && [[ "$id" != "$SLURM_JOBID" ]]; then
            toSend="$toSend\n${line:0:40} Pending"
            pending=$((pending + 1))
            if [ $unholdCounter -gt 0 ]; then 
                if [[ "$deps" == null ]]; then 
                    unholdCounter=$((unholdCounter - 1))

                    echo Trying to release: $line 

                    # this copied from ssbatch
                    if [[ "$inputs" == "none" ]]; then
                        echo No inputs
                        resAjust="$resAjust#This job does not have input.\n"
                        ref=${ref//\//-}

                        rows=`( wc -l $smartSlurmJobRecordDir/stats/$program.$ref.mem.noInput 2>/dev/null || echo 0 ) | awk '{print $1}'`
                        #echo rows  $rows
                        # empty or more than 60 minutes but less than 4 records
                        if test `find $smartSlurmJobRecordDir/stats/$program.$ref.mem.noInput -mmin +20 2>/dev/null` && [ $rows -lt 200 ] || [ $rows -eq 0 ]; then
                            mkdir -p $smartSlurmJobRecordDir/stats/
                            #cat /home/*/smartSlurm/stats/myJobRecord.txt > $smartSlurmJobRecordDir/jobRecord.txt
                            #test `find $smartSlurmJobRecordDir/jobRecord.txt -mmin +20` && echo jobRecord.txt synced within 21 hour. No need to re-sync || cat $HOME/smartSlurm/myJobRecord.txt > $smartSlurmJobRecordDir/jobRecord.txt

                            # todo: could use single file here
                            grep COMPLETED $smartSlurmJobRecordDir/jobRecord.txt 2>/dev/null | awk -F, -v a=$program -v b=$ref '{ if($12 == a && $13 == b) {print $7 }}' | uniq > $smartSlurmJobRecordDir/stats/$program.$ref.mem.noInput

                            grep COMPLETED $smartSlurmJobRecordDir/jobRecord.txt 2>/dev/null | awk -F, -v a=$program -v b=$ref '{ if($12 == a && $13 == b) {print $8 }}' | uniq > $smartSlurmJobRecordDir/stats/$program.$ref.time.noInput

                            if [ -s $smartSlurmJobRecordDir/stats/$program.$ref.mem.noInput ] && [ -s $smartSlurmJobRecordDir/stats/$program.$ref.time.noInput ]; then

                                #OUT="$(mktemp -d)"
                                paste $smartSlurmJobRecordDir/stats/$program.$ref.time.noInput $smartSlurmJobRecordDir/stats/$program.$ref.mem.noInput | column -s $'\t' -t | sed '$ d' > $smartSlurmJobRecordDir/stats/$program.$ref.timeMem.noInput

                              
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
                    else 

                    #inputs1=(numbers3.txt numbers1.txt)

# Store the output of the command in a variable
#output=$(du --apparent-size -c -L "${inputs1[@]}")
                        IFS=','
                        read -ra input_array <<< "${inputs#,}"
                        inputSize=`{ du --apparent-size -c -L ${input_array[@]} 2>/dev/null || echo notExist; } | tail -n 1 | cut -f 1`

                        if [[ "$inputSize" == "notExist" ]]; then
                            resAjust="$resAjust#Some or all input files not exist: $inputs\n"
                            echo Error! missingInputFile: ${inputs//,/ }
                            
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

                                #[ test `find $smartSlurmJobRecordDir/jobRecord.txt -mmin -20` ] && echo jobRecord.txt synced within 20 hour. No need to re-sync || cat $HOME/smartSlurm/myJobRecord.txt > $smartSlurmJobRecordDir/jobRecord.txt

                                #jobStatistics.sh $program ${ref//\//-} 4 1>&2

                                #filter by program and reference
                                # todo: maybe able to replace / in ref at begaining of the script?
                                ref=${ref//\//-}
                                grep COMPLETED $smartSlurmJobRecordDir/jobRecord.txt 2>/dev/null | awk -F"," -v a=$program -v b=$ref '{ if($12 == a && $13 == b && $2!=0) {print $2, $7 }}' | sort -r  -k1,1 -k2,2 | sort -u -k1,1 > $smartSlurmJobRecordDir/stats/$program.$ref.mem

                                grep COMPLETED $smartSlurmJobRecordDir/jobRecord.txt 2>/dev/null | awk -F"," -v a=$program -v b=$ref '{ if($12 == a && $13 == b && $2!=0) {print $2, $8 }}' | sort -r  -k1,1 -k2,2 | sort -u -k1,1 > $smartSlurmJobRecordDir/stats/$program.$ref.time

                                echo "Got mem data from jobRecord.txt (content of $smartSlurmJobRecordDir/stats/$program.$ref.mem):"
                                echo `cat $smartSlurmJobRecordDir/stats/$program.$ref.mem`

                                echo "Got time data from jobRecord.txt (content of $smartSlurmJobRecordDir/stats/$program.$ref.time):"
                                echo `cat $smartSlurmJobRecordDir/stats/$program.$ref.time`

                                if [[ $(wc -l <$smartSlurmJobRecordDir/stats/$program.$ref.mem) -lt 3 ]]; then
                                    echo There are less than 3 records. No way to fit a curve. User defaut values...
                                    resAjust="$resAjust#There are less than 3 records. No way to fit a curve.\n"

                                else



                           

                                    gnuplot -e 'set key outside; set key reverse; set key invert; set term png; set output "'"$smartSlurmJobRecordDir/stats/$program.$ref.mem.png"'"; set title "Input Size vs. Memory Usage"; set xlabel "Input Size(K)"; set ylabel "Memory Usage(M)"; f(x)=a*x+b; fit f(x) "'"$smartSlurmJobRecordDir/stats/$program.$ref.mem"'" u 1:2 via a, b; t(a,b)=sprintf("f(x) = %.2fx + %.2f", a, b); plot "'"$smartSlurmJobRecordDir/stats/$program.$ref.mem"'" u 1:2,f(x) t t(a,b); print "Finala=", a; print "Finalb=",b; stats "'"$smartSlurmJobRecordDir/stats/$program.$ref.mem"'" u 1 ' 2>&1 | grep 'Final\| M' | awk 'NF<5{print $1, $2}' |sed 's/:/=/' | sed 's/ //g' > $smartSlurmJobRecordDir/stats/$program.$ref.mem.stat ; echo STDFIT=`cat fit.log | grep FIT_STDFIT | tail -n 1 | awk '{print $8}'` >> $smartSlurmJobRecordDir/stats/$program.$ref.mem.stat

                                    echo RSquare="$(gnuplot -e 'stats "'"$smartSlurmJobRecordDir/stats/$program.$ref.mem"'" using 1:2;' 2>&1| grep Correlation | cut -d' ' -f7 | awk '{print $1 * $1 }')" >> $smartSlurmJobRecordDir/stats/$program.$ref.mem.stat
                                
                                    echo SCount=$(wc -l $smartSlurmJobRecordDir/stats/$program.$ref.mem | cut -d' ' -f1) >> $smartSlurmJobRecordDir/stats/$program.$ref.mem.stat

                                    sed -i 's/\x0//g' $smartSlurmJobRecordDir/stats/$program.$ref.mem.stat

                                    gnuplot -e 'set key outside; set key reverse; set key invert; set term png; set output "'"$smartSlurmJobRecordDir/stats/$program.$ref.time.png"'"; set title "Input Size vs. Time Usage"; set xlabel "Input Size(K)"; set ylabel "Time Usage(Min)"; f(x)=a*x+b; fit f(x) "'"$smartSlurmJobRecordDir/stats/$program.$ref.time"'" u 1:2 via a, b; t(a,b)=sprintf("f(x) = %.2fx + %.2f", a, b); plot "'"$smartSlurmJobRecordDir/stats/$program.$ref.time"'" u 1:2,f(x) t t(a,b); print "Finala=", a; print "Finalb=",b; stats "'"$smartSlurmJobRecordDir/stats/$program.$ref.time"'" u 1 ' 2>&1 | grep 'Final\| M' | awk 'NF<5{print $1, $2}' |sed 's/:/=/' | sed 's/ //g' > $smartSlurmJobRecordDir/stats/$program.$ref.time.stat ; echo STDFIT=`cat fit.log | grep FIT_STDFIT | tail -n 1 | awk '{print $8}'` >> $smartSlurmJobRecordDir/stats/$program.$ref.time.stat

                                    echo RSquare="$(gnuplot -e 'stats "'"$smartSlurmJobRecordDir/stats/$program.$ref.time"'" using 1:2;' 2>&1| grep Correlation | cut -d' ' -f7 | awk '{print $1 * $1 }')" >> $smartSlurmJobRecordDir/stats/$program.$ref.time.stat

                                    sed -i 's/\x0//g' $smartSlurmJobRecordDir/stats/$program.$ref.time.stat
                                   
                                    echo There are more than 3 $program $ref jobs already run for this program, statics is ready for current job:
                                    # echo Memeory statisics:
                                    # echo "inputsize mem(M)"
                                    # cat $smartSlurmJobRecordDir/stats/$program.$ref.mem.stat
                                    # echo
                                    # echo Time statistics:
                                    # echo "inputsize time(minute)"
                                    # cat $smartSlurmJobRecordDir/stats/$program.$ref.time.stat

                                    # mv $OUT/mem.txt $smartSlurmJobRecordDir/stats/$program.$ref.mem
                                    # mv $OUT/time.txt $smartSlurmJobRecordDir/stats/$program.$ref.time

                                    # convert $OUT/mem.pdf -background White -flatten $smartSlurmJobRecordDir/stats/$program.$ref.mem.pdf
                                    # convert $OUT/time.pdf -background White -flatten $smartSlurmJobRecordDir/stats/$program.$ref.time.pdf
                                    # pdftoppm $smartSlurmJobRecordDir/stats/$program.$ref.mem.pdf  -png > $smartSlurmJobRecordDir/stats/$program.$ref.mem.png
                                    # pdftoppm $smartSlurmJobRecordDir/stats/$program.$ref.time.pdf  -png > $smartSlurmJobRecordDir/stats/$program.$ref.time.png

                                    echo
                                    echo You can see the plot using commands:
                                    echo display $smartSlurmJobRecordDir/stats/$program.$ref.mem.png
                                    echo display $smartSlurmJobRecordDir/stats/$program.$ref.time.png

                                    # cd -



                                    #echo got files in $smartSlurmJobRecordDir/stats:
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
                                        #echo got estimation $output
                                    fi

                                fi
                                #rm -r $OUT 2>/dev/null
                            fi

                        fi
                    fi

                    # wait for until get estimation
                    if [ ! -z "$min" ] && [ ! -z "$mem" ]; then
                    #    echo did not get new min and mem. directly release job. 
                    #    scontrol release $id 
                    #else 

                        # this part was from adjust downsteamjobs
                        echo -e "$resAjust"

                        echo -e "$resAjust\n" >> $smartSlurmLogDir/$name.out

                        
                        

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
                    
                    # didn't get estimate, but already have 3 successful jobs, release one job anyway
                    # because jobRecord need to be unique by programName + reference + inputSize + memory
                    # If the first 5 jobs have same unique value, there is only one record in jobRecords.txt
                    elif [[ "$succ" -gt 2 ]]; then 
                        unholdCounter=0; 
                        echo Fail to estimate new resource. But alreay have 3 success jobs, directly release one job anyway. 
                        scontrol release $id 
                        continue
                    fi 
                fi
            fi            
            

        elif [ -f $smartSlurmLogDir/$name.failed.requeued.1.time ]; then 
            toSend="$toSend\n${line:0:40} Requeued"
            requeue=$((requeue + 1))    
        else
            toSend="$toSend\n${line:0:40} Unknow"
            unknown=$((unknown + 1))
        fi
        if [ "$id" == "$SLURM_JOBID" ]; then 
            #check if statics is available for new ten hoding jobs
            unholdCounter=5 # only take care of jobs without dependency
        fi 
    fi
done

current=$((succ + fail + requeue))
total=$((succ + fail + running + pending + +requeue + unknown))
s="$current/$total Succ:$succ/$total Requeue:$requeue/$total Running:$running/$total Pending:$pending/$total Fail:$fail/$total Unknown:$unknown/$total"

[ -f $smartSlurmLogDir/allJobs.txt ] && echo -e "$s\n$toSend" > $smartSlurmLogDir/summary.$SLURMJOB_ID


# if [ $((running + pending)) -le 5 ]; then
#     echo -e "$toSend" | mail -s $s $USER
#     [ "$USER" != ld32 ] && echo -e "$toSend" | mail -s $s ld32
# fi

[ ! -f $smartSlurmLogDir/job_$SLURM_JOBID.memCPU.txt ] && echo Not found $smartSlurmLogDir/job_$SLURM_JOBID.memCPU.txt && exit 


# todo: should make the plot wider instead of shink it: 
# https://stackoverflow.com/questions/13869439/gnuplot-how-to-increase-the-width-of-my-graph

rowTotal=`wc -l $smartSlurmLogDir/job_$SLURM_JOBID.memCPU.txt | cut -d' ' -f1`
if [ "$rowTotal" -gt 50 ]; then 
    maxMem=0; maxCpu=0; 
    rate=`echo "scale=2;$rowTotal/50"|bc`
    IFS=$'\n'; rowCount1=0; rowCount2=0
    echo > $smartSlurmLogDir/job_$SLURM_JOBID.memCPU1.txt
    for t in `cat $smartSlurmLogDir/job_$SLURM_JOBID.memCPU.txt`; do
        mem=`echo $t | cut -d' ' -f2`
        cpu=`echo $t | cut -d' ' -f5`
        [ "$mem" -gt $maxMem ] && maxMem=$mem && mem1=`echo $t | cut -d' ' -f3,4`
        [ "$cpu" -gt $maxCpu ] && maxCpu=$cpu
        rowCount1=$((rowCount1 + 1))
        rowMax=`echo "scale=2;$rowCount2*$rate"|bc`
        rowMax=${rowMax%.*}; [ -z "$rowMax" ] && rowMax=1; 
        if [ "$rowMax" -le "$rowCount1" ]; then 
            rowCount2=$((rowCount2 + 1))
            echo $rowCount2 $maxMem $mem1 $maxCpu >> $smartSlurmLogDir/job_$SLURM_JOBID.memCPU1.txt
            maxMem=0; maxCpu=0;
        fi 
    done
else 
    cp $smartSlurmLogDir/job_$SLURM_JOBID.memCPU.txt $smartSlurmLogDir/job_$SLURM_JOBID.memCPU1.txt
fi 

# time vs. memory for current job
gnuplot -e "set key outside; set key reverse; set key invert; set datafile separator ' '; set style data histogram; set style histogram rowstacked gap 2; set style fill solid border rgb 'black'; set xtics rotate by -45; set terminal png size 800,600; set output '$smartSlurmLogDir/job_$SLURM_JOBID.mem.png'; set title 'Time vs. Mem for job $SLURM_JOBID'; set xlabel 'Time'; set ylabel 'Mem (M)'; plot '$smartSlurmLogDir/job_$SLURM_JOBID.memCPU1.txt' using 2:xtic(1) title 'Used' lc rgb 'green', '' using 3:xtic(1) title 'Wasted' lc rgb 'red', '' using 4:xtic(1) title 'Saved' lc rgb 'yellow'"

# time vs. CPU usage for current job
gnuplot -e "set key outside; set key reverse; set key invert; set datafile separator ' '; set style data histogram; set style histogram rowstacked gap 2; set style fill solid border rgb 'black'; set xtics rotate by -45; set terminal png size 800,600; set output '$smartSlurmLogDir/job_$SLURM_JOBID.cpu.png'; set title 'Time vs. CPU Usage for job $SLURM_JOBID'; set xlabel 'Time'; set ylabel 'CPU Usage (%)'; plot '$smartSlurmLogDir/job_$SLURM_JOBID.memCPU1.txt' using 5:xtic(1) title 'Used' lc rgb 'green'"



