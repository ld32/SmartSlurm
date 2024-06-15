#!/bin/sh

echoerr() { echo "$@" 1>&2; }

usage() { echoerr -e "Usage: \n$0 bowtie2 hg19 inputs\nReturn mem in M and time in minutes."; exit 1; }

#set -x  

echoerr Running: $0 $@

program=$1
ref=$2; ref=${ref//\//-}
inputs=$3
flag=$4 
defaultMem=$5
defaultMin=$6
adjust=$7

mem=""; min=""
[ -z "$defaultMin" ] && usage && exit 1

[ -f $smartSlurmJobRecordDir/stats/extraMem.$program.$ref ] && maxExtra=`sort -n $smartSlurmJobRecordDir/stats/extraMem.$program.$ref | tail -n1 | cut -d' ' -f1` && extraMem=$(( $maxExtra * 2 )) || extraMem=$(( $defaultExtraMem * 2 ))

[ -z "$adjust" ] && resAjust="#Original mem $defaultMem M, Original time: $defaultMin mins\n"

if [ $inputs == none ]; then 
    inputSize=0
    resAjust="$resAjust#This job does not have input.\n"

    rows=`( wc -l $smartSlurmJobRecordDir/stats/$program.$ref.memTime.noInput 2>/dev/null || echo 0 ) | awk '{print $1}'`
    #echoerr rows  $rows
    # empty or more than 60 minutes but less than 4 records
    if test `find $smartSlurmJobRecordDir/stats/$program.$ref.memTime.noInput -mmin +2 2>/dev/null` && [ $rows -lt 200 ] || [ $rows -eq 0 ]; then
        mkdir -p $smartSlurmJobRecordDir/stats/

        grep COMPLETED $smartSlurmJobRecordDir/jobRecord.txt 2>/dev/null | awk -F, -v a=$program -v b=$ref '{ if($12 == a && $13 == b) {print $7, $8 }}' | uniq > $smartSlurmJobRecordDir/stats/$program.$ref.memTime.noInput

        # make plot and calculate statistics
        gnuplot -e 'set key outside; set key reverse; set key invert; set term png; set output "'"$smartSlurmJobRecordDir/stats/$program.$ref.stat.noInput.png"'"; set title "Time vs. Memory Usage"; set xlabel "Time(Min)"; set ylabel "Memory(M)"; f(x)=a*x+b; fit f(x) "'"$smartSlurmJobRecordDir/stats/$program.$ref.timeMem.noInput"'" u 1:2 via a, b; t(a,b)=sprintf("f(x) = %.2fx + %.2f", a, b); plot "'"$smartSlurmJobRecordDir/stats/$program.$ref.timeMem.noInput"'" u 1:2,f(x) t t(a,b); print "Finala=", a; print "Finalb=",b; stats "'"$smartSlurmJobRecordDir/stats/$program.$ref.timeMem.noInput"'" u 1 ' 2>&1 | grep 'Final\| M' | awk 'NF<5{print $1, $2}' | sed 's/:/=/' | sed 's/ //g' > $smartSlurmJobRecordDir/stats/$program.$ref.timeMem.stat.noInput

        rows=`{ wc -l $smartSlurmJobRecordDir/stats/$program.$ref.memTime.noInput 2>/dev/null || echo 0; } | cut -f 1 -d " "`
    fi

    # at least 3 records
    if [ $rows -ge 3 ]; then

        cutoffRow=$(( ($row - 1)  / 10  + 1)) # 90th percentile

        mem=`cut -d' ' -f1 $smartSlurmJobRecordDir/stats/$program.$ref.memTime.noInput | sort -nr | tr '\n' ' ' | cut -f $cutoffRow -d " "`

        mem=$((${mem/\.*/} + extraMem))

        min=`cut -d' ' -f2 $smartSlurmJobRecordDir/stats/$program.$ref.memTime.noInput | sort -nr | tr '\n' ' ' | cut -f $cutoffRow -d " "`

        min=$((${min/\.*/} + defaultExtraTime))

        resAjust="$resAjust#Got estimation based on program.reference: $program.$ref.\n"
        resAjust="$resAjust#Give ${extraMem} M extra memory and $defaultExtraTime more minutes. \n#So use this to submit the job: $mem M  ${min} min"

    else
        resAjust="$resAjust#There are less than 3 job records. Use default mem and time."
    fi


else 
    inputSize=`{ du --apparent-size -c -L ${inputs//,/ } 2>/dev/null || echo notExist; } | tail -n 1 | cut -f 1`
 
    if [[ "$inputSize" == "notExist" ]]; then
        inputSize=missingInputFile
        resAjust="$resAjust#Some or all input files not exist: $inputs\n"
        echoerr Error! missingInputFile: ${inputs//,/ }
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

            #resAjust="$resAjust\n`cat $smartSlurmJobRecordDir/stats/$program.${ref//\//-}.mem.stat`\n"
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

        if [[ "$output" == "outOfRange" ]] && test `find $smartSlurmJobRecordDir/stats/$program.${ref//\//-}.mem.stat -mmin +2` || [ ! -f $smartSlurmJobRecordDir/stats/$program.${ref//\//-}.mem.stat ]; then
            resAjust="$resAjust#Do not have a formula, or it is old and out of range. Let us build one...\n"

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
                gnuplot -e 'set key outside; set key reverse; set key invert; set term png; set output "'"$smartSlurmJobRecordDir/stats/$program.$ref.mem.png"'"; set title "Input Size vs. Memory Usage"; set xlabel "Input Size(K)"; set ylabel "Memory Usage(M)"; f(x)=a*x+b; fit f(x) "'"$smartSlurmJobRecordDir/stats/$program.$ref.mem"'" u 1:2 via a, b; t(a,b)=sprintf("f(x) = %.2fx + %.2f", a, b); plot "'"$smartSlurmJobRecordDir/stats/$program.$ref.mem"'" u 1:2,f(x) t t(a,b); print "Finala=", a; print "Finalb=",b; stats "'"$smartSlurmJobRecordDir/stats/$program.$ref.mem"'" u 1 ' 2>&1 | grep 'Final\| M' | awk 'NF<5{print $1, $2}' |sed 's/:/=/' | sed 's/ //g' > $smartSlurmJobRecordDir/stats/$program.$ref.mem.stat ; echo STDFIT=`cat fit.log | grep FIT_STDFIT | tail -n 1 | awk '{print $8}'` >> $smartSlurmJobRecordDir/stats/$program.$ref.mem.stat

                echo RSquare="$(gnuplot -e 'stats "'"$smartSlurmJobRecordDir/stats/$program.$ref.mem"'" using 1:2;' 2>&1| grep Correlation | cut -d' ' -f7 | awk '{print $1 * $1 }')" >> $smartSlurmJobRecordDir/stats/$program.$ref.mem.stat
                
                echo SCount=$(wc -l $smartSlurmJobRecordDir/stats/$program.$ref.mem | cut -d' ' -f1) >> $smartSlurmJobRecordDir/stats/$program.$ref.mem.stat

                sed -i 's/\x0//g' $smartSlurmJobRecordDir/stats/$program.$ref.mem.stat

                gnuplot -e 'set key outside; set key reverse; set key invert; set term png; set output "'"$smartSlurmJobRecordDir/stats/$program.$ref.time.png"'"; set title "Input Size vs. Time Usage"; set xlabel "Input Size(K)"; set ylabel "Time Usage(Min)"; f(x)=a*x+b; fit f(x) "'"$smartSlurmJobRecordDir/stats/$program.$ref.time"'" u 1:2 via a, b; t(a,b)=sprintf("f(x) = %.2fx + %.2f", a, b); plot "'"$smartSlurmJobRecordDir/stats/$program.$ref.time"'" u 1:2,f(x) t t(a,b); print "Finala=", a; print "Finalb=",b; stats "'"$smartSlurmJobRecordDir/stats/$program.$ref.time"'" u 1 ' 2>&1 | grep 'Final\| M' | awk 'NF<5{print $1, $2}' |sed 's/:/=/' | sed 's/ //g' > $smartSlurmJobRecordDir/stats/$program.$ref.time.stat ; echo STDFIT=`cat fit.log | grep FIT_STDFIT | tail -n 1 | awk '{print $8}'` >> $smartSlurmJobRecordDir/stats/$program.$ref.time.stat

                echo RSquare="$(gnuplot -e 'stats "'"$smartSlurmJobRecordDir/stats/$program.$ref.time"'" using 1:2;' 2>&1| grep Correlation | cut -d' ' -f7 | awk '{print $1 * $1 }')" >> $smartSlurmJobRecordDir/stats/$program.$ref.time.stat

                sed -i 's/\x0//g' $smartSlurmJobRecordDir/stats/$program.$ref.time.stat

                echoerr
                echoerr You can see the plot using commands:
                echoerr display $smartSlurmJobRecordDir/stats/$program.$ref.mem.png
                echoerr display $smartSlurmJobRecordDir/stats/$program.$ref.time.png

                #echoerr got files in $smartSlurmJobRecordDir/stats:
                #ls -lrt $smartSlurmJobRecordDir/stats
                if [ -f $smartSlurmJobRecordDir/stats/$program.${ref//\//-}.mem.stat ]; then
                    if [ -f .command.sh ] && [ -f .command.run ]; then 
                        output=`estimateMemTime.sh $program ${ref//\//-} $inputSize 2>> ../../../.nextflow.log`
                    else #if [[ "$parentCmd" == */bin/snakemake* ]]; then
                        output=`estimateMemTime.sh $program ${ref//\//-} $inputSize 2>> .smartSlurm.log`
                    fi 

                    #resAjust="$resAjust\n`cat $smartSlurmJobRecordDir/stats/$program.${ref//\//-}.mem.stat`\n"
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
        fi
    fi
fi 

if [ -z $mem ] || [ -z $min ]; then 
    mem=$defaultMem    
    min=$defaultMin
else 
    
    #[ "$mem" -lt 100 ] && mem=100 && resAjust="$resAjust\n#Mem is reset to 100M. "
    #[ "$min" -lt 10 ] && min=10 && resAjust="$resAjust\n#Time is reset to 10min. "

    if [[ $adjust == "adjust" ]]; then 
        echo $inputSize $mem $min $extraMem >> $smartSlurmLogDir/$flag.adjust
    fi  
fi
echo -e "$resAjust" >> $smartSlurmLogDir/$flag.out       
echo $inputSize $mem $min $extraMem

echoerr Got $inputSize $mem $min $extraMem