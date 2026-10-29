#!/bin/sh

echoerr() { echo "$@" 1>&2; }

usage() { echoerr -e "Usage: \n$0 bowtie2 hg19 inputs\nReturn mem in M and time in minutes."; exit 1; }

#set -x

# v2 jobRecords: schema SSOT (jrEmitV2 compat shim reads v1 or v2; JR_AWK_COLMAP
# gives name->index). Readers below parse v2 by name via the shim.
. "$(dirname "$0")/jobRecordSchema.sh" 2>/dev/null

output="Running: $0"
for param in "$@"; do
    if [[ "$param" == *\ * ]]; then
        output="$output \"$param\""
    else
        output="$output $param"
    fi
done
echoerr "$output"

program=$1
ref=$2; ref=${ref//\//-}
inputs=$3
flag=$4 
defaultMem=$5
defaultMin=$6
adjust=$7
pipeline=${8:-unknown}

# v2 grouping: fit per (pipeline, step, ref). refGroup allows optional size-class
# pooling of similar-sized refs (P4c ref-class); for now it is the ref itself.
# statKey names the per-group stats files. refRe matches the refs in this group.
refGroup=$(jrRefClass "$ref" 2>/dev/null); [ -n "$refGroup" ] || refGroup="$ref"
refRe=$(jrRefMatchRe "$ref" 2>/dev/null); [ -n "$refRe" ] || refRe="^$ref\$"
statKey="$pipeline.$program.$refGroup"

# jrExtractBoth <memOutfile> <timeOutfile> : build BOTH the inputSize->mem and
# inputSize->time extracts for this (pipeline, step, refGroup) in a SINGLE pass over
# jobRecord.txt (one jrEmitV2 + one awk emitting both). Fallback ladder: if the
# pipeline-scoped history has fewer than 3 points (too thin to fit -- e.g. right
# after migration, where old rows are pipeline=legacy), re-emit pooling across
# pipelines for the same step + ref group. Then dedupe each by input size.
jrExtractBoth() {
    local memOut=$1 timeOut=$2 jr="$smartSlurmJobRecordDir/jobRecord.txt"
    local rawM rawT; rawM=$(mktemp); rawT=$(mktemp)
    # $1 of the inner filter is the pipeline to match ("" = pool across pipelines).
    _jrEmitPair() {
        jrEmitV2 "$jr" 2>/dev/null | awk -F, -v p="$1" -v a="$program" -v rre="$refRe" \
            -v mo="$rawM" -v to="$rawT" "$JR_AWK_COLMAP"'
            $(col["status"])=="COMPLETED" && (p=="" || $(col["pipeline"])==p) && \
            $(col["step"])==a && $(col["ref"]) ~ rre && $(col["inputSize"]) ~ /^[0-9]+$/ && $(col["inputSize"])+0>0 {
                print $(col["inputSize"]), $(col["memUsed"])  > mo
                print $(col["inputSize"]), $(col["timeUsed"]) > to
            }'
    }
    : > "$rawM"; : > "$rawT"
    _jrEmitPair "$pipeline"
    if [ "$(wc -l < "$rawM" 2>/dev/null || echo 0)" -lt 3 ]; then
        : > "$rawM"; : > "$rawT"
        _jrEmitPair ""
    fi
    sort -r -k1,1 -k2,2 "$rawM" | sort -u -k1,1 > "$memOut"
    sort -r -k1,1 -k2,2 "$rawT" | sort -u -k1,1 > "$timeOut"
    rm -f "$rawM" "$rawT"
}

mem=""; min=""
[ -z "$defaultMin" ] && usage && exit 1

[ -f $smartSlurmJobRecordDir/stats/extraMem.$program.$ref ] && maxExtra=`sort -n $smartSlurmJobRecordDir/stats/extraMem.$program.$ref | tail -n1 | cut -d' ' -f1` && extraMem=$(( $maxExtra * 2 )) || extraMem=$defaultExtraMem

[ -z "$adjust" ] && resAjust="#Original mem $defaultMem M, Original time: $defaultMin mins\n"

if [ $inputs == none ]; then 
    inputSize=0
    resAjust="$resAjust#This job does not have input.\n"

    rows=`( wc -l $smartSlurmJobRecordDir/stats/$statKey.memTime.noInput 2>/dev/null || echo 0 ) | awk '{print $1}'`
    #echoerr rows  $rows
    # empty or more than 60 minutes but less than 4 records
    if test `find $smartSlurmJobRecordDir/stats/$statKey.memTime.noInput -mmin +2 2>/dev/null` && [ $rows -lt 200 ] || [ $rows -eq 0 ]; then
        mkdir -p $smartSlurmJobRecordDir/stats/

        # No-input history (percentile source). Mirror jrExtractBoth's fallback
        # ladder: query this pipeline first, and if it has fewer than 3 records,
        # re-query pooling across ALL pipelines (same step + ref group) so a job
        # without -I still gets a percentile from past usage instead of the bare
        # default (e.g. right after migration, where old rows are pipeline=legacy).
        noInAwk='$(col["status"])=="COMPLETED" && (p=="" || $(col["pipeline"])==p) && $(col["step"])==a && $(col["ref"]) ~ rre { print $(col["memUsed"]), $(col["timeUsed"]) }'
        nif="$smartSlurmJobRecordDir/stats/$statKey.memTime.noInput"
        jrEmitV2 "$smartSlurmJobRecordDir/jobRecord.txt" 2>/dev/null | awk -F, -v p="$pipeline" -v a="$program" -v rre="$refRe" "$JR_AWK_COLMAP""$noInAwk" | uniq > "$nif"
        if [ "$( { wc -l < "$nif" 2>/dev/null || echo 0; } )" -lt 3 ]; then
            jrEmitV2 "$smartSlurmJobRecordDir/jobRecord.txt" 2>/dev/null | awk -F, -v p="" -v a="$program" -v rre="$refRe" "$JR_AWK_COLMAP""$noInAwk" | uniq > "$nif"
        fi

        # make plot and calculate statistics
        gnuplot -e 'set key outside; set key reverse; set key invert; set fit logfile "/dev/null"; set fit quiet; set term png; set output "'"$smartSlurmJobRecordDir/stats/$statKey.stat.noInput.png"'"; set title "Time vs. Memory Usage"; set xlabel "Time(Min)"; set ylabel "Memory(M)"; f(x)=a*x+b; fit f(x) "'"$smartSlurmJobRecordDir/stats/$statKey.memTime.noInput"'" u 1:2 via a, b; t(a,b)=sprintf("f(x) = %.2fx + %.2f", a, b); plot "'"$smartSlurmJobRecordDir/stats/$statKey.memTime.noInput"'" u 1:2,f(x) t t(a,b); print "Finala=", a; print "Finalb=",b; print "STDFIT=", FIT_STDFIT; stats "'"$smartSlurmJobRecordDir/stats/$statKey.memTime.noInput"'" u 1 ' 2>&1 | grep -E 'Final|STDFIT|Mean:|Minimum:|Maximum:|Median:' | awk 'NF<5{print $1, $2}' | sed 's/:/=/' | sed 's/ //g' > $smartSlurmJobRecordDir/stats/$statKey.memTime.stat.noInput

        rows=`{ wc -l $smartSlurmJobRecordDir/stats/$statKey.memTime.noInput 2>/dev/null || echo 0; } | cut -f 1 -d " "`
    fi

    # at least 3 records
    if [ $rows -ge 3 ]; then

        cutoffRow=$(( ($rows - 1)  / 10  + 1)) # 90th percentile

        mem=`cut -d' ' -f1 $smartSlurmJobRecordDir/stats/$statKey.memTime.noInput | sort -nr | tr '\n' ' ' | cut -f $cutoffRow -d " "`

        mem=$((${mem/\.*/} + extraMem))

        min=`cut -d' ' -f2 $smartSlurmJobRecordDir/stats/$statKey.memTime.noInput | sort -nr | tr '\n' ' ' | cut -f $cutoffRow -d " "`

        min=$((${min/\.*/} + defaultExtraTime))

        resAjust="$resAjust#Got estimation based on program.reference: $program.$ref.\n"
        resAjust="$resAjust#Give ${extraMem} M extra memory and $defaultExtraTime more minutes. \n#So use this to submit the job: $mem M  ${min} min"

    else
        resAjust="$resAjust#There are less than 3 job records. Use default mem and time."
    fi
else 
    if [[ "$inputs" =~ ,jobSize:[0-9]+$ ]]; then 
        inputSize=${inputs#,jobSize:}
    else 
        inputSize=`{ du --apparent-size -c -L ${inputs//,/ } 2>/dev/null || echo notExist; } | tail -n 1 | cut -f 1`
    fi 
    if [[ "$inputSize" == "notExist" ]]; then
        inputSize=missingInputFile
        resAjust="$resAjust#Some or all input files not exist: $inputs\n"
        echoerr Error! missingInputFile: ${inputs//,/ }
    else
        #inputSize=$(($inputSize/1024)); # convert to M
        resAjust="$resAjust#InputSize: $inputSize\n"

        #rm ~/.rcbio/$program.$ref.mem.stat # for testing
        if [ -f $smartSlurmJobRecordDir/stats/$statKey.mem.stat ]; then
            
            if [ -f .command.sh ] && [ -f .command.run ]; then 
                output=`calculateMemTime.sh $program ${ref//\//-} $inputSize $pipeline 2>> ../../../.nextflow.log`
            else #if [[ "$parentCmd" == */bin/snakemake* ]]; then
                output=`calculateMemTime.sh $program ${ref//\//-} $inputSize $pipeline`
            fi 

            #resAjust="$resAjust\n`cat $smartSlurmJobRecordDir/stats/$statKey.mem.stat`\n"
            resAjust="$resAjust\n#Output from calculateMemTime.sh: $output \n"
            
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

        if [[ "$output" == "outOfRange" ]] && test `find $smartSlurmJobRecordDir/stats/$statKey.mem.stat -mmin +2` || [ ! -f $smartSlurmJobRecordDir/stats/$statKey.mem.stat ]; then
            resAjust="$resAjust#Do not have a formula, or it is old and out of range. Let us build one...\n"

            #filter by program and reference
            # todo: maybe able to replace / in ref at begaining of the script?
            ref=${ref//\//-}
            # single pass builds both extracts (see jrExtractBoth)
            jrExtractBoth "$smartSlurmJobRecordDir/stats/$statKey.mem" "$smartSlurmJobRecordDir/stats/$statKey.time"

            echoerr "Got mem data from jobRecord.txt (content of $smartSlurmJobRecordDir/stats/$statKey.mem):"
            echoerr `cat $smartSlurmJobRecordDir/stats/$statKey.mem`

            echoerr "Got time data from jobRecord.txt (content of $smartSlurmJobRecordDir/stats/$statKey.time):"
            echoerr `cat $smartSlurmJobRecordDir/stats/$statKey.time`

            if [[ $(wc -l <$smartSlurmJobRecordDir/stats/$statKey.mem) -lt 3 ]]; then
                echoerr There are less than 3 records. No way to fit a curve. User defaut values...
                resAjust="$resAjust#There are less than 3 records. No way to fit a curve.\n"

            else
                gnuplot -e 'set key outside; set key reverse; set key invert; set fit logfile "/dev/null"; set fit quiet; set term png; set output "'"$smartSlurmJobRecordDir/stats/$statKey.mem.png"'"; set title "Input Size vs. Memory Usage"; set xlabel "Input Size(K)"; set ylabel "Memory Usage(M)"; f(x)=a*x+b; fit f(x) "'"$smartSlurmJobRecordDir/stats/$statKey.mem"'" u 1:2 via a, b; t(a,b)=sprintf("f(x) = %.2fx + %.2f", a, b); plot "'"$smartSlurmJobRecordDir/stats/$statKey.mem"'" u 1:2,f(x) t t(a,b); print "Finala=", a; print "Finalb=",b; print "STDFIT=", FIT_STDFIT; stats "'"$smartSlurmJobRecordDir/stats/$statKey.mem"'" u 1 ' 2>&1 | grep -E 'Final|STDFIT|Mean:|Minimum:|Maximum:|Median:' | awk 'NF<5{print $1, $2}' |sed 's/:/=/' | sed 's/ //g' > $smartSlurmJobRecordDir/stats/$statKey.mem.stat

                #echo RSquare="$(gnuplot -e 'stats "'"$smartSlurmJobRecordDir/stats/$statKey.mem"'" using 1:2;' 2>&1| grep Correlation | cut -d' ' -f7 | awk '{print $1 * $1 }')" >> $smartSlurmJobRecordDir/stats/$statKey.mem.stat
                echo RSquare="$(gnuplot -e 'stats "'"$smartSlurmJobRecordDir/stats/$statKey.mem"'" using 1:2;' 2>&1| grep -i Correlation | awk '{print $NF * $NF}')" >> $smartSlurmJobRecordDir/stats/$statKey.mem.stat
                 
                echo SCount=$(wc -l $smartSlurmJobRecordDir/stats/$statKey.mem | cut -d' ' -f1) >> $smartSlurmJobRecordDir/stats/$statKey.mem.stat

                sed -i 's/\x0//g' $smartSlurmJobRecordDir/stats/$statKey.mem.stat

                gnuplot -e 'set key outside; set key reverse; set key invert; set fit logfile "/dev/null"; set fit quiet; set term png; set output "'"$smartSlurmJobRecordDir/stats/$statKey.time.png"'"; set title "Input Size vs. Time Usage"; set xlabel "Input Size(K)"; set ylabel "Time Usage(Min)"; f(x)=a*x+b; fit f(x) "'"$smartSlurmJobRecordDir/stats/$statKey.time"'" u 1:2 via a, b; t(a,b)=sprintf("f(x) = %.2fx + %.2f", a, b); plot "'"$smartSlurmJobRecordDir/stats/$statKey.time"'" u 1:2,f(x) t t(a,b); print "Finala=", a; print "Finalb=",b; print "STDFIT=", FIT_STDFIT; stats "'"$smartSlurmJobRecordDir/stats/$statKey.time"'" u 1 ' 2>&1 | grep -E 'Final|STDFIT|Mean:|Minimum:|Maximum:|Median:' | awk 'NF<5{print $1, $2}' |sed 's/:/=/' | sed 's/ //g' > $smartSlurmJobRecordDir/stats/$statKey.time.stat

                #echo RSquare="$(gnuplot -e 'stats "'"$smartSlurmJobRecordDir/stats/$statKey.time"'" using 1:2;' 2>&1| grep Correlation | cut -d' ' -f7 | awk '{print $1 * $1 }')" >> $smartSlurmJobRecordDir/stats/$statKey.time.stat
                echo RSquare="$(gnuplot -e 'stats "'"$smartSlurmJobRecordDir/stats/$statKey.time"'" using 1:2;' 2>&1| grep -i Correlation | awk '{print $NF * $NF}')" >> $smartSlurmJobRecordDir/stats/$statKey.time.stat
 
                sed -i 's/\x0//g' $smartSlurmJobRecordDir/stats/$statKey.time.stat

                echoerr
                echoerr You can see the plot using commands:
                echoerr display $smartSlurmJobRecordDir/stats/$statKey.mem.png
                echoerr display $smartSlurmJobRecordDir/stats/$statKey.time.png

                #echoerr got files in $smartSlurmJobRecordDir/stats:
                #ls -lrt $smartSlurmJobRecordDir/stats
                if [ -f $smartSlurmJobRecordDir/stats/$statKey.mem.stat ]; then
                    if [ -f .command.sh ] && [ -f .command.run ]; then 
                        output=`calculateMemTime.sh $program ${ref//\//-} $inputSize $pipeline 2>> ../../../.nextflow.log`
                    else #if [[ "$parentCmd" == */bin/snakemake* ]]; then
                        output=`calculateMemTime.sh $program ${ref//\//-} $inputSize $pipeline 2>> .smartSlurm.log`
                    fi 

                    #resAjust="$resAjust\n`cat $smartSlurmJobRecordDir/stats/$statKey.mem.stat`\n"
                    resAjust="$resAjust\n#Output from calculateMemTime.sh: $output \n"
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
fi

# ---------------------------------------------------------------------------
# Failure floor (issue B2): never submit at or near a proven-insufficient
# allocation. For this (program,ref), find the largest allocation that OOM'd/OOT'd
# at an input size <= this one (by monotonicity, a smaller input that needed more
# proves this larger input needs more too), then escalate it with the SAME helper
# auto-rerun uses. The final request is max(estimate, floor), so initial submit,
# auto-rerun, and a fresh manual rerun can never drop below a known failure.
command -v nextMemAfterFail >/dev/null 2>&1 || { _sjrd="$smartSlurmJobRecordDir"; . "$(dirname "$0")/../config/config.txt"; [ -n "$_sjrd" ] && smartSlurmJobRecordDir="$_sjrd"; }
jr=$smartSlurmJobRecordDir/jobRecord.txt
if [ -f "$jr" ] && echo "$inputSize" | grep -qE '^[0-9]+$'; then
    recentFails=`jrEmitV2 "$jr" 2>/dev/null | awk -F, -v a="$program" -v rre="$refRe" -v t="$(date -d '1 hour ago' +%s)" "$JR_AWK_COLMAP"'$(col["status"])!="COMPLETED" && ($(col["epoch"])+0)>t && $(col["step"])==a && $(col["ref"]) ~ rre' | wc -l`
    maxFailMem=`jrEmitV2 "$jr" 2>/dev/null | awk -F, -v a="$program" -v rre="$refRe" -v sz="$inputSize" "$JR_AWK_COLMAP"'$(col["step"])==a && $(col["ref"]) ~ rre && $(col["inputSize"]) ~ /^[0-9]+$/ && $(col["inputSize"])<=sz && $(col["status"]) ~ /OOM/ {if($(col["memAllocated"])>m)m=$(col["memAllocated"])} END{print m+0}'`
    maxFailTime=`jrEmitV2 "$jr" 2>/dev/null | awk -F, -v a="$program" -v rre="$refRe" -v sz="$inputSize" "$JR_AWK_COLMAP"'$(col["step"])==a && $(col["ref"]) ~ rre && $(col["inputSize"]) ~ /^[0-9]+$/ && $(col["inputSize"])<=sz && $(col["status"]) ~ /OOT/ {if($(col["timeAllocated"])>t)t=$(col["timeAllocated"])} END{print t+0}'`
    if [ "${maxFailMem:-0}" -gt 0 ]; then
        memFloor=`nextMemAfterFail "$maxFailMem" "$recentFails"`
        if [ "${memFloor:-0}" -gt "${mem:-0}" ]; then
            resAjust="$resAjust\n#Failure floor: raising mem from $mem to $memFloor M (proven-insufficient OOM alloc was $maxFailMem M)."
            mem=$memFloor
        fi
    fi
    if [ "${maxFailTime:-0}" -gt 0 ]; then
        timeFloor=`nextTimeAfterFail "$maxFailTime" "$recentFails"`
        if [ "${timeFloor:-0}" -gt "${min:-0}" ]; then
            resAjust="$resAjust\n#Failure floor: raising time from $min to $timeFloor min (proven-insufficient OOT alloc was $maxFailTime min)."
            min=$timeFloor
        fi
    fi
fi

if [[ $adjust == "adjust" ]]; then
    echo $inputSize $mem $min $extraMem >> $smartSlurmLogDir/$flag.adjust
fi  

echo -e "$resAjust" >> $smartSlurmLogDir/$flag.out       

echo $inputSize $mem $min $extraMem

# for testing
#echo $inputSize 8000 1 1

echoerr Got inputsize: $inputSize mem: $mem time: $min extraMem: $extraMem extralTime: $defaultExtraTime
