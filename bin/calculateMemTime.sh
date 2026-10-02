#!/bin/sh

echoerr() { echo "$@" 1>&2; }

usage() { echoerr -e "Usage: \ncalculateMemTime.sh bowtie2 hg19 inputs\nReturn mem in M and time in minutes."; exit 1; }

#set -x  

echoerr Running: $0 $@

# if [ -z "$smartSlurmJobRecordDir" ]; then 
#     if [ -f ~/.smartSlurm/config/config.txt ]; then
#         source ~/.smartSlurm/config/config.txt
#     else
#         source $(dirname $0)/../config/config.txt || { echoerr Config list file not found: config.txt; exit 1; }
#     fi
# fi   

software=$1
ref=$2 
inputSize=$3
 

echoerr Calculating mem:
memFormu=memFormu:
if [ -s $smartSlurmJobRecordDir/stats/$software.$ref.mem.stat ]; then   
    
    unset Finala Finalb Maximum STDFIT SCount
    .  $smartSlurmJobRecordDir/stats/$software.$ref.mem.stat # Finala=0.03 Finalb=5.0 Mean=250.0000 Minimum=200.0000 Maximum=300.0000 Median=250.0000 

    #echoerr content: $software.$ref.mem.stat.final:

    #cat ~/.rcbio/$software.$ref.mem.stat.final 1>&2

    # Robust numeric defaults: a missing/empty stat key must not break arithmetic
    # (e.g. STDFIT is empty for perfectly collinear history).
    : "${Finala:=0}" "${Finalb:=0}" "${Maximum:=0}" "${STDFIT:=0}" "${SCount:=0}"

    echoerr Finala: $Finala Finalb: $Finalb Maximum: $Maximum  STDFIT: $STDFIT

    echoerr "mem formula: ( $Finala x $inputSize + $Finalb + $STDFIT * 2 ) x 1.0"

    # Compute within the observed range, or extrapolate only with enough data;
    # otherwise report outOfRange and let the caller use defaults.
    if (( $(echo "$Maximum + 0.01 >= $inputSize" |bc -l) )) || (( $(echo "$SCount > 9" |bc -l) )); then
        # Round and floor to 100 M with awk. awk (not bc + ${mem%.*} + [[ -le ]])
        # because bc prints sub-1 negatives as "-.3", which ${mem%.*} turns into
        # "-", breaking the integer test and silently passing garbage downstream.
        # A negative/low prediction (flat or noise-driven slope) floors cleanly.
        mem=`awk -v a="$Finala" -v x="$inputSize" -v b="$Finalb" -v s="$STDFIT" \
             'BEGIN{ v = a*x + b + 2*s; if (v < 100) v = 100; printf "%d", v + 0.5 }'`
        memFormu=$memFormu${Finala}X${inputSize}+$Finalb+$STDFIT*2
    else
        echoerr outOfRange
        echo outOfRange
        exit
    fi
    echoerr
else
    mem=0
fi

echoerr Got  $mem 

echoerr Calculating time: 
timeFormu=timeFormu:

if [ -s $smartSlurmJobRecordDir/stats/$software.$ref.time.stat ]; then
    
    unset Finala Finalb Maximum STDFIT SCount
    .  $smartSlurmJobRecordDir/stats/$software.$ref.time.stat # Finala=0.03 Finalb=5.0 Mean=250.0000 Minimum=200.0000 Maximum=300.0000 Median=250.0000 

    #echoerr content: $software.$ref.time.stat.final:

    #cat ~/.rcbio/$software.$ref.time.stat.final 1>&2

    : "${Finala:=0}" "${Finalb:=0}" "${Maximum:=0}" "${STDFIT:=0}" "${SCount:=0}" "${RSquare:=0}"

    echoerr Finala: $Finala Finalb: $Finalb Maximum: $Maximum  STDFIT: $STDFIT

    echoerr "time formula: ( $Finala x $inputSize + $Finalb + $STDFIT * 2 ) x 1.0"

    # Round and floor to 10 min with awk (same robustness as the mem branch).
    time=`awk -v a="$Finala" -v x="$inputSize" -v b="$Finalb" -v s="$STDFIT" \
          'BEGIN{ v = a*x + b + 2*s; if (v < 10) v = 10; printf "%d", v + 0.5 }'`
    timeFormu=$timeFormu:${Finala}X${inputSize}+$Finalb+$STDFIT*2
else
    time=0
fi
echoerr Got $time

# +1 to round up the number to integer, for example 0.8 becomes 2, 3.5 becomes 5
# memory in M and time in minutes
: "${RSquare:=0}"
output="$mem $time \n\n$memFormu\n$timeFormu\nRSquare=$RSquare" 
      
echo $output