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

# Target first-pass completion rate (issue B1). Used as the percentile for the
# size-independent (constant) model. Overridable via config/env; default 93%.
: "${targetFirstPassPct:=93}"

# pctile <datafile "size value" per line> <percentile> -> integer p-th percentile
# of the value column (rounded). Empty output if the file has no usable rows.
pctile() {
    cut -d' ' -f2 "$1" 2>/dev/null | sort -n | awk -v p="$2" '
        { v[NR] = $1 }
        END { if (NR == 0) exit
              i = int(p/100*NR + 0.9999); if (i < 1) i = 1; if (i > NR) i = NR
              printf "%d", v[i] + 0.5 }'
}

# relMarginFrac <datafile "size value"> <slope a> <intercept b> <percentile>
# -> p-th percentile of the RELATIVE residuals (value-(a*size+b))/(a*size+b),
# clamped >= 0, as a decimal fraction. A distribution-free, scale-invariant
# upper-prediction margin for the line model (issue B1): the estimate is then
# (a*x+b)*(1+frac), so ~p%% of observed jobs fall under it regardless of size.
# Relative (not absolute) because resource-usage noise is typically multiplicative
# (bigger jobs vary more in absolute terms). Empty if the file has no usable rows.
relMarginFrac() {
    awk -v a="$2" -v b="$3" '{ f=a*$1+b; if (f>0) print ($2-f)/f }' "$1" 2>/dev/null | sort -n | awk -v p="$4" '
        { r[NR] = $1 }
        END { if (NR == 0) exit
              i = int(p/100*NR + 0.9999); if (i < 1) i = 1; if (i > NR) i = NR
              m = r[i]; if (m < 0) m = 0; printf "%.6f", m }'
}

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

    # Issue A1/A2: choose the model per dimension.
    #  - Positive slope  -> the fitted LINE captures real size dependence (and is
    #    monotone non-decreasing, so no bigger-input-uses-less inversion).
    #  - Slope <= 0       -> a non-positive slope cannot be trusted to extrapolate
    #    and usually means the resource is size-INDEPENDENT (e.g. constant memory,
    #    only time scales) or the slope is just scheduler/perf noise. Fall back to
    #    a size-independent CONSTANT model: the target percentile of observed usage.
    if awk -v a="$Finala" 'BEGIN{ exit !(a > 0) }'; then
        # Scale-invariant upper-prediction margin targeting targetFirstPassPct (B1):
        # estimate = line*(1+frac), frac = p-th percentile of relative residuals.
        # Fall back to an STDFIT-derived fraction if the raw data is gone.
        frac=`relMarginFrac "$smartSlurmJobRecordDir/stats/$software.$ref.mem" "$Finala" "$Finalb" "$targetFirstPassPct"`
        [ -z "$frac" ] && frac=`awk -v s="$STDFIT" -v a="$Finala" -v x="$inputSize" -v b="$Finalb" 'BEGIN{ d=a*x+b; if(d<=0){print 0}else{printf "%.6f", 2*s/d} }'`
        echoerr "mem model: LINE ( $Finala x $inputSize + $Finalb ) * (1 + frac[$frac] @p${targetFirstPassPct})"
        # Within observed range, or extrapolate only with enough data; else outOfRange.
        if (( $(echo "$Maximum + 0.01 >= $inputSize" |bc -l) )) || (( $(echo "$SCount > 9" |bc -l) )); then
            # Round and floor to 100 M with awk (robust to negative/sub-1 results).
            mem=`awk -v a="$Finala" -v x="$inputSize" -v b="$Finalb" -v fr="$frac" \
                 'BEGIN{ v = (a*x + b) * (1 + fr); if (v < 100) v = 100; printf "%d", v + 0.5 }'`
            memFormu=$memFormu${Finala}X${inputSize}+$Finalb,p${targetFirstPassPct}*$frac
        else
            echoerr outOfRange
            echo outOfRange
            exit
        fi
    else
        # Constant model: p-th percentile of observed memUsed (size-independent, so
        # no outOfRange guard). Falls back to Maximum if the raw data file is gone.
        mem=`pctile "$smartSlurmJobRecordDir/stats/$software.$ref.mem" "$targetFirstPassPct"`
        [ -z "$mem" ] && mem=$Maximum
        mem=`awk -v v="$mem" 'BEGIN{ if (v < 100) v = 100; printf "%d", v + 0.5 }'`
        memFormu=${memFormu}const.p${targetFirstPassPct}=${mem}
        echoerr "mem model: CONSTANT (slope $Finala <= 0); p${targetFirstPassPct} of observed = $mem"
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

    # Issue A1/A2: same per-dimension model choice as mem. Positive slope -> line;
    # otherwise size-independent constant model (p-th percentile of observed time).
    if awk -v a="$Finala" 'BEGIN{ exit !(a > 0) }'; then
        frac=`relMarginFrac "$smartSlurmJobRecordDir/stats/$software.$ref.time" "$Finala" "$Finalb" "$targetFirstPassPct"`
        [ -z "$frac" ] && frac=`awk -v s="$STDFIT" -v a="$Finala" -v x="$inputSize" -v b="$Finalb" 'BEGIN{ d=a*x+b; if(d<=0){print 0}else{printf "%.6f", 2*s/d} }'`
        echoerr "time model: LINE ( $Finala x $inputSize + $Finalb ) * (1 + frac[$frac] @p${targetFirstPassPct})"
        time=`awk -v a="$Finala" -v x="$inputSize" -v b="$Finalb" -v fr="$frac" \
              'BEGIN{ v = (a*x + b) * (1 + fr); if (v < 10) v = 10; printf "%d", v + 0.5 }'`
        timeFormu=$timeFormu:${Finala}X${inputSize}+$Finalb,p${targetFirstPassPct}*$frac
    else
        time=`pctile "$smartSlurmJobRecordDir/stats/$software.$ref.time" "$targetFirstPassPct"`
        [ -z "$time" ] && time=$Maximum
        time=`awk -v v="$time" 'BEGIN{ if (v < 10) v = 10; printf "%d", v + 0.5 }'`
        timeFormu=${timeFormu}:const.p${targetFirstPassPct}=${time}
        echoerr "time model: CONSTANT (slope $Finala <= 0); p${targetFirstPassPct} of observed = $time"
    fi
else
    time=0
fi
echoerr Got $time

# +1 to round up the number to integer, for example 0.8 becomes 2, 3.5 becomes 5
# memory in M and time in minutes
: "${RSquare:=0}"
output="$mem $time \n\n$memFormu\n$timeFormu\nRSquare=$RSquare" 
      
echo $output