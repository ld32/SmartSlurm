#!/bin/sh

echoerr() { echo "$@" 1>&2; }

usage() { echoerr -e "Usage: \nestimateMemTime.sh bowtie2 hg19 inputs\nReturn mem in M and time in minutes."; exit 1; }

#set -x  

echoerr Running: estimateMemTime.sh $@

if [ -f ~/.smartSlurm/config.txt ]; then 
    source ~/.smartSlurm/config.txt
else     
    source $(dirname $0)/../config/config.txt || { echo Config list file not found: config.txt; exit 1; }
fi

software=$1
ref=$2 
inputSize=$3
 
#echoerr content of .rcbio
#ls -l ~/.rcbio 1>&2  

echoerr Estimating mem: 
if [ -s $jobRecordDir/stats/$software.$ref.mem.stat ]; then   
    .  $jobRecordDir/stats/$software.$ref.mem.stat # Finala=0.03 Finalb=5.0 Mean=250.0000 Minimum=200.0000 Maximum=300.0000 Median=250.0000 

    #echoerr content: $software.$ref.mem.stat.final:

    #cat ~/.rcbio/$software.$ref.mem.stat.final 1>&2

    Finala=`printf "%.15f\n" $Finala`
    Finalb=`printf "%.15f\n" $Finalb`
    Maximum=`printf "%.15f\n" $Maximum`
    STDFIT=`printf "%.15f\n" $STDFIT`
    
    echoerr Finala: $Finala Finalb: $Finalb Maximum: $Maximum  STDFIT: $STDFIT

    echoerr "mem formula: ( $Finala x $inputSize + $Finalb + $STDFIT * 2 ) x 1.0"

    if (( $(echo "$Maximum + 0.01 > $inputSize" |bc -l) )); then 
        mem=`echo "( $Finala * $inputSize + $Finalb + $STDFIT * 2 ) * 1.0" |bc `
    else
        echoerr outOfRange 
        echo outOfRange
        exit  
    fi
    echoerr
else
    mem=0
fi

echoerr Estimating time: 

if [ -s $jobRecordDir/stats/$software.$ref.time.stat ]; then
    .  $jobRecordDir/stats/$software.$ref.time.stat # Finala=0.03 Finalb=5.0 Mean=250.0000 Minimum=200.0000 Maximum=300.0000 Median=250.0000 

    #echoerr content: $software.$ref.time.stat.final:

    #cat ~/.rcbio/$software.$ref.time.stat.final 1>&2

    Finala=`printf "%.15f\n" $Finala`
    Finalb=`printf "%.15f\n" $Finalb`
    Maximum=`printf "%.15f\n" $Maximum`
    STDFIT=`printf "%.15f\n" $STDFIT`
    
    echoerr Finala: $Finala Finalb: $Finalb Maximum: $Maximum  STDFIT: $STDFIT

    echoerr "time formula: ( $Finala x $inputSize + $Finalb + $STDFIT * 2 ) x 1.0"

    time=`echo "( $Finala * $inputSize + $Finalb + $STDFIT * 2 ) * 1.0" |bc `
    
    #echoerr "time formula: ( $Finala x $inputSize + $Finalb ) x 1.0"
else 
    time=0
fi
echoerr Got  $mem $time

# +1 to round up the number to integer, for example 0.8 becomes 2, 3.5 becomes 5
# memory in M and time in minutes
output="$((${mem%.*} + 1)) $((${time%.*} + 1))" 
      
echo $output