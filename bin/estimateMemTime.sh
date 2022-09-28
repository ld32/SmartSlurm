#!/bin/sh

echoerr() { echo "$@" 1>&2; }

usage() { echoerr -e "Usage: \nestimateMemTime.sh bowtie2 hg19 inputs\nReturn mem in M and time in minutes."; exit 1; }

#set -x 
echoerr Running: estimateMemTime.sh $@

software=$1
ref=$2 
inputSize=$3


#echoerr content of .rcbio
#ls -l ~/.rcbio 1>&2  

echoerr Estimating mem:
        
.  ~/.smartSlurm/$software.$ref.mem.stat.final # Finala=0.03 Finalb=5.0 Mean=250.0000 Minimum=200.0000 Maximum=300.0000 Median=250.0000 

#echoerr content: $software.$ref.mem.stat.final:

#cat ~/.rcbio/$software.$ref.mem.stat.final 1>&2

Finala=`printf "%.15f\n" $Finala`
Finalb=`printf "%.15f\n" $Finalb`
Maximum=`printf "%.15f\n" $Maximum`
echoerr Finala: $Finala Finalb: $Finalb Maximum: $Maximum

echoerr "mem formula: ( $Finala x $inputSize + $Finalb ) x 1.5 + 10"

if (( $(echo "$Maximum > $inputSize" |bc -l) )); then 
    mem=`echo "( $Finala * $inputSize + $Finalb ) * 1.5 + 10" |bc ` # at least 10M
else
    echoerr outOfRange 
    echo outOfRange
    exit  
fi
echoerr
echoerr Estimating time: 

.  ~/.smartSlurm/$software.$ref.time.stat.final # Finala=0.03 Finalb=5.0 Mean=250.0000 Minimum=200.0000 Maximum=300.0000 Median=250.0000 

#echoerr content: $software.$ref.time.stat.final:

#cat ~/.rcbio/$software.$ref.time.stat.final 1>&2

Finala=`printf "%.15f\n" $Finala`
Finalb=`printf "%.15f\n" $Finalb`
Maximum=`printf "%.15f\n" $Maximum`
echoerr Finala: $Finala Finalb: $Finalb Maximum: $Maximum

time=`echo "( $Finala * $inputSize + $Finalb ) * 1.5 + 10" |bc` # at least 10 minutes 

echoerr "time formula: ( $Finala x $inputSize + $Finalb ) x 1.5 + 10"

echoerr Got  $mem $time

# +1 to round up the number to integer
# memory in M and time in minutes
echo $((${mem%.*} + 1)) $((${time%.*} + 1)) 
      
