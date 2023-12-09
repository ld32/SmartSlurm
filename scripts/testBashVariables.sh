#!/bin/sh
#loopStart:i

for i in {1..1}; do
    fastq="fastq/subset.1m.Adelman004_TT029_DMSO_A_S1_L001_R1_001.fastq"

    #@1,0,step1,,,sbatch -p short -t 5
    rowCout=$(wc -l $fastq | cut -d' ' -f1); \
    rowCout=$((rowCout/4/2)); rowCout=$((rowCout*4)); fastq=${fastq%.fastq}; echo $fastq;  echo ${fastq#fa}

done
