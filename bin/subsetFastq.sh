#!/bin/sh

#set -x 
usage(){
    echo -e "Usage:\nsubsetFastq.sh <path to original .fastq.gz file, such as: abc/x.fastq.gz, or the folder containing .fastq.gz files, such as: abc/fastq> <number of reads, such as: 500000 or 1000000>"; 
    exit 1
}

readCount="${@: -1}"

test -n "$readCount" -a "$readCount" -ge 0 2>/dev/null || usage

rowCount=$(( $readCount * 4 ))

for i in $@; do 
    [ -f $i ] && [[ "$i" == *fastq.gz ]] && echo Working on $i ... && zcat $i  | head -n $rowCount | gzip > ${i##*/}; 
    for j in $i/*.fastq.gz; do 
        [ -f $j ] && [[ "$j" == *fastq.gz ]] && echo Working on $j ...  && zcat $j  | head -n $rowCount | gzip > ${j##*/}; 
    done
done     