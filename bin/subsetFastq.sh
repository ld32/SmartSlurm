#!/bin/sh

usage(){
    echo -e "Usage:\nsubsetFastq.sh <path to original fastq file in .fastq.gz format> <number of reads, for example 500000 or 1000000>"; 
    exit 1
}

inputDir=$1
readCount=$2

[ -d "$inputDir" ] || usage

test -n "$readCount" -a "$readCount" -ge 0 2>/dev/null || usage

rowCount=$(( $readCount * 4 ))

# for i in *.fastq.tar.gz; do 
#     zcat $i  | head -n 2000000 > $cwd/${i%.gz}; done 
# done 

# cd -
# for i in *.fastq; do 
#     tar -cvzf $i.gz $i; 
#     rm $i
# done



for i in $inputDir/*.fastq.gz; do 
    echo Working on $i ... 
    zcat $i  | head -n $rowCount | gzip > ${i##*/}; 
done 