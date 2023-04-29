#!/bin/sh

#set -x 

usage(){
    echo -e "Usage:\nsubsetFastq.sh <path to original fastq file in .fastq.gz format> <number of reads, for example 500000 or 1000000>"; 
    exit 1
}

readCount="${@: -1}"

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

for i in $@; do 
    [ -f $i ] && [[ "$i" == *fastq.gz ]] && echo Working on $i ... && zcat $i  | head -n $rowCount | gzip > ${i##*/}; 
    for j in $i/*.fastq.gz; do 
        [ -f $j ] && [[ "$j" == *fastq.gz ]] && echo Working on $j ...  && zcat $j  | head -n $rowCount | gzip > ${j##*/}; 
    done
done     