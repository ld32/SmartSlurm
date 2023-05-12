#!/bin/sh

#set -x 
usage(){
    echo -e "Usage:\nsubsetFastq.sh <path to original .fastq.gz file, such as: abc/*.fastq.gz, or the folder containing .fastq.gz files, such as: abc/fastq> <number of reads, such as: 500000 or 1000000>"; 
    exit 1
}

readCount="${@: -1}"

test -n "$readCount" -a "$readCount" -ge 0 2>/dev/null || usage

rowCount=$(( $readCount * 4000000 ))
million=`echo "scale=6;$readCount/1000000"|bc| sed 's/^\./0./' | sed 's/0\{1,\}$//'`
[[ "$million" == 0* ]] && million="$million."
folder=fastq.${million}subset

mkdir -p $folder
for i in $@; do 
    if [ -f $i ] && [[ "$j" == *fastq.gz ]]; then 
        file=$folder/subset.$million${i##*/}
        [ -f "$file" ] && echo Current working folder already has file $file Skipping it... && continue 
        echo Working on $j ...  && zcat $j  | head -n $rowCount | gzip > $file;    
    elif [ -d $i ]; then 
        for j in $i/*.fastq.gz; do 
            file=$folder/subset.$million${j##*/}
            [ -f "$file" ] && echo Current working folder already has file $file Skipping it... && continue 
            echo Working on $j ...  && zcat $j  | head -n $rowCount | gzip > $file; 
        done
    fi    
done     
