#!/bin/sh

#set -x 
usage(){
    echo -e "Usage:\nsubsetFastq.sh <path to original .fastq.gz file, such as: abc/*.fastq.gz, or the folder containing .fastq.gz files, such as: abc/fastq> <number of reads, such as: 500000 or 1000000>"; 
    exit 1
}

readCount="${@: -1}"

test -n "$readCount" -a "$readCount" -ge 0 2>/dev/null || usage

rowCount=$(( $readCount * 4 ))
million=`echo "scale=6;$readCount/1000000"|bc| sed 's/^\./0./' | sed 's/0\{1,\}$//'`

million=${million%.}
#[[ "$million" == 0* ]] && million="$million."
folder=fastq.${million}m.subset

mkdir -p $folder
for i in $@; do 
    if [ -f $i ] && [[ "$i" == *fastq.gz ]]; then 
        file=$folder/subset.${million}m.${i##*/}
        [ -f "$file" ] && echo Current working folder already has file $file Skipping it... && continue 
        echo Working on $i ...  && zcat $i  | head -n $rowCount | gzip > $file;    
    elif [ -d $i ]; then 
        for j in $i/*.fastq.gz; do 
            file=$folder/subset.${million}m.${j##*/}
            [ -f "$file" ] && echo Current working folder already has file $file Skipping it... && continue 
            echo Working on $j ...  && zcat $j  | head -n $rowCount | gzip > $file; 
        done
    fi  

    if [ -f $i ] && [[ "$i" == *fastq.bz2 ]]; then 
        file=$folder/subset.${million}m.${i##*/}
        [ -f "$file" ] && echo Current working folder already has file $file Skipping it... && continue 
        echo Working on $j ...  && bzcat $i  | head -n $rowCount | bzip2 > $file;    
    elif [ -d $i ]; then 
        for j in $i/*.fastq.bz2; do 
            file=$folder/subset.${million}m.${j##*/}
            [ -f "$file" ] && echo Current working folder already has file $file Skipping it... && continue 
            echo Working on $j ...  && bzcat $j  | head -n $rowCount | bzip2 > $file; 
        done
    fi   
done

