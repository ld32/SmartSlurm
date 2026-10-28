#!/bin/sh
number=$1
for s in 1 2; do
    in=reads$s.txt

    #@1,0,trim,,in,sbatch -p short -c 1 --mem 2G -t 50:0
    cp $in trim.$number.$s.txt

    #@2,1,count,,in,sbatch -p short -c 1 --mem 2G -t 50:0
    wc -l trim.$number.$s.txt > count.$number.$s.txt

done

#@3,2,merge,,,sbatch -p short -c 1 --mem 2G -t 50:0
cat count.$number.*.txt > merged.$number.txt

