#!/bin/sh

number=$1

[ -z "$number" ] && echo -e "Error: number is missing.\nUsage: `basename $0` <numbert>" && exit 1

for i in {1..1}; do
    
    input=numbers$i.txt
    
    #@1,0,findNumber.Checkpoint,,input,sbatch --qos=testbump -p short -c 1 --mem 2G -t 50:0 
    findNumber.sh $number $input > $number.$i.txt; sleep 1200;

done

##@2,1,mergeNumber,,,sbatch -p short -c 1 --mem 4G -t 50:0 
#cat $number.*.txt > all$number.txt