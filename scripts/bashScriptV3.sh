#!/bin/sh

number=$1

[ -z "$number" ] && echo -e "Error: number is missing.\nUsage: `basename $0` <numbert>" && exit 1

for i in {1..5}; do
    
    input=numbers$i.txt
    
    #@1,0,findNumber,,input,sbatch -p short -c 1 --mem 4G -t 50:0 
    findNumber.sh $number $input > $number.$i.txt

done

# test what happens if input does not exist whe job is submitted
input=all$number.txt; 
#@2,1,mergeNumber,,input,sbatch -p short -c 1 --mem 4G -t 50:0 
cat $number.*.txt > all$number.txt; findNumber.sh $number $input > final$number.txt 