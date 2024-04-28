#!/bin/sh

number=$1

[ -z "$number" ] && echo -e "Error: number is missing.\nUsage: bashScript <numbert>" && exit 1

for i in {1..3}; do
    
    input=numbers$i.txt
    
    #@1,0,findNumber,,input,sbatch -p short -c 1 --mem 2G -t 50:0 
    findNumber.sh 1234 $input > $number.$i.txt
  
done

#@2,1,mergeNumber,,,sbatch -p short -c 1 --mem 2G -t 50:0 
cat $number.*.txt > all$number.txt