#!/bin/sh
# Clean 2-step pipeline: two parallel findNumber jobs feed one mergeNumber job.
# Blocks are properly terminated by a blank line before the loop's `done`.
number=$1
[ -z "$number" ] && echo "Error: number is missing." && exit 1

for i in {1..2}; do
    input=numbers$i.txt
    #@1,0,findNumber,,input,sbatch -p short -c 1 --mem 2G -t 50:0
    findNumber.sh $number $input > $number.$i.txt

done

#@2,1,mergeNumber,,,sbatch -p short -c 1 --mem 4G -t 50:0
cat $number.*.txt > all$number.txt

