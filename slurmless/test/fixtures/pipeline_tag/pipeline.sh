#!/bin/sh
#@pipeline rnaseq_qc
# This pipeline hardcodes its identity with #@pipeline (not the filename, not a
# per-run flag). The directive must set the pipeline id AND must not be mistaken
# for a step marker -- the two real steps below still convert and run.
number=$1
[ -z "$number" ] && echo "Error: number is missing." && exit 1

input=numbers1.txt
#@1,0,findNumber,,input,sbatch -p short -c 1 --mem 2G -t 50:0
findNumber.sh $number $input > $number.out.txt

#@2,1,mergeNumber,,,sbatch -p short -c 1 --mem 4G -t 50:0
cat $number.out.txt > all$number.txt
