#!/bin/sh

outputs=""
for i in {1..2}; do
    input=bigText$i.txt
    output=1234.$i.txt
    #@1,0,useMemTimeWithInput.sh,i,,input,sbatch -p short -c 1 --mem 2G -t 50:0 
    useMemTimeWithInput.sh $input; grep 1234 $input > $output
    outputs=$outputs,$output
    
    output=5678.$i.txt
    #@2,0,useMemTimeWithInput.sh,i,,input,sbatch -p short -c 1 --mem 2G -t 50:0
    useMemTimeWithInput.sh $input; grep 5678 $input > $output
    outputs=$outputs,$output
done

input1=bigText1.txt
output=all.txt
#@3,1.2,useMemTimeWithInput.sh,,,input
useMemTimeWithInput.sh $input; cat 1234.*.txt 5678.*.txt > $output

