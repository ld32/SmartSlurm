#!/bin/sh

outputs=""
for i in {1..1}; do
    input=bigText$i.txt
    output=1234.$i.txt
    #@1,0,useSomeMemTimeAccordingInputSize.sh,i,,input,sbatch -p short -c 1 --mem 1M -t 2:0:0
    ls  $input; grep 1234 $input > $output
    ##outputs=$outputs,$output

    output=5678.$i.txt
    input=bigTexta$i.txt
    #@2,1,useSomeMemTimeAccordingInputSize.sh,i,,input,sbatch -p short -c 1 --mem 10M -t 6:0
    ls useSomeMemTimeAccordingInputSize.sh $input

done
