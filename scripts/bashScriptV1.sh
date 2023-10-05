#!/bin/sh

for i in {1..1}; do
    input=bigText$i.txt
    output=1234.$i.txt
    useSomeMemTimeAccordingInputSize.sh $input; grep 1234 $input > $output

    output=5678.$i.txt
    useSomeMemTimeAccordingInputSize.sh $input; grep 5678 $input > $output
done

input=bigText1.txt
output=all.txt
useSomeMemTimeAccordingInputSize.sh $input; cat 1234.*.txt 5678.*.txt > $output