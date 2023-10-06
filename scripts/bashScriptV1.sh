#!/bin/sh

for i in {1..1}; do
    input=bigText$i.txt
    output=1234.$i.txt
    useMemTimeWithInput.sh $input; grep 1234 $input > $output

    output=5678.$i.txt
    useMemTimeWithInput.sh $input; grep 5678 $input > $output
done

input=bigText1.txt
output=all.txt
useMemTimeWithInput.sh $input; cat 1234.*.txt 5678.*.txt > $output