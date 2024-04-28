#!/bin/sh

number=$1

[ -z "$number" ] && echo -e "Error: number is missing.\nUsage: bashScript <numbert>" && exit 1

for i in {1..5}; do
    
    input=numbers$i.txt
    
    findNumber.sh 1234 $input > $number.$i.txt
  
done

cat $number.*.txt > all$number.txt