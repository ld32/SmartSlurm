#!/bin/sh
number=$1
for s in {1..2}; do
    input=numbers$s.txt
    #@1,0,align,,input,sbatch -p short -c 1 --mem 2G -t 50:0
    findNumber.sh $number $input > $number.$s.txt

done
# a later step that loops over GROUP ids disjoint from sample ids -> axis differs
for g in {7..8}; do
    #@2,1,perGroup,,,sbatch -p short -c 1 --mem 2G -t 50:0
    echo "group $g" > group$g.txt

done
