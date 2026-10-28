#!/bin/sh
# Per-sample outer loop; an inner per-chromosome loop makes step 2 a FAN-OUT:
# one (sample,step) grid cell backs multiple jobs -> should render as '@'.
number=$1
for s in {1..2}; do
    input=numbers$s.txt

    #@1,0,prep,,input,sbatch -p short -c 1 --mem 2G -t 50:0
    findNumber.sh $number $input > $number.$s.txt

    for c in {1..3}; do
        #@2,1,callChrom,,,sbatch -p short -c 1 --mem 2G -t 50:0
        echo "sample $s chrom $c" > $number.$s.$c.chrom.txt

    done

done

#@3,2,summarize,,,sbatch -p short -c 1 --mem 4G -t 50:0
cat $number.*.chrom.txt > summary$number.txt

