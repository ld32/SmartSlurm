#!/bin/sh
# SmartSlurm test pipeline (legacy blank-line blocks / --wrap).
# Bodies avoid quoted-whitespace args, which the --wrap flattener collapses.
# (For richer text ops -- awk, tr ' ', quoted patterns -- use the #@end version.)
manifest=$1
[ -z "$manifest" ] && { echo "usage: $0 <manifest>"; exit 1; }
samples=$(cut -d' ' -f1 "$manifest")
while read sid rest; do echo "$rest" > phrase.$sid.txt; done < "$manifest"

for s in $samples; do
    in=phrase.$s.txt

    #@1,0,normalize,,in,sbatch -p short -c 1 --mem 100M -t 5:00
    tr A-Z a-z < $in > norm.$s.txt

    #@2,1,tokenize,,in,sbatch -p short -c 1 --mem 100M -t 5:00
    xargs -n1 < norm.$s.txt | sort > tokens.$s.txt

    #@3,2,count,,in,sbatch -p short -c 1 --mem 100M -t 5:00
    uniq -c tokens.$s.txt > counts.$s.txt

done

#@4,3,summary,,,sbatch -p short -c 1 --mem 100M -t 5:00
cat counts.*.txt > all_counts.txt

