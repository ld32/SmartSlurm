#!/bin/sh
# Reference CORRECT rAP pipeline (ground truth) — what a faithful conversion of the
# correct variants should be equivalent to. Fan-out `lower` per word position, then
# a single `assemble` that waits for all of them.
samples="$(cut -d' ' -f1 words.txt)"
for s in $samples; do
    #@1,0,lower,,none,sbatch -p short -c 1 --mem 100M -t 5:00
    lower.sh $s
    #@end
done
#@2,1,assemble,,none,sbatch -p short -c 1 --mem 100M -t 5:00
assemble.sh
#@end
