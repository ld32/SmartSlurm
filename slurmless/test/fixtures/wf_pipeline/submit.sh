#!/bin/bash
# Hand-rolled SLURM driver for a trivial text-processing pipeline. Nothing here
# knows about SmartSlurm or the harness — stock sbatch + afterok dependencies.
#
#   per sample s (from phrases.txt):  normalize(s) --afterok--> wordcount(s)
#   merge  --afterok--> ALL wordcount jobs          (fan-in)
#   report --afterok--> merge
set -u
OPTS="-p short -c 1 --mem 100M -t 5:00"
samples=$(cut -d' ' -f1 phrases.txt)

merge_deps=""
for s in $samples; do
  j1=$(sbatch $OPTS -J "norm_$s" -o "norm_$s.out" normalize.sh "$s" | awk '{print $4}')
  j2=$(sbatch $OPTS -J "wc_$s"   -o "wc_$s.out" --dependency=afterok:"$j1" wordcount.sh "$s" | awk '{print $4}')
  merge_deps="$merge_deps:$j2"
done

# fan-in: merge waits for every per-sample wordcount job
jm=$(sbatch $OPTS -J merge  -o merge.out  --dependency=afterok"$merge_deps" merge.sh  | awk '{print $4}')
# linear tail
sbatch $OPTS -J report -o report.out --dependency=afterok:"$jm" report.sh
