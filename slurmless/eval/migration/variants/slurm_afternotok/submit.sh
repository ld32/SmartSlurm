#!/bin/sh
# Hand-rolled SLURM with a failure-recovery branch: fetch each word; if fetch of
# position 5 FAILS, an afternotok cleanup recovers it; assemble waits for the good
# fetches AND the recovery.
set -u
deps=""; f5=""
for p in $(cut -d' ' -f1 words.txt); do
  j=$(sbatch -J "fetch_$p" -o "fetch_$p.out" fetch.sh "$p" | awk '{print $4}')
  if [ "$p" = 5 ]; then f5="$j"; else deps="$deps:$j"; fi
done
r=$(sbatch -J recover_5 -o recover_5.out --dependency=afternotok:"$f5" recover.sh 5 | awk '{print $4}')
sbatch -J assemble -o assemble.out --dependency=afterok"$deps":"$r" assemble.sh
