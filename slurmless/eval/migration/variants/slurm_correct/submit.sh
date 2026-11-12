#!/bin/sh
# Hand-rolled SLURM submission: one `lower` job per word position, then an
# `assemble` job that waits (afterok) for all of them.
set -u
deps=""
for p in $(cut -d' ' -f1 words.txt); do
  j=$(sbatch -J "lower_$p" -o "lower_$p.out" lower.sh "$p" | awk '{print $4}')
  deps="$deps:$j"
done
sbatch -J assemble -o assemble.out --dependency=afterok"$deps" assemble.sh
