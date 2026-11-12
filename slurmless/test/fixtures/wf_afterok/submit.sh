#!/bin/bash
# A hand-rolled SLURM submission script — the kind the Migration Assistant has to
# run UNMODIFIED in the sandbox. It uses nothing but stock sbatch: submit stepA,
# capture its job id, then submit stepB with an afterok dependency on stepA.
set -eu

a=$(sbatch -J stepA -o a.out stepA.sh | awk '{print $4}')
echo "submitted stepA as $a"

b=$(sbatch -J stepB -o b.out --dependency=afterok:"$a" stepB.sh | awk '{print $4}')
echo "submitted stepB as $b (afterok:$a)"
