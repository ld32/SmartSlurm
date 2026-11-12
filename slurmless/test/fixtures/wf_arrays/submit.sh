#!/bin/bash
# Unmodified hand-rolled SLURM script using job arrays and array dependencies.
set -u
sub() { sbatch "$@" | awk '{print $4}'; }

# success array: tasks 1..3 all succeed; -o uses %A/%a per-task log names
AR=$(sub  -J arr      --array=1-3 -o arr_%A_%a.out        job.sh arr)
# afterok on the whole array -> waits for ALL tasks, then runs
sub       -J afterall --dependency=afterok:"$AR"          job.sh postall  >/dev/null
# aftercorr array -> task N runs after AR's task N (all succeed -> all run)
sub       -J corr     --array=1-3 --dependency=aftercorr:"$AR" job.sh corr >/dev/null

# failing array: task 2 exits 1
BAD=$(sub -J bad      --array=1-3                         job.sh bad failtask2)
# afterok on a FAILED array -> cancelled (aggregate state is failure)
sub       -J afterbad --dependency=afterok:"$BAD"         job.sh postbad  >/dev/null
# aftercorr on the failed array -> task 2 cancelled; tasks 1 & 3 run
sub       -J corrbad  --array=1-3 --dependency=aftercorr:"$BAD" job.sh corrbad >/dev/null
