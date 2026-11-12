#!/bin/bash
# Submit three jobs, then scancel one before it runs — a script that manages its
# own jobs with scancel. B depends on A, so cancelling A should orphan B; C is
# independent and should still run.
set -u
A=$(sbatch -J A job.sh A 0 | awk '{print $4}')
B=$(sbatch -J B --dependency=afterok:"$A" job.sh B 0 | awk '{print $4}')
C=$(sbatch -J C job.sh C 0 | awk '{print $4}')
scancel "$A"                       # cancel A while it is still pending
echo "$A $B $C" > ids.txt
