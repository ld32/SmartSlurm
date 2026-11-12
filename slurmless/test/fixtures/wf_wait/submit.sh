#!/bin/bash
# `sbatch --wait` blocks until the job finishes and exits with the job's code.
# (set -u only: a failing --wait must not abort the script.)
set -u
sbatch -J okp  --wait job.sh okp 0;  echo "okp rc=$?"  >  waitlog.txt
sbatch -J badp --wait job.sh badp 1; echo "badp rc=$?" >> waitlog.txt
