#!/bin/bash
# Generic payload: drop a marker proving this job RAN, then exit with the given
# code (default 0). Used to build dependency-grammar scenarios.
name="$1"; code="${2:-0}"
echo "ran $name as job ${SLURM_JOB_ID:-?}" > "ran.$name"
exit "$code"
