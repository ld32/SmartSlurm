#!/bin/bash
# Generic array/step payload. Writes a per-task marker "<tag>.<task>.txt" (task 0
# for a non-array job), then exits 0 — unless told to fail a specific task.
#   usage: job.sh <tag> [failtask2]
tag="$1"; rule="${2:-}"
t="${SLURM_ARRAY_TASK_ID:-0}"
echo "$tag task $t (job ${SLURM_JOB_ID:-?})" > "$tag.$t.txt"
[ "$rule" = failtask2 ] && [ "$t" = 2 ] && exit 1
exit 0
