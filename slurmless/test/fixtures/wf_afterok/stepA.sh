#!/bin/bash
# stepA: the upstream job. Writes a marker that stepB depends on.
echo "A-was-here (job $SLURM_JOB_ID)" > marker.txt
