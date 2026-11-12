#!/bin/bash
# stepB: the downstream job. It reads stepA's marker — so it can only succeed if
# the afterok dependency was honored (stepA ran to completion first). If the
# scheduler ran B before A (or concurrently), marker.txt would be missing and B
# would fail — which is exactly what the test guards against.
set -eu
cat marker.txt > fromB.txt
echo "B-done (job $SLURM_JOB_ID)" >> fromB.txt
