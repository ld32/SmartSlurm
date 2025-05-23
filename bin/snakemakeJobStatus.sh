#!/usr/bin/env bash

set -euo pipefail

jobid="$1"

# Fall back to squeue if sacct returns nothing (e.g. job hasn't started yet)
sacct_out=$(sacct -j "$jobid" --format=State --noheader 2>/dev/null | grep -Ev '^\s*$' | head -n1 | awk '{print $1}')

if [[ -z "$sacct_out" ]]; then
    # Check with squeue: job is pending/running
    if squeue -j "$jobid" -h | grep -q .; then
        exit 0    # running/pending
    else
        exit 0    # race condition: possibly still starting, treat as running
    fi
fi

# Normalize State string
sacct_out=$(echo "$sacct_out" | tr '[:upper:]' '[:lower:]')

if [[ "$sacct_out" == *"completed"* ]]; then
    exit 2    # Success (finished)
elif [[ "$sacct_out" == *"failed"* ]] || \
     [[ "$sacct_out" == *"cancelled"* ]] || \
     [[ "$sacct_out" == *"timeout"* ]] || \
     [[ "$sacct_out" == *"out_of_memory"* ]]; then
    exit 1    # Failure
else
    exit 0    # Still running or pending
fi