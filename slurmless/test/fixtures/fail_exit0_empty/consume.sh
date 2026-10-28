#!/bin/bash
# Fails when handed an empty file — this is the SYMPTOM step.
set -e
[ -s "$1" ] || { echo "consume: input '$1' is empty" >&2; exit 3; }
wc -l "$1"
