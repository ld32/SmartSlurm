#!/bin/bash
# Step 3 (afterok: ALL wordcount jobs — fan-in): gather per-sample counts and sum.
set -eu
cat wc.*.txt > allcounts.txt
awk '{s+=$1} END{print s+0}' allcounts.txt > grandtotal.txt
