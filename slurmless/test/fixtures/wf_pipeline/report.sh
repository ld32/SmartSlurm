#!/bin/bash
# Step 4 (afterok: merge): format the final report.
set -eu
printf 'samples: %s\n' "$(ls wc.*.txt | wc -l | tr -d ' ')" > report.txt
printf 'total words: %s\n' "$(cat grandtotal.txt)" >> report.txt
