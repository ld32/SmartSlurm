#!/bin/sh
n=$(ls w.*.txt | wc -l)
out=""; p=1
while [ "$p" -le "$n" ]; do out="$out $(cat "w.$p.txt")"; p=$((p + 1)); done
printf '%s\n' "${out# }" > quote.txt
