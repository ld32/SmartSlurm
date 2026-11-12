#!/bin/sh
# Per-word step: look up this position's word and lower-case it.
p="$1"
w=$(grep "^$p " words.txt | cut -d' ' -f2-)
printf '%s\n' "$w" | tr 'A-Z' 'a-z' > "w.$p.txt"
