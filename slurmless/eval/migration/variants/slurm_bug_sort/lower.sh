#!/bin/sh
p="$1"
w=$(grep "^$p " words.txt | cut -d' ' -f2-)
printf '%s\n' "$w" | tr 'A-Z' 'a-z' > "w.$p.txt"
