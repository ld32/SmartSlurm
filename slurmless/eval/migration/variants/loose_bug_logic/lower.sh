#!/bin/sh
# BUG (logic): no lower-case transform -> all-caps phrase.
p="$1"
w=$(grep "^$p " words.txt | cut -d' ' -f2-)
printf '%s\n' "$w" > "w.$p.txt"
