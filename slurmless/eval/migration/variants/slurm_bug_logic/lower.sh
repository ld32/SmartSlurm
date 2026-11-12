#!/bin/sh
# BUG (logic): writes the word WITHOUT lower-casing it, so the phrase comes out in
# capitals. A faithful port preserves this.
p="$1"
w=$(grep "^$p " words.txt | cut -d' ' -f2-)
printf '%s\n' "$w" > "w.$p.txt"
