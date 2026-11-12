#!/bin/sh
# Fetch one word. Position 5 deliberately FAILS (simulating a flaky source) and
# writes nothing -- the cleanup step recovers it.
p="$1"
w=$(grep "^$p " words.txt | cut -d' ' -f2-)
if [ "$p" = 5 ]; then echo "fetch $p failed" >&2; exit 1; fi
printf '%s\n' "$w" | tr 'A-Z' 'a-z' > "w.$p.txt"
