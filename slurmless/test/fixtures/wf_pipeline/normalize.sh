#!/bin/bash
# Step 1: normalize this sample's phrase -> lowercase, squeeze/trim spaces.
set -eu
id="$1"
phrase=$(sed -n "s/^$id //p" phrases.txt)
printf '%s\n' "$phrase" | tr 'A-Z' 'a-z' | tr -s ' ' | sed 's/^ *//; s/ *$//' > "norm.$id.txt"
