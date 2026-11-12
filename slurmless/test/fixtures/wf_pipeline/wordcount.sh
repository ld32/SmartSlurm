#!/bin/bash
# Step 2 (afterok: normalize): count words in the normalized phrase.
set -eu
id="$1"
wc -w < "norm.$id.txt" | tr -d ' ' > "wc.$id.txt"
