#!/bin/bash
# minimal tool: grep a number out of an input file (writes to stdout)
set -e
[ -z "$2" ] && { echo "Error: usage findNumber.sh <n> <file>" >&2; exit 1; }
grep "$1" "$2"
