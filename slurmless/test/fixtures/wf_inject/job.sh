#!/bin/bash
# Marker payload. usage: job.sh <name> [code]
name="$1"; code="${2:-0}"
echo "ran $name" > "ran.$name"
exit "$code"
