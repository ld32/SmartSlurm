#!/bin/sh
# Reference wiring (NOT shown to the agent): the intended correct order.
for p in $(cut -d' ' -f1 words.txt); do sh lower.sh "$p"; done
sh assemble.sh
