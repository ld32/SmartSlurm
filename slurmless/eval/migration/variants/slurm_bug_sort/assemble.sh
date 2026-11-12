#!/bin/sh
# BUG (ordering): assembles the words ALPHABETICALLY instead of by position, so the
# phrase comes out scrambled. A faithful port must reproduce this (wrong) behavior.
cat w.*.txt | sort | tr '\n' ' ' | sed 's/ $//' > quote.txt
printf '\n' >> quote.txt
