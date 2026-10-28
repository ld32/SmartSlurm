#!/bin/sh
# Reproduces "symptom != cause": step 1 cannot read its input, writes an EMPTY
# file, but exits 0. Step 2 consumes the empty file and fails. Current rAP marks
# step 1 SUCCESS and points the user at step 2 — the wrong step. The overhaul
# (M1 producer-check + M2 tracer) will flip step 1 to FAILED.
number=$1

#@1,0,makeTable,,,sbatch -p short -c 1 --mem 2G -t 50:0
badProducer.sh > all$number.txt

#@2,1,consume,,,sbatch -p short -c 1 --mem 2G -t 50:0
consume.sh all$number.txt

