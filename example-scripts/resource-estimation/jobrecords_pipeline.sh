#!/bin/sh
# =============================================================================
# rAP EXAMPLE — "jobrecords": resource estimation / auto-rerun testbed
# =============================================================================
# One `load` job per sample; each job's MEMORY and WALL-TIME scale with its
# input size (via sizedLoad.sh). Run it over several sizes so SmartSlurm records
# mem/time at each size and can fit a curve; set FAIL_PROB=1 to force an OOM/OOT
# overshoot and exercise the resource auto-increment + rerun path.
#
# PREREQUISITE: copy sizedLoad.sh into the run directory (the jobs call it).
#
# CLEAN run (all succeed, seeds the jobRecord fit):
#   runAsPipeline --script "jobrecords_pipeline.sh jobrecords_manifest.txt" \
#                 --sbatch-options "sbatch -p short -c 1 --mem 64M -t 2:00" --tmp noTmp
#
# RERUN test (force an OOM on the big sample -> auto-increment + requeue):
#   FAIL_PROB=1 FAIL_MODE=oom OVERSHOOT_X=4 runAsPipeline --script ... (as above)
#
# MANIFEST  (one sample per line:  <id> <rows>)  -- rows is the size knob.
# =============================================================================
manifest=$1
[ -z "$manifest" ] && { echo "usage: $0 <manifest>"; exit 1; }

# knobs (override from your shell); baked into each job at submit time
mkr=${MEM_KB_PER_ROW:-20}; spr=${SEC_PER_ROW:-0.02}; noi=${NOISE_PCT:-5}
fp=${FAIL_PROB:-0}; fm=${FAIL_MODE:-both}; ox=${OVERSHOOT_X:-4}

# prep on the submit host: a data file of N rows per sample (N drives resources)
samples=$(cut -d' ' -f1 "$manifest")
while read sid rows; do seq 1 "$rows" > data.$sid.txt; done < "$manifest"

for s in $samples; do
    in=data.$s.txt
    #@1,0,load,,in,sbatch -p short -c 1 --mem 64M -t 2:00
    MEM_KB_PER_ROW=$mkr SEC_PER_ROW=$spr NOISE_PCT=$noi FAIL_PROB=$fp FAIL_MODE=$fm OVERSHOOT_X=$ox bash sizedLoad.sh $in > load.$s.txt 2>&1
    #@end
done
