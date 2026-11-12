#!/bin/bash
# A typical linear workflow with a cleanup hook. The harness injects a failure
# into 'align' (resource-sim: e.g. OOM) to exercise failure handling:
#   prep -> align -> post   (afterok chain)
#   cleanup                 (afternotok:align -> runs only if align fails)
set -u
sub() { sbatch "$@" | awk '{print $4}'; }
A=$(sub -J prep                                   job.sh prep 0)
B=$(sub -J align   --dependency=afterok:"$A"      job.sh align 0)
C=$(sub -J post    --dependency=afterok:"$B"      job.sh post 0)
D=$(sub -J cleanup --dependency=afternotok:"$B"   job.sh cleanup 0)
