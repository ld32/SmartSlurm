#!/bin/bash
# Unmodified hand-rolled SLURM script exercising the full --dependency grammar.
# Base jobs: A succeeds, B fails (exit 1), G succeeds. Everything else depends on
# them in various ways. (set -u, not -e: payloads run later, at drain time.)
set -u
sub() { sbatch "$@" | awk '{print $4}'; }   # submit, echo the captured job id

A=$(sub -J A -o A.out job.sh A 0)            # succeeds
B=$(sub -J B -o B.out job.sh B 1)            # fails
G=$(sub -J G -o G.out job.sh G 0)            # succeeds

# afterany: runs once B terminates, regardless of its exit state
sub -J C --dependency=afterany:"$B"             job.sh C 0 >/dev/null
# afternotok: cleanup that fires because B FAILED
sub -J D --dependency=afternotok:"$B"           job.sh D 0 >/dev/null
# afternotok on a SUCCEEDED job -> can never be satisfied -> cancelled
sub -J E --dependency=afternotok:"$A"           job.sh E 0 >/dev/null
# AND within one term: both A and G must succeed
sub -J F --dependency=afterok:"$A":"$G"         job.sh F 0 >/dev/null
# OR ('?'): A ok OR B ok -> satisfied by A even though B failed
sub -J H --dependency="afterok:$A?afterok:$B"   job.sh H 0 >/dev/null
# afterok on a FAILED job, hold instead of kill -> stays pending (HELD)
sub -J I --kill-on-invalid-dep=no --dependency=afterok:"$B" job.sh I 0 >/dev/null
# afterok on a FAILED job, default -> cancelled
sub -J J --dependency=afterok:"$B"              job.sh J 0 >/dev/null
# after: runs once B has started/terminated (+time offset is a no-op here)
sub -J K --dependency=after:"$B"+5              job.sh K 0 >/dev/null
# singleton: the second 'solo' waits for the first to finish
S1=$(sub -J solo -o S1.out job.sh solo1 0)
sub -J solo --dependency=singleton              job.sh solo2 0 >/dev/null
