#!/bin/bash
# ---------------------------------------------------------------------------
# mock memCpuMonitor.sh  --  replay model-derived resource usage, Slurm-free.
#
# The real bin/memCpuMonitor.sh samples the live process tree every 60s and
# writes $smartSlurmLogDir/job_<JOBID>.memCPU.txt, whose:
#   - line 1 field 6 = START epoch        (cleanUp.sh derives elapsed MIN from it)
#   - every line field 2 = mem used (M)   (cleanUp.sh takes the peak as memUsed)
#
# In a sandbox the payloads are tiny and sub-60s, so the real monitor records a
# near-zero peak and ~1 min -- useless for curve fitting. When a resource model
# is active (SS_MODEL set and mock sbatch has written $id.resource), this shim
# instead writes a synthetic memCPU file whose START is BACK-DATED so cleanUp.sh
# computes the modeled minutes, and whose sample peak is the modeled memory.
#
# When no model is active it EXECs the real monitor unchanged, so functional
# tests and the Tier-2 real-capture fixtures still exercise genuine sampling.
#
# Args (same as the real tool): <flag> <reservedMem> <reservedTime> <defaultMem> <core>
# ---------------------------------------------------------------------------
set -uo pipefail

_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
_real="$SMARTSLURM_ROOT/bin/memCpuMonitor.sh"
res="${MOCK_SLURM_STATE:-/tmp/mockslurm}/jobs/${SLURM_JOBID:-0}.resource"

# No active model (or no cached resource for this job) -> be the real monitor.
if [ -z "${SS_MODEL:-}" ] || [ ! -f "$res" ]; then
  [ -x "$_real" ] && exec "$_real" "$@"
  exit 0
fi

read -r memUsed timeUsed maxRSS _state < "$res"
: "${memUsed:=1}" "${timeUsed:=1}"

logdir="${smartSlurmLogDir:-smartSlurmLog}"
mkdir -p "$logdir"
f="$logdir/job_${SLURM_JOBID}.memCPU.txt"

now=$(date +%s)
# Back-date START so cleanUp.sh's ceil((FINISH-START)/60) == timeUsed, with a
# ~30s margin to absorb the real gap between this write and cleanUp reading it.
offset=$(( timeUsed*60 - 30 )); [ "$offset" -lt 1 ] && offset=1
start=$(( now - offset ))

# line 1: fields "0 0 0 0 0 <START> ..."  (field 6 = START, as the real tool)
printf '0 0 0 0 0 %s jobID:%s time:%s\n' "$start" "${SLURM_JOBID}" "$(date)" > "$f"
# one sample: field 2 = peak mem (M)   ->  "<counter> <mem> <remaining> <saved> <cpu>"
printf '1 %s 0 0 0\n' "$memUsed" >> "$f"
