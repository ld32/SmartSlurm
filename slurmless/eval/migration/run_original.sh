#!/bin/bash
# run_original.sh <variantDir> : establish the REF (reference behavior) for a variant
# by running its original as intended, in the slurmless sandbox. Writes
# <variantDir>/ref/quote.txt. SLURM variants (submit.sh) run in wf-mode; loose
# variants run reference_driver.sh (the intended correct wiring).
set -uo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
v="${1:?usage: run_original.sh <variantDir>}"
vdir="$(cd "$v" && pwd)"
root="$(cd "$here/../.." && pwd)"                 # eval/migration -> repo root
: "${SMARTSLURM_ROOT:?set SMARTSLURM_ROOT to a SmartSlurm checkout}"
# shellcheck disable=SC1090
source "$root/test/lib/harness.sh"
ss_activate >/dev/null 2>&1
ws="${SS_SCRATCH:-/tmp/ss-scratch}/eval_orig_$(basename "$vdir")"
ss_new_workspace "$ws" >/dev/null
cp "$vdir"/*.sh "$vdir"/*.txt . 2>/dev/null; chmod +x ./*.sh 2>/dev/null || true
if   [ -f submit.sh ];           then ss_wf_activate; bash ./submit.sh >/dev/null 2>&1; ss_wf_drain
elif [ -f reference_driver.sh ]; then bash ./reference_driver.sh >/dev/null 2>&1
else echo "run_original: no submit.sh or reference_driver.sh in $vdir" >&2; exit 2; fi
mkdir -p "$vdir/ref"; cp quote.txt "$vdir/ref/quote.txt" 2>/dev/null || true
echo "REF[$(basename "$vdir")]: $(cat "$vdir/ref/quote.txt" 2>/dev/null)"
