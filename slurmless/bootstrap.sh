#!/bin/bash
# ---------------------------------------------------------------------------
# bootstrap.sh  --  one command to stand slurmless up and prove it works.
#
#   ./bootstrap.sh [BLOCK]     # default BLOCK = clean_2step (the baseline)
#
# It finds the SmartSlurm checkout to test (the parent checkout, $SMARTSLURM_ROOT, or
# a fresh clone of upstream master), activates the Slurm-free environment, and
# runs one test block so you get a green/red answer immediately.
#
# Env:
#   SMARTSLURM_ROOT   use this SmartSlurm checkout (skips discovery/clone)
#   SS_NO_SH_REPOINT  do not repoint /bin/sh -> bash (see lib/harness.sh)
#   SS_NO_INSTALL     do not apt-get install gnuplot
# ---------------------------------------------------------------------------
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BLOCK="${1:-clean_2step}"
UPSTREAM="https://github.com/ld32/smartSlurm"

# 1) Resolve a SmartSlurm checkout.
if [ -z "${SMARTSLURM_ROOT:-}" ]; then
  d="$HERE"
  while [ "$d" != "/" ]; do
    [ -e "$d/bin/runAsPipeline" ] && { SMARTSLURM_ROOT="$d"; break; }
    d="$(dirname "$d")"
  done
fi
if [ -z "${SMARTSLURM_ROOT:-}" ]; then
  echo "[bootstrap] no SmartSlurm checkout found above $HERE and SMARTSLURM_ROOT is unset."
  clone_dir="$(dirname "$HERE")/SmartSlurm"
  if [ -e "$clone_dir/bin/runAsPipeline" ]; then
    SMARTSLURM_ROOT="$clone_dir"
  else
    echo "[bootstrap] cloning upstream master into $clone_dir ..."
    git clone --depth 1 "$UPSTREAM" "$clone_dir" || {
      echo "[bootstrap] clone failed. Set SMARTSLURM_ROOT to a SmartSlurm checkout and retry." >&2
      exit 1; }
    SMARTSLURM_ROOT="$clone_dir"
  fi
fi
export SMARTSLURM_ROOT
echo "[bootstrap] SMARTSLURM_ROOT=$SMARTSLURM_ROOT"

# 2) Run one block (run_tests.sh calls ss_activate for us).
echo "[bootstrap] running block: $BLOCK"
"$HERE/test/run_tests.sh" "$BLOCK"
rc=$?
if [ "$rc" -eq 0 ]; then
  echo "[bootstrap] OK — slurmless is working. Next: ./test/run_tests.sh --list"
else
  echo "[bootstrap] block '$BLOCK' did not pass (exit $rc). See the report above." >&2
fi
exit $rc
