#!/bin/bash
# ---------------------------------------------------------------------------
# check_sacct_format.sh  --  detect drift between the mock sacct and a real one.
#
# The mock sacct (test/mockSlurm/sacct) is the most Slurm-version-sensitive file
# in the harness: cleanUp.sh parses sacct's tabular output by COLUMN, so if a
# Slurm upgrade changes sacct's column set, names, or ordering for the format
# cleanUp.sh requests, the mock silently diverges from the cluster.
#
# Run this:
#   - in the harness repo (validates the mock alone), or
#   - ON A SLURM CLUSTER with a completed job id, to DIFF the mock against the
#     real sacct for the exact --format cleanUp.sh uses:
#         test/mockSlurm/check_sacct_format.sh <realJobId>
#
# Exit 0 = match (or mock-only sanity ok); exit 1 = drift detected. CI-friendly.
# ---------------------------------------------------------------------------
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MOCK_SACCT="$HERE/sacct"

# The exact --format cleanUp.sh requests (keep in sync with bin/cleanUp.sh). The
# %width suffixes do not change the column SET, only its print width.
FMT='JobID,Submit,Start,End,MaxRSS,State,NodeList%30,Partition,ReqTRES%30,TotalCPU,Elapsed%14,Timelimit%14,ExitCode'

# canonical column names (strip %width), in order
canon_cols() { echo "$1" | tr ',' '\n' | sed 's/%.*//'; }
EXPECTED="$(canon_cols "$FMT" | paste -sd, -)"

echo "cleanUp.sh --format : $FMT"
echo "expected columns    : $EXPECTED"
echo

# --- mock side: the mock's header for this format must match the column set ----
mock_hdr="$(MOCK_SLURM_STATE="$(mktemp -d)" "$MOCK_SACCT" --format="$FMT" -j 1 2>/dev/null | head -1)"
mock_cols="$(echo "$mock_hdr" | tr -s ' ' '\n' | sed '/^$/d' | paste -sd, -)"
echo "mock sacct columns  : $mock_cols"
rc=0
if [ "$mock_cols" != "$EXPECTED" ]; then
  echo "DRIFT: mock sacct header does not match the requested format." >&2
  echo "       (did test/mockSlurm/sacct's colval() fall out of sync with FMT?)" >&2
  rc=1
else
  echo "ok - mock sacct emits exactly the requested columns"
fi

# --- real side (optional): diff against the cluster's sacct -------------------
realId="${1:-}"
if [ -n "$realId" ] && command -v sacct >/dev/null 2>&1 && [[ "$(readlink -f "$(command -v sacct)")" != "$MOCK_SACCT" ]]; then
  echo
  echo "real sacct version  : $(sacct --version 2>/dev/null || echo '?')"
  real_hdr="$(sacct --format="$FMT" -j "$realId" 2>/dev/null | head -1)"
  real_cols="$(echo "$real_hdr" | tr -s ' ' '\n' | sed '/^$/d' | paste -sd, -)"
  echo "real sacct columns  : $real_cols"
  if [ "$real_cols" != "$EXPECTED" ]; then
    echo "DRIFT: real sacct columns differ from the format cleanUp.sh requests." >&2
    echo "       Update test/mockSlurm/sacct (colval/row emission) and the pin in" >&2
    echo "       test/mockSlurm/README.md to the new Slurm version." >&2
    rc=1
  else
    echo "ok - real sacct matches the expected column set for this format"
  fi
else
  echo
  echo "(no real sacct diff: pass a completed job id on a Slurm cluster to enable it)"
fi

exit $rc
