#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# jr_curate_adversary3.sh -- third-pass adversarial guard for bin/jrCurate.sh
# (SmartSlurm feat/curation-fixes @ f2dc474). Run:
#
#   SMARTSLURM_ROOT=/path/to/smartslurm bash test/jr_curate_adversary3.sh
#
# FINDING C (post-apply gate FALSE FAILURE on a jobID collision):
#   cmd_import's post-apply integrity gate re-reads the store and flags a failure
#   for ANY row whose jobID is in the assignment map but whose pipeline != the
#   chosen name:
#       NR>1 && ($5 in want){ seen[$5]=1; if($1!=want[$5]) wrong++ }   (jrCurate.sh:278)
#   But the relabel only touches rows that are pipeline==legacy. If a NON-legacy
#   row happens to share a jobID with an assigned legacy row (job-id counter reset,
#   multi-cluster pooled store, merged history), that non-legacy row legitimately
#   keeps its OWN name -- yet the gate counts it as "wrong", so cmd_import dies
#   "post-apply check failed ... store left as-is -- inspect it" with rc!=0, even
#   though the relabel was applied EXACTLY correctly.
#
#   Thanks to the round-2 no-restore fix the store is NOT clobbered (no data loss),
#   but the user is told curation FAILED and pointed at the backup; a user who
#   trusts the message and restores the backup reverts a correct curation.
#
#   Required: when every assigned legacy row is correctly relabeled, import must
#   SUCCEED (rc=0) regardless of an unrelated same-jobID non-legacy row, which must
#   keep its own name.
#
# This guard is EXPECTED TO BE RED on f2dc474 (it documents the open finding). It
# goes green once the gate scopes its check to legacy-origin rows (e.g. verify no
# assigned jobID still reads "legacy" AND each assigned jobID appears with its new
# name, instead of demanding EVERY same-jobID row carry the new name).
# ---------------------------------------------------------------------------
set -u
ROOT="${SMARTSLURM_ROOT:-}"
if [ -z "$ROOT" ] || [ ! -r "$ROOT/bin/jrCurate.sh" ]; then
    echo "jr_curate_adversary3: set SMARTSLURM_ROOT to a checkout with bin/jrCurate.sh" >&2
    exit 2
fi
. "$ROOT/bin/jobRecordSchema.sh"
CURATE(){ bash "$ROOT/bin/jrCurate.sh" "$@"; }
PASS=0; FAIL=0
ok(){ PASS=$((PASS+1)); printf '  ok  : %s\n' "$1"; }
no(){ FAIL=$((FAIL+1)); printf '  FAIL: %s\n' "$1"; }
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
V2(){ printf '%s,%s,%s,%s,%s,100,1000,2000,10,20,1,COMPLETED,u,n,%s,2026-01-01,%s\n' "$1" "$2" "$3" "$4" "$5" "$6" "$7"; }

echo "FINDING C: a non-legacy row sharing a jobID with an assigned legacy row must not false-fail import"
d="$TMP/C"; mkdir -p "$d"; jr="$d/jobRecord.txt"
printf '%s\n' "$JR_HEADER_V2" > "$jr"
V2 legacy align hg38 sA 500 1735700500 2.1.align.sA >> "$jr"   # legacy row, jobID 500
V2 rnaseq count hg38 sB 500 1735700501 3.2.count.sB >> "$jr"   # NON-legacy, SAME jobID 500
CURATE export "$d" >/dev/null 2>&1
sed -i 's/^legacyPipe1 | -/legacyPipe1 | realp/' "$d/jobRecord.curate.manifest"
CURATE import "$d/jobRecord.curate.manifest" "$d" >"$d/out" 2>&1; rc=$?

leg=$(awk -F, '$5==500 && $2=="align"{print $1}' "$jr")
rna=$(awk -F, '$5==500 && $2=="count"{print $1}' "$jr")

# Invariant 1 (data safety -- holds on f2dc474 thanks to no-restore): both rows present, correctly labeled.
if [ "$leg" = realp ] && [ "$rna" = rnaseq ]; then
    ok "no data loss: legacy row relabeled (realp) AND non-legacy row kept its name (rnaseq)"
else
    no "relabel wrong / data lost: legacy->'$leg' (want realp), non-legacy->'$rna' (want rnaseq)"
fi

# Invariant 2 (the finding -- RED on f2dc474): a correct relabel must report success.
if [ "$rc" = 0 ]; then
    ok "import reported success for a correct relabel despite the same-jobID non-legacy row"
else
    no "FALSE FAILURE: relabel was correct but import exited $rc and told the user it failed -- msg: $(grep -i 'post-apply' "$d/out" | head -1)"
fi

echo "----------------------------------------------------------------"
printf 'jr_curate_adversary3: PASS=%d FAIL=%d\n' "$PASS" "$FAIL"
[ "$FAIL" -eq 0 ]
