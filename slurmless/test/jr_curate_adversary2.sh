#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# jr_curate_adversary2.sh -- second-pass adversarial guards for bin/jrCurate.sh
# (SmartSlurm feat/curation-fixes). Run:
#
#   SMARTSLURM_ROOT=/path/to/smartslurm bash test/jr_curate_adversary2.sh
#
# Each guard asserts a CONTRACT invariant, not an exit code. exit 0 = all hold.
#
# GUARD A (SILENT DATA LOSS -- headline): cmd_import takes its backup OUTSIDE the
#   store lock, then both failure branches blindly `cp -p "$bak" "$jr"`. A record
#   that a concurrent writer (cleanUp's jrAppendLine, which appends even without
#   the lock) commits AFTER the backup but BEFORE the restore is reverted away --
#   a COMPLETED jobRecord vanishes with no error, and the tool even prints
#   "store unchanged". Required: a failed/contended import must NEVER delete a
#   record that was not in the backup.
#
# GUARD B (FALSE REJECTION): the manifest is "|"-separated and import rejects any
#   line with NF>5 ("a value contains the | separator"). But the read-only
#   steps/refs hint columns are copied verbatim from the store, and jrRow/jrEmitV2
#   sanitize ", \t\r" but NOT "|". A legacy row whose ref or program legitimately
#   contains "|" produces a manifest line with NF>5, so import rejects the WHOLE
#   manifest -- including the line the user correctly edited -- and the store can
#   never be curated. Required: a pipe that lives only in a read-only hint column
#   (never in the edited newName) must not block the import.
# ---------------------------------------------------------------------------
set -u
ROOT="${SMARTSLURM_ROOT:-}"
if [ -z "$ROOT" ] || [ ! -r "$ROOT/bin/jobRecordSchema.sh" ] || [ ! -r "$ROOT/bin/jrCurate.sh" ]; then
    echo "jr_curate_adversary2: set SMARTSLURM_ROOT to a checkout with bin/jrCurate.sh + bin/jobRecordSchema.sh" >&2
    exit 2
fi
. "$ROOT/bin/jobRecordSchema.sh"
CURATE(){ bash "$ROOT/bin/jrCurate.sh" "$@"; }
PASS=0; FAIL=0
ok(){ PASS=$((PASS+1)); printf '  ok  : %s\n' "$1"; }
no(){ FAIL=$((FAIL+1)); printf '  FAIL: %s\n' "$1"; }
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT

# ======================================================================
echo "GUARD A: a contended/failed import must not clobber a concurrently-appended record"
# ======================================================================
d="$TMP/A"; mkdir -p "$d"; jr="$d/jobRecord.txt"
printf '%s\n' "$JR_HEADER_V2" > "$jr"
printf 'legacy,align,hg38,sA,101,100,1000,2000,10,20,1,COMPLETED,u,n,1735700101,2026-01-01,2.1.align.sA\n' >> "$jr"
CURATE export "$d" >/dev/null 2>&1
sed -i 's/^legacyPipe1 | -/legacyPipe1 | curated1/' "$d/jobRecord.curate.manifest"
R='rnaseq,count,hg38,sB,777,200,2000,4000,5,10,2,COMPLETED,u,n,1735700777,2026-01-02,3.2.count.sB'
(   # concurrent writer: hold the lock, wait for import's backup, append R, hold past timeout
  mkdir "$jr.lock" 2>/dev/null
  for i in $(seq 1 300); do ls "$jr".back.*.txt >/dev/null 2>&1 && break; sleep 0.01; done
  printf '%s\n' "$R" >> "$jr"
  sleep 2.5
  rmdir "$jr.lock" 2>/dev/null
) &
holder=$!
( JR_LOCK_TRIES=30 JR_LOCK_STALE=30 CURATE import "$d/jobRecord.curate.manifest" "$d" ) >"$TMP/impA.out" 2>&1
wait "$holder" 2>/dev/null
if awk -F, 'NR>1 && $5==777{f=1} END{exit !f}' "$jr"; then
    ok "concurrent COMPLETED record (jobID 777) survived a contended import"
else
    no "SILENT DATA LOSS: concurrent record 777 clobbered by restore-from-backup; import said: $(cat "$TMP/impA.out")"
fi

# ======================================================================
echo "GUARD B: a pipe in a read-only hint column must not false-reject an unedited line"
# ======================================================================
d="$TMP/B"; mkdir -p "$d"; jr="$d/jobRecord.txt"
printf '%s\n' "$JR_HEADER_V2" > "$jr"
# legacy rows whose ref legitimately contains '|' (survives jrRow/jrEmitV2 sanitize)
printf 'legacy,align,hg38|alt,sA,201,100,1000,2000,10,20,1,COMPLETED,u,n,1735700201,2026-01-01,2.1.align.sA\n' >> "$jr"
printf 'legacy,trim,hg38|alt,sA,202,100,1000,2000,10,20,1,COMPLETED,u,n,1735700202,2026-01-01,1.0.trim.sA\n' >> "$jr"
CURATE export "$d" >/dev/null 2>&1
sed -i 's/^legacyPipe1 | -/legacyPipe1 | realpipe/' "$d/jobRecord.curate.manifest"
CURATE import "$d/jobRecord.curate.manifest" "$d" >"$TMP/impB.out" 2>&1
r201=$(awk -F, '$5==201{print $1}' "$jr")
if [ "$r201" = realpipe ]; then
    ok "unedited line with a '|' in the ref hint still imported (201 -> realpipe)"
else
    no "FALSE REJECTION: a '|' in the read-only ref column blocked the import (201 still '$r201'); import said: $(cat "$TMP/impB.out")"
fi

echo "----------------------------------------------------------------"
printf 'jr_curate_adversary2: PASS=%d FAIL=%d\n' "$PASS" "$FAIL"
[ "$FAIL" -eq 0 ]
