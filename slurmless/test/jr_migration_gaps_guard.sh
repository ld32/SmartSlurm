#!/usr/bin/env bash
# Regression guard for jobRecord v1->v2 migration paths the existing suite
# (jr_conversion_stress.sh + the migrate_v1 block) does NOT exercise. Each check
# was proven to have TEETH: it fails when the matching code path is neutered, and
# passes on the current fix/jobrecord-v2-recognition tree. Run it under BOTH awks:
#
#   SMARTSLURM_ROOT=/path/to/smartslurm test/jr_migration_gaps_guard.sh
#   d=$(mktemp -d); ln -sf "$(command -v mawk)" "$d/awk"; \
#     PATH="$d:$PATH" SMARTSLURM_ROOT=/path/to/smartslurm test/jr_migration_gaps_guard.sh
#
# Covered gaps (adversary pass, 2026-11):
#   1. Case-2 dup/collision header diverted to .unrecognized during MIGRATION
#      (jrLint reporting was guarded; the migrate-time diversion+warning was not).
#   2. _jrLock stale-lock reclaim by jrMigrateStore (crashed holder must not wedge).
#   3. jrMigrateStore skips (does not corrupt/hang) when a FRESH lock is held.
#   4. A leftover .v2.PID/.rw.PID crash temp does not break the next migration.
#   5. Estimate-equivalence across migration: the fit inputs the estimator reads
#      (inputSize/memUsed/timeUsed per step/ref) are byte-identical v1-on-the-fly
#      vs after persisting the migration (no silent estimation drift).
set -u
ROOT="${SMARTSLURM_ROOT:-}"
if [ -z "$ROOT" ] || [ ! -r "$ROOT/bin/jobRecordSchema.sh" ]; then
    echo "jr_migration_gaps_guard: set SMARTSLURM_ROOT to a SmartSlurm checkout" >&2; exit 2
fi
. "$ROOT/bin/jobRecordSchema.sh"
PASS=0; FAIL=0; FAILS=()
ok(){ PASS=$((PASS+1)); }
no(){ FAIL=$((FAIL+1)); FAILS+=("$1"); printf '  FAIL: %s\n' "$1"; }
chk(){ if eval "$1"; then ok; else no "$2 [cond: $1]"; fi; }
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT

v1r21(){ printf '%s,%s,0,0,4096,120,%s,95,%s,u,,%s,%s,1.0.%s.%s,%s,0,0,0,0,2026-01-01,%s\n' \
  "$1" "${2:-5000}" "${3:-3800}" "${4:-COMPLETED}" "${5:-bwa}" "${6:-hg38}" "${5:-bwa}" "$1" "${7:-2}" "$((1735700000+$1))"; }
H_V1_21="1jobID,2inputSize,3memDefault,4timeDefault,5memAllocated,6timeAllocated,7memUsed,8timeUsed,9status,10user,11saccMem,12program,13ref,14flag,15cores,16extraMem,17extraTime,18ratioM,19ratioT,20date,21epoch"

echo "== gap 1: duplicate/collision header is diverted to .unrecognized DURING migration =="
# Realistic trigger: two legacy stores concatenated -> a second "1jobID" header
# line beyond the real one. It is byte-identical to the header, excluded from data,
# and invisible to the recog/bad accounting -- the migration must surface it loudly
# and preserve it in a .unrecognized file (not drop it, not count it as data).
f="$TMP/dup.txt"
{ echo "$H_V1_21"; v1r21 100; echo "$H_V1_21"; v1r21 101; } > "$f"
err=$(jrMigrateStore "$f" 2>&1 1>/dev/null)
chk "[ \"\$(jrSchemaVersion '$f')\" = 2 ]"               "gap1: store migrates to v2"
chk "echo \"\$err\" | grep -qi 'identical to the jobRecord header'" "gap1: migration WARNS about the extra header-shaped line"
uf=$(ls "$f".unrecognized.*.txt 2>/dev/null | head -1)
chk "[ -n '$uf' ] && grep -q '1jobID' '$uf'"            "gap1: the extra header line is saved to a .unrecognized file"
# the two real data rows still survive as v2
got=$(awk -F, 'NR>1 && NF==17 && $1 !~ /^[0-9]+$/{print $5}' "$f" | sort -n | xargs)
chk "[ '$got' = '100 101' ]"                            "gap1: both real data rows survive (want '100 101' got '$got')"

echo "== gap 2: jrMigrateStore reclaims a STALE lock (crashed holder must not wedge) =="
f="$TMP/stale.txt"; v1r21 200 > "$f"
mkdir "$f.lock"; touch -d '120 seconds ago' "$f.lock"
JR_LOCK_STALE=30 JR_LOCK_TRIES=5 jrMigrateStore "$f" 2>/dev/null
chk "[ \"\$(jrSchemaVersion '$f')\" = 2 ]"              "gap2: stale lock reclaimed -> store migrated to v2"
chk "[ ! -d '$f.lock' ]"                                "gap2: lock released after reclaim"

echo "== gap 3: jrMigrateStore SKIPS (no corrupt, no hang) when a FRESH lock is held =="
f="$TMP/held.txt"; v1r21 300 > "$f"; cp "$f" "$f.orig"
mkdir "$f.lock"   # fresh (not stale)
JR_LOCK_STALE=99999 JR_LOCK_TRIES=3 jrMigrateStore "$f" 2>/dev/null
chk "diff -q '$f' '$f.orig' >/dev/null"                 "gap3: held-lock migrate left the store byte-identical"
chk "[ \"\$(jrSchemaVersion '$f')\" = 1 ]"              "gap3: store still v1 (migration correctly skipped)"
rmdir "$f.lock" 2>/dev/null

echo "== gap 4: a leftover crash temp (.v2.PID / .rw.PID) does not break the next migration =="
f="$TMP/crash.txt"; v1r21 400 > "$f"
: > "$f.v2.99999"; : > "$f.rw.88888"           # litter from a killed writer
jrMigrateStore "$f" 2>/dev/null
got=$(awk -F, 'NR>1{print $5}' "$f" | xargs)
chk "[ \"\$(jrSchemaVersion '$f')\" = 2 ] && [ '$got' = 400 ]" "gap4: migration succeeds despite leftover temp litter"

echo "== gap 5: estimate-equivalence across migration (fit inputs unchanged) =="
# The estimator reads (inputSize, memUsed, timeUsed) per step/ref via jrEmitV2.
# A v1 store read on the fly MUST yield the identical fit inputs as the same store
# after it is persisted to v2 -- otherwise migration silently shifts estimates.
f="$TMP/est.txt"
for i in 1 2 3 4 5; do v1r21 $((500+i)) $((1000*i)) $((300*i)) COMPLETED bwa hg38 2 >> "$f"; done
fitpull(){ jrEmitV2 "$1" | awk -F, 'NR==1{for(i=1;i<=NF;i++)c[$i]=i;next}
  $(c["step"])=="bwa" && $(c["ref"])=="hg38" && $(c["status"])=="COMPLETED" && $(c["inputSize"])~/^[0-9]+$/ \
  {print $(c["inputSize"]),$(c["memUsed"]),$(c["timeUsed"])}' | sort -n; }
pre=$(fitpull "$f"); jrMigrateStore "$f" 2>/dev/null; post=$(fitpull "$f")
chk "[ -n \"\$pre\" ]"                                   "gap5: pre-migration fit inputs are non-empty"
chk "[ \"\$pre\" = \"\$post\" ]"                        "gap5: fit inputs identical pre/post migration (no estimate drift)"

echo
echo "================  RESULT  ================"
echo "PASS=$PASS  FAIL=$FAIL"
[ "$FAIL" -eq 0 ] && echo "ALL GREEN" || { printf '%s\n' "${FAILS[@]}"; exit 1; }
