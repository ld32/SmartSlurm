#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# jr_pipeline_identity_guard.sh -- adversarial regression guards for the
# runAsPipeline pipeline-name identity/validation feature and its interaction
# with the up-front jobRecord v1->v2 migration (the feat/pipeline-and-curation
# merge). Run:
#
#   SMARTSLURM_ROOT=/path/to/smartslurm bash test/jr_pipeline_identity_guard.sh
#
# Validate under BOTH awks:
#   SD=$(mktemp -d); ln -sf "$(command -v mawk)" "$SD/awk"
#   PATH="$SD:$PATH" SMARTSLURM_ROOT=... bash test/jr_pipeline_identity_guard.sh
#
# SECTION A (must be GREEN -- these had NO harness coverage before this guard;
#   a mutation that removes the all-digit/reserved/header-collision rejection in
#   runAsPipeline left the whole jR suite green, yet such a name makes jrEmitV2
#   silently drop EVERY row of that pipeline. Exit code reflects Sections A and B.):
#     A1 an invalid --pipeline (all-digit / reserved / header-collision) is
#        rejected with rc!=0 AND the store is left UNMIGRATED (validation must run
#        BEFORE the up-front migration -- a rejected run must not mutate anything).
#     A2 an invalid #@pipeline DIRECTIVE is likewise rejected before any mutation.
#     A3 a valid --pipeline NAME flows into the jobRecord pipeline column; migrated
#        legacy rows stay pipeline=legacy (cross-pipeline isolation through migration).
#     A4 a valid #@pipeline DIRECTIVE flows into the jobRecord pipeline column.
#
# SECTION B (hard asserts -- dryrun / lint must be read-only on the store):
#     B1 `runAsPipeline --mode dryrun` MUST NOT mutate the store. The per-step
#        `ssbatch ... dryrun` path used to call jrEnsureHeader (ssbatch:91) BEFORE
#        checking the trailing `dryrun` arg, so a dry run migrated a v1 store
#        v1->v2 and wrote a .back file -- contradicting runAsPipeline's own
#        gate+comment ("a --lint / --mode dryrun must not mutate the store").
#        Fixed on feat/pipeline-and-curation: ssbatch skips jrEnsureHeader when
#        its last arg is `dryrun` (real runs still ensure the header, and
#        runAsPipeline migrates up front). This assert pins that fix across the
#        full matrix of real store shapes (absent, empty, v1 21-col/20-col, mixed
#        width, headerless-v2, canonical-v2, mixed v1+v2, hybrid) -- each must be
#        byte-for-byte untouched by a dry run, and an absent store must stay absent.
#     B2 `--lint` is the control -- always read-only (it exits before the driver).
# ---------------------------------------------------------------------------
set -u
ROOT="${SMARTSLURM_ROOT:-}"
if [ -z "$ROOT" ] || [ ! -x "$ROOT/bin/runAsPipeline" ] || [ ! -r "$ROOT/bin/jobRecordSchema.sh" ]; then
    echo "set SMARTSLURM_ROOT to a checkout with bin/runAsPipeline + bin/jobRecordSchema.sh" >&2
    exit 2
fi
HERE="$(cd "$(dirname "$0")/.." && pwd)"
sch="$ROOT/bin/jobRecordSchema.sh"
jrver() { ( . "$sch"; jrSchemaVersion "$1" ); }

SB="$(mktemp -d)"; trap 'rm -rf "$SB"' EXIT
export SS_SCRATCH="$SB/scratch"; mkdir -p "$SS_SCRATCH"
. "$HERE/test/lib/harness.sh"; ss_activate >/dev/null 2>&1

PASS=0; FAIL=0
ok()   { PASS=$((PASS+1)); echo "  ok  : $1"; }
bad()  { FAIL=$((FAIL+1)); echo "  FAIL: $1"; }

seed_v1() {   # one 21-col legacy COMPLETED row for jobID $1 (default 1001)
    local n="${1:-1001}"
    printf '%s,%s,0,0,4096,120,%s,95,COMPLETED,tester,,bwa,none,1.0.bwa.%s,1,0,0,0,0,2026-01-01,%s\n' \
        "$n" 5000 3800 "$n" "$((1735700000+n))"
}
setup_ws() {  # clean workspace with the clean_2step fixture; sets $jr, leaves the
              # store ABSENT (caller places the shape under test, or seeds v1)
    ss_new_workspace "$SS_SCRATCH/wk" >/dev/null 2>&1
    cp "$HERE"/test/fixtures/clean_2step/* .; chmod +x ./*.sh
    jr="$smartSlurmJobRecordDir/jobRecord.txt"
    rm -f "$jr" "$jr".back.*.txt "$jr".unrecognized.*.txt
}
fresh_ws() { setup_ws; seed_v1 1001 > "$jr"; }   # + a fresh 21-col v1 store
rap() { bash "$ROOT/bin/runAsPipeline" "$@" -q >"$SS_SCRATCH/o.log" 2>&1; }

# Canonical store-shape constructors, byte-matched to run_tests.sh migrate_v1
# (the real layouts seen in the wild). v1row21/20 = legacy; v2row = canonical v2.
v1row21() { printf '%s,%s,0,0,4096,120,%s,95,%s,tester,,%s,%s,%s,1,0,0,0,0,2026-01-01,%s\n' \
              "$1" "$2" "$3" "${6:-COMPLETED}" "$4" "${5:-none}" "1.0.$4.$1" "$(( 1735700000 + $1 ))"; }
v1row20() { printf '%s,%s,0,0,4096,120,%s,95,%s,tester,,%s,%s,%s,1,0,0,0,0,2026-01-01\n' \
              "$1" "$2" "$3" "${6:-COMPLETED}" "$4" "${5:-none}" "1.0.$4.$1"; }
v2row()   { local id=$1 st=$2 rf=$3 sz=$4 mu=$5 sm=${6:-}
    printf 'rnaseq,%s,%s,%s,%s,%s,%s,%s,95,120,1,COMPLETED,tester,n,%s,2026-01-05,1.0.%s.%s\n' \
      "$st" "$rf" "$sm" "$id" "$sz" "$mu" "$mu" "$(( 1735900000 + id ))" "$st" "$id"; }
HDR2='pipeline,step,ref,sample,jobID,inputSize,memUsed,memAllocated,timeUsed,timeAllocated,cores,status,user,node,epoch,isoDate,flag'

# Build store shape $1 into $jr (store already removed by setup_ws). "absent"
# leaves no file at all.
mk_shape() {
    case "$1" in
      absent)        ;;  # no file
      empty)         : > "$jr" ;;
      v1_21col)      v1row21 1001 5000 3800 bwa hg38 > "$jr" ;;
      v1_20col)      v1row20 1002 6000 7900 bowtie2 mm10 > "$jr" ;;
      v1_mixedwidth) { v1row21 1001 5000 3800 bwa hg38; v1row20 1002 6000 7900 bowtie2 mm10; } > "$jr" ;;
      headerless_v2) { v2row 1001 bwa hg38 5000 3800 sampA; v2row 1002 bwa hg38 6000 4200 sampB; } > "$jr" ;;
      canonical_v2)  { echo "$HDR2"; v2row 1001 bwa hg38 5000 3800 sampA; v2row 1002 bwa hg38 6000 4200 sampB; } > "$jr" ;;
      mixed_v1v2)    { v1row21 2001 5000 3800 bwa hg38; v1row20 2002 6000 7900 bwa hg38
                       printf 'rnaseq,align,hg38,sampC,2003,7000,9000,10000,300,600,4,COMPLETED,tester,,1735800000,2026-01-03,4.align.sampC\n'; } > "$jr" ;;
      hybrid_v2first){ echo "$HDR2"; v2row 801 bwa hg38 5000 3800 sampA; v1row21 802 6000 7900 bwa hg38; } > "$jr" ;;
      *)             echo "mk_shape: unknown shape '$1'" >&2; return 2 ;;
    esac
}

# Dry-run read-only assertion for the shape currently in $jr: a `runAsPipeline
# --mode dryrun` must leave the store byte-for-byte unchanged -- no migration, no
# .back, no .unrecognized, and an absent store must STAY absent (pre-fix,
# jrEnsureHeader created a header file here). Asserts only on store state, not on
# the dry run's exit code (a dry run may legitimately warn while mutating nothing).
check_readonly() {
    local label="$1" pre post prebk postbk postun ex_pre ex_post
    if [ -e "$jr" ]; then ex_pre=1; pre=$(md5sum "$jr" | cut -d' ' -f1); else ex_pre=0; pre=ABSENT; fi
    prebk=$(ls "$jr".back.*.txt 2>/dev/null | wc -l)
    rap --script "./pipeline.sh 777" --tmp noTmp --mode dryrun \
        --sbatch-options "sbatch -p short -c 1 --mem 4G -t 50:0"
    if [ -e "$jr" ]; then ex_post=1; post=$(md5sum "$jr" | cut -d' ' -f1); else ex_post=0; post=ABSENT; fi
    postbk=$(ls "$jr".back.*.txt 2>/dev/null | wc -l)
    postun=$(ls "$jr".unrecognized.*.txt 2>/dev/null | wc -l)
    if [ "$ex_pre" = "$ex_post" ] && [ "$pre" = "$post" ] && [ "$postbk" -eq "$prebk" ] && [ "$postun" -eq 0 ]; then
        ok "B1[$label] dryrun read-only (md5 $pre, backups $postbk, unrecognized 0)"
    else
        bad "B1[$label] dryrun MUTATED store: exists $ex_pre->$ex_post, md5 $pre->$post, backups $prebk->$postbk, unrecognized=$postun -- ssbatch dryrun guard regressed (must skip jrEnsureHeader when last arg is 'dryrun')"
    fi
}

echo "SECTION A -- pipeline-name validation + identity flow (must be green)"

# --- A1: invalid --pipeline rejected before migration ------------------------
for bad_name in 123 legacy unknown pipeline 1jobID; do
    fresh_ws
    rap --script "./pipeline.sh 777" --tmp noTmp --pipeline "$bad_name" \
        --sbatch-options "sbatch -p short -c 1 --mem 4G -t 50:0"
    rc=$?; after=$(jrver "$jr"); bk=$(ls "$jr".back.*.txt 2>/dev/null | wc -l)
    if [ "$rc" -ne 0 ] && [ "$after" = 1 ] && [ "$bk" -eq 0 ]; then
        ok "A1 --pipeline '$bad_name' rejected (rc=$rc) and store left UNMIGRATED (v1, no backup)"
    else
        bad "A1 --pipeline '$bad_name': rc=$rc store-ver=$after backups=$bk (expected rc!=0, ver=1, 0 backups)"
    fi
done

# --- A2: invalid #@pipeline directive rejected before migration --------------
for bad_name in 999 legacy; do
    fresh_ws
    sed -i "2i #@pipeline $bad_name" pipeline.sh
    rap --script "./pipeline.sh 777" --tmp noTmp \
        --sbatch-options "sbatch -p short -c 1 --mem 4G -t 50:0"
    rc=$?; after=$(jrver "$jr"); bk=$(ls "$jr".back.*.txt 2>/dev/null | wc -l)
    if [ "$rc" -ne 0 ] && [ "$after" = 1 ] && [ "$bk" -eq 0 ]; then
        ok "A2 #@pipeline '$bad_name' rejected (rc=$rc) and store left UNMIGRATED"
    else
        bad "A2 #@pipeline '$bad_name': rc=$rc store-ver=$after backups=$bk (expected rc!=0, ver=1)"
    fi
done

# --- A3: valid --pipeline flows into the pipeline column; legacy preserved ----
fresh_ws
rap --script "./pipeline.sh 777" --tmp noTmp --pipeline myrna \
    --sbatch-options "sbatch -p short -c 1 --mem 4G -t 50:0"
rc=$?
named=$(awk -F, '$1=="myrna"{c++} END{print c+0}' "$jr")
leg=$(awk -F, '$1=="legacy" && $5==1001{c++} END{print c+0}' "$jr")
relabeled=$(awk -F, '$5==1001 && $1=="myrna"{c++} END{print c+0}' "$jr")
if [ "$rc" -eq 0 ] && [ "$named" -ge 1 ] && [ "$leg" -eq 1 ] && [ "$relabeled" -eq 0 ]; then
    ok "A3 --pipeline myrna: new rows tagged myrna ($named), legacy row 1001 preserved as legacy, not relabeled"
else
    bad "A3 --pipeline myrna: rc=$rc named=$named legacyKept=$leg relabeled=$relabeled (want rc=0,named>=1,legacyKept=1,relabeled=0)"
fi

# --- A4: valid #@pipeline directive flows into the pipeline column ------------
fresh_ws
sed -i '2i #@pipeline rnaDir' pipeline.sh
rap --script "./pipeline.sh 777" --tmp noTmp \
    --sbatch-options "sbatch -p short -c 1 --mem 4G -t 50:0"
rc=$?
dnamed=$(awk -F, '$1=="rnaDir"{c++} END{print c+0}' "$jr")
dleg=$(awk -F, '$1=="legacy" && $5==1001{c++} END{print c+0}' "$jr")
if [ "$rc" -eq 0 ] && [ "$dnamed" -ge 1 ] && [ "$dleg" -eq 1 ]; then
    ok "A4 #@pipeline rnaDir: new rows tagged rnaDir ($dnamed), legacy row preserved"
else
    bad "A4 #@pipeline rnaDir: rc=$rc named=$dnamed legacyKept=$dleg (want rc=0,named>=1,legacyKept=1)"
fi

echo
echo "SECTION B -- dryrun / lint are read-only on the store (hard asserts)"
# --- B1: --mode dryrun must not mutate the store, across every store shape ----
# Broadened from one shape (21-col v1) to the full matrix of real-world layouts,
# so a regression on headerless-v2 / 20-col / mixed / hybrid / absent / empty
# cannot slip past a single-shape assert. Each shape must be byte-for-byte
# untouched by a dry run (pre-fix, a dry run migrated v1/mixed/hybrid shapes and
# even CREATED a store file where none existed).
for shape in absent empty v1_21col v1_20col v1_mixedwidth headerless_v2 canonical_v2 mixed_v1v2 hybrid_v2first; do
    setup_ws
    mk_shape "$shape"
    check_readonly "$shape"
done
# --- B2: --lint must be safe (control) ---------------------------------------
fresh_ws
before_sum=$(md5sum "$jr" | cut -d' ' -f1)
rap --script "./pipeline.sh 777" --tmp noTmp --lint \
    --sbatch-options "sbatch -p short -c 1 --mem 4G -t 50:0"
after=$(jrver "$jr"); bk=$(ls "$jr".back.*.txt 2>/dev/null | wc -l)
after_sum=$(md5sum "$jr" | cut -d' ' -f1)
if [ "$after" = 1 ] && [ "$bk" -eq 0 ] && [ "$before_sum" = "$after_sum" ]; then
    ok "B2 --lint left the store untouched (control: lint exits before the driver runs)"
else
    bad "B2 --lint MUTATED the store (ver 1->$after, backups=$bk) -- lint regression"
fi

echo
echo "----------------------------------------------------------------"
echo "jr_pipeline_identity_guard: Section-A/B2 PASS=$PASS FAIL=$FAIL"
[ "$FAIL" -eq 0 ] && echo "SECTION A GREEN" || echo "SECTION A RED"
exit $(( FAIL > 0 ? 1 : 0 ))
