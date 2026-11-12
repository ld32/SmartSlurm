#!/usr/bin/env bash
# Stress / property test for the jobRecord v1->v2 conversion + prune utilities.
# Goal: never again get hybrid (mixed v1/v2) stores or self-destroying prunes.
#
#   SMARTSLURM_ROOT=/path/to/smartslurm test/jr_conversion_stress.sh [SEED]
#
# Exercises jrSchemaVersion / jrMigrateStore / jrEmitV2 / jrRepair and the REAL
# cleanUp prune branch-selection, over (a) a hand-built adversarial matrix and
# (b) a deterministic fuzzer. Checks strong invariants, not just row counts.
# Exit 0 = ALL GREEN; exit 1 = at least one invariant violated (details printed).
# The migrate_v1 harness block invokes this with fixed seeds; run it directly
# with more seeds to fuzz harder. Point SMARTSLURM_ROOT at the bin/ under test.
set -u
ROOT="${SMARTSLURM_ROOT:-}"
if [ -z "$ROOT" ] || [ ! -r "$ROOT/bin/jobRecordSchema.sh" ]; then
    echo "jr_conversion_stress: set SMARTSLURM_ROOT to a SmartSlurm checkout (needs bin/jobRecordSchema.sh)" >&2
    exit 2
fi
SCH="$ROOT/bin/jobRecordSchema.sh"
. "$SCH"
SEED="${1:-12345}"; RANDOM=$SEED
PASS=0; FAIL=0; FAILS=()
ok(){ PASS=$((PASS+1)); }
no(){ FAIL=$((FAIL+1)); FAILS+=("$1"); printf '  FAIL: %s\n' "$1"; }
chk(){ if eval "$1"; then ok; else no "$2 [cond: $1]"; fi; }

TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT

# ---- row/header generators -------------------------------------------------
# v1 21-col: jobID,inputSize,memDef,minDef,allocMem,allocTime,memUsed,timeUsed,
#            status,user,saccMem,program,ref,flag,core,extraMem,extraTime,ratioM,ratioT,date,epoch
v1r21(){ local id=$1 sz=${2:-5000} m=${3:-3800} st=${4:-COMPLETED} pg=${5:-bwa} rf=${6:-hg38} fl=${7:-}
  : "${fl:=1.0.$pg.$id}"
  printf '%s,%s,0,0,4096,120,%s,95,%s,u,,%s,%s,%s,1,0,0,0,0,2026-01-0%s,%s\n' \
  "$id" "$sz" "$m" "$st" "$pg" "$rf" "$fl" "$((RANDOM%9+1))" "$((1735700000+id))"; }
# v1 20-col: same but NO trailing epoch
v1r20(){ local id=$1 sz=${2:-5000} m=${3:-3800} st=${4:-COMPLETED} pg=${5:-bwa} rf=${6:-hg38} fl=${7:-}
  : "${fl:=1.0.$pg.$id}"
  printf '%s,%s,0,0,4096,120,%s,95,%s,u,,%s,%s,%s,1,0,0,0,0,2026-01-0%s\n' \
  "$id" "$sz" "$m" "$st" "$pg" "$rf" "$fl" "$((RANDOM%9+1))"; }
# v2 17-col via the shipped builder (the only sanctioned v2 writer)
v2r(){ local id=$1 sz=${2:-5000} m=${3:-3800} st=${4:-COMPLETED} pg=${5:-bwa} rf=${6:-hg38} fl=${7:-} pp=${8:-rnaseq} sm=${9:-}
  : "${fl:=1.0.$pg.$id}"
  jrRow pipeline="$pp" step="$pg" ref="$rf" sample="$sm" \
  jobID="$id" inputSize="$sz" memUsed="$m" memAllocated="$m" \
  timeUsed=95 timeAllocated=120 cores=1 status="$st" user=u node=n \
  epoch="$((1735700000+id))" isoDate=2026-01-03 flag="$fl"; }
H_V2="pipeline,step,ref,sample,jobID,inputSize,memUsed,memAllocated,timeUsed,timeAllocated,cores,status,user,node,epoch,isoDate,flag"
H_V1_21="1jobID,2inputSize,3memDefault,4timeDefault,5memAllocated,6timeAllocated,7memUsed,8timeUsed,9status,10user,11saccMem,12program,13ref,14flag,15cores,16extraMem,17extraTime,18ratioM,19ratioT,20date,21epoch"
H_V1_18="1jobID,2inputSize,3memDefault,4timeDefault,5memAllocated,6timeAllocated,7memUsed,8timeUsed,9status,10user,11program,12ref,13flag,14cores,15extraMem,16extraTime,17date"

# jobIDs present (by the v2 jobID column) in a migrated file, sorted + trimmed
v2ids(){ awk -F, '$1!="pipeline" && $1 !~ /^[0-9]+$/ && NF==17 {print $5}' "$1" | sort -n | xargs; }
norm(){ printf '%s\n' $1 | sort -n | xargs; }   # normalize an id list (trim/sort)
# SAFETY-CRITICAL post-condition: NO raw v1-shaped row survives (numeric jobID,
# 20/21 cols). A v1 row read through the v2 column map is what the prune deletes.
# (Inert junk that readers skip is tolerated; a clean v1/mixed migration has none
# anyway, since jrEmitV2 emits only header + mapped/passed rows.)
no_v1_rows(){ ! awk -F, 'NR>1 && $1 ~ /^[0-9]+$/ && (NF==20||NF==21){f=1} END{exit !f}' "$1"; }

# INVARIANT runner: given a store file $1 and the space-sep list of jobIDs that
# MUST survive $2 (the convertible rows), migrate and assert the full contract.
invariant(){
  local f="$1" expect="$2" label="$3"
  local before_ids; before_ids=$(norm "$expect")
  jrMigrateStore "$f" 2>/dev/null
  chk "[ \"\$(jrSchemaVersion '$f')\" = 2 ]"          "$label: ends as v2"
  chk "head -1 '$f' | grep -qx '$H_V2'"               "$label: v2 header present"
  chk "no_v1_rows '$f'"                               "$label: no raw v1 row survives (prune-safe)"
  local got; got=$(v2ids "$f")
  chk "[ \"$got\" = \"$before_ids\" ]"                "$label: exactly the convertible jobIDs survive (want '$before_ids' got '$got')"
  # idempotency: a second migrate changes nothing and makes no new backup
  local n1 n2 sum1 sum2
  n1=$(ls "$f".back.*.txt 2>/dev/null | wc -l); sum1=$(md5sum "$f"|cut -d' ' -f1)
  jrMigrateStore "$f" 2>/dev/null
  n2=$(ls "$f".back.*.txt 2>/dev/null | wc -l); sum2=$(md5sum "$f"|cut -d' ' -f1)
  chk "[ '$sum1' = '$sum2' ]"                         "$label: idempotent (byte-identical on 2nd migrate)"
  chk "[ '$n1' = '$n2' ]"                             "$label: idempotent (no extra backup)"
}

# field-mapping fidelity: ALL 17 v2 columns from one v1 row, with a DISTINCT
# sentinel per v1 position, so any mis-map or swap is visible (not just the few
# columns previously checked). Covers the columns the estimator depends on
# (memAllocated/timeAllocated/cores) and, critically, EPOCH -- which a {n,} interval
# regex silently blanked on mawk. Run under both awks (see the cross-awk harness).
echo "== field mapping fidelity (all 17 columns) =="
f="$TMP/fm.txt"
# v1 21-col positions: 1jobID 2inputSize 3memDef 4timeDef 5memAlloc 6timeAlloc
#   7memUsed 8timeUsed 9status 10user 11saccMem 12program 13ref 14flag 15cores
#   16extraMem 17extraTime 18ratioM 19ratioT 20date 21epoch
printf '7001,IN,MD,TD,MA5,TA6,MU7,TU8,COMPLETED,USER,SACC,minimap2,mm39,3.1.minimap2.sampX,CORE15,EM,ET,RM,RT,DATE20,1735712345\n' > "$f"
jrMigrateStore "$f" 2>/dev/null
row=$(awk -F, 'NR==2' "$f")
chk "[ \"\$(echo '$row'|cut -d, -f1)\"  = legacy ]"     "fm col1  pipeline=legacy"
chk "[ \"\$(echo '$row'|cut -d, -f2)\"  = minimap2 ]"   "fm col2  step <- program(\$12)"
chk "[ \"\$(echo '$row'|cut -d, -f3)\"  = mm39 ]"       "fm col3  ref <- \$13"
chk "[ \"\$(echo '$row'|cut -d, -f4)\"  = sampX ]"      "fm col4  sample <- flag"
chk "[ \"\$(echo '$row'|cut -d, -f5)\"  = 7001 ]"       "fm col5  jobID <- \$1"
chk "[ \"\$(echo '$row'|cut -d, -f6)\"  = IN ]"         "fm col6  inputSize <- \$2"
chk "[ \"\$(echo '$row'|cut -d, -f7)\"  = MU7 ]"        "fm col7  memUsed <- \$7"
chk "[ \"\$(echo '$row'|cut -d, -f8)\"  = MA5 ]"        "fm col8  memAllocated <- \$5"
chk "[ \"\$(echo '$row'|cut -d, -f9)\"  = TU8 ]"        "fm col9  timeUsed <- \$8"
chk "[ \"\$(echo '$row'|cut -d, -f10)\" = TA6 ]"        "fm col10 timeAllocated <- \$6"
chk "[ \"\$(echo '$row'|cut -d, -f11)\" = CORE15 ]"     "fm col11 cores <- \$15"
chk "[ \"\$(echo '$row'|cut -d, -f12)\" = COMPLETED ]"  "fm col12 status <- \$9"
chk "[ \"\$(echo '$row'|cut -d, -f13)\" = USER ]"       "fm col13 user <- \$10"
chk "[ -z \"\$(echo '$row'|cut -d, -f14)\" ]"           "fm col14 node <- empty (v1 had none)"
chk "[ \"\$(echo '$row'|cut -d, -f15)\" = 1735712345 ]" "fm col15 EPOCH <- \$21 (mawk-safe; not blanked)"
chk "[ \"\$(echo '$row'|cut -d, -f16)\" = DATE20 ]"     "fm col16 isoDate <- \$20"
chk "[ \"\$(echo '$row'|cut -d, -f17)\" = 3.1.minimap2.sampX ]" "fm col17 flag <- \$14"

echo "== adversarial matrix =="
# 1 pure v1 headerless 21
f="$TMP/a1.txt"; { v1r21 101; v1r21 102; } > "$f"; invariant "$f" "101 102" "pure-v1-21-headerless"
# 2 pure v1 headerless 20
f="$TMP/a2.txt"; { v1r20 201; v1r20 202; } > "$f"; invariant "$f" "201 202" "pure-v1-20-headerless"
# 3 pure v1 with 21-col header
f="$TMP/a3.txt"; { echo "$H_V1_21"; v1r21 301; v1r20 302; } > "$f"; invariant "$f" "301 302" "v1-21hdr-mixed-width"
# 4 pure v1 with stale 18-col header  (header must be skipped, not converted)
f="$TMP/a4.txt"; { echo "$H_V1_18"; v1r21 401; } > "$f"; invariant "$f" "401" "v1-18hdr"
# 5 pure v2 (sanctioned)
f="$TMP/a5.txt"; { echo "$H_V2"; v2r 501; v2r 502; } > "$f"; invariant "$f" "501 502" "pure-v2"
# 6 MIXED v1-first then appended v2 (the ld32 old-writer vector)
f="$TMP/a6.txt"; { v1r21 601; v1r20 602; v2r 603; } > "$f"; invariant "$f" "601 602 603" "mixed-v1then-v2-headerless"
# 7 MIXED with v1 header then v1 rows then appended v2
f="$TMP/a7.txt"; { echo "$H_V1_21"; v1r21 701; v2r 702; } > "$f"; invariant "$f" "701 702" "mixed-v1hdr-then-v2"
# 8 *** HYBRID the other way: v2 header + v2 rows, then a v1 row appended ***
#    (old node writing v1 into a freshly-migrated v2 store during rollout)
f="$TMP/a8.txt"; { echo "$H_V2"; v2r 801; v1r21 802; } > "$f"; invariant "$f" "801 802" "HYBRID-v2then-appended-v1"
# 9 HYBRID: v2 header, v1 row FIRST among data, then v2
f="$TMP/a9.txt"; { echo "$H_V2"; v1r20 901; v2r 902; } > "$f"; invariant "$f" "901 902" "HYBRID-v2hdr-v1first"
# 10 interleaved
f="$TMP/a10.txt"; { v1r21 1001; v2r 1002; v1r20 1003; v2r 1004; } > "$f"; invariant "$f" "1001 1002 1003 1004" "interleaved"
# 11 non-numeric inputSize carried (missingInputFile) -- must survive, not drop
f="$TMP/a11.txt"; v1r21 1101 missingInputFile 0 COMPLETED bwa hg38 "1.0.bwa.1101" > "$f"
invariant "$f" "1101" "nonnumeric-inputSize"
chk "awk -F, 'NR==2{exit !(\$6==\"missingInputFile\")}' '$f'" "nonnumeric-inputSize: value carried verbatim"
# 12 blank lines + trailing whitespace interspersed
f="$TMP/a12.txt"; { v1r21 1201; echo ""; v2r 1202; echo "   "; } > "$f"; invariant "$f" "1201 1202" "blank-lines"
# 13 dotless flag (no sample) and fanout flag
f="$TMP/a13.txt"; { v1r21 1301 5000 3800 COMPLETED bwa hg38 "bwa"; v1r21 1302 5000 3800 COMPLETED bwa hg38 "2.1.bwa.sampA.chr1"; } > "$f"
invariant "$f" "1301 1302" "dotless+fanout-flags"

echo "== prune self-destruction (the real cleanUp branch selection) =="
# Replicate cleanUp.sh lines 283-288 EXACTLY, then assert a COMPLETED v2 row is
# never deleted by the v1-positional branch running on v2 data.
prune(){ # $1=file ; mutates in place like cleanUp
  local f=$1 t; t=$(date -d '7 days ago' +%s)
  if [ "$(jrSchemaVersion "$f")" = 2 ]; then
    awk -F, -v t="$t" 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; print; next} ($(c["epoch"])+0 > t) || ($(c["status"])=="COMPLETED")' "$f" > "$f.new" && mv "$f.new" "$f"
  else
    awk -F"," -v t="$t" '{ if($21 > t || $9 == "COMPLETED") print $0 }' "$f" > "$f.new" && mv "$f.new" "$f"
  fi
}
# In production, cleanUp calls jrMigrateStore BEFORE this prune (line 67). Model
# that: migrate, then prune, then assert COMPLETED rows all survived.
for name in a6 a7 a8 a9 a10; do
  f="$TMP/p_$name.txt"
  case $name in
    a6) { v1r21 601 5000 3800 COMPLETED; v1r20 602 5000 3800 COMPLETED; v2r 603 5000 3800 COMPLETED; } > "$f"; want="601 602 603";;
    a7) { echo "$H_V1_21"; v1r21 701 5000 3800 COMPLETED; v2r 702 5000 3800 COMPLETED; } > "$f"; want="701 702";;
    a8) { echo "$H_V2"; v2r 801 5000 3800 COMPLETED; v1r21 802 5000 3800 COMPLETED; } > "$f"; want="801 802";;
    a9) { echo "$H_V2"; v1r20 901 5000 3800 COMPLETED; v2r 902 5000 3800 COMPLETED; } > "$f"; want="901 902";;
    a10) { v1r21 1001 5000 3800 COMPLETED; v2r 1002 5000 3800 COMPLETED; v1r20 1003 5000 3800 COMPLETED; } > "$f"; want="1001 1002 1003";;
  esac
  jrMigrateStore "$f" 2>/dev/null
  prune "$f"
  got=$(v2ids "$f")
  wsorted=$(norm "$want")
  chk "[ \"$got\" = \"$wsorted\" ]" "prune-survival[$name]: all COMPLETED rows kept (want '$wsorted' got '$got')"
done

echo "== concurrency: N parallel migrations on one store =="
f="$TMP/conc.txt"; { for i in $(seq 1 30); do v1r21 $((2000+i)); done; v2r 2999; } > "$f"
for i in $(seq 1 12); do ( jrMigrateStore "$f" 2>/dev/null ) & done; wait
chk "[ \"\$(jrSchemaVersion '$f')\" = 2 ]"  "concurrency: store is v2"
chk "no_v1_rows '$f'"                       "concurrency: no raw v1 row survives"
got=$(v2ids "$f"); want=$(norm "$(seq 2001 2030; echo 2999)")
chk "[ \"$got\" = \"$want\" ]"             "concurrency: all 31 rows survive exactly once (got count $(echo $got|wc -w))"
nb=$(ls "$f".back.*.txt 2>/dev/null | wc -l)
chk "[ '$nb' -ge 1 ]"                      "concurrency: at least one backup made"
# no leftover temp/lock
chk "! ls '$f'.v2.* >/dev/null 2>&1"       "concurrency: no leftover .v2.PID temp"
chk "[ ! -d '$f'.lock ]"                   "concurrency: lock released"

echo "== concurrency: N parallel APPENDS to one store (locked, no lost/corrupt rows) =="
# K workers each append M unique rows via jrAppendRow. Every row must land exactly
# once, the store must stay valid v2 with one header, and no row may be corrupted
# by an interleaved write. This is the write-path analogue of the migration race.
f="$TMP/app.txt"; : > "$f"; K=10; M=15
for w in $(seq 1 $K); do
  ( for j in $(seq 1 $M); do
      jrAppendRow "$f" pipeline=rnaseq step=bwa ref=hg38 jobID=$((w*1000+j)) \
        inputSize=$((1000+j)) memUsed=$((100+j)) status=COMPLETED flag="1.0.bwa.s$w" 2>/dev/null
    done ) &
done; wait
exp=$(norm "$(for w in $(seq 1 $K); do for j in $(seq 1 $M); do echo $((w*1000+j)); done; done)")
got=$(v2ids "$f")
chk "[ \"\$(jrSchemaVersion '$f')\" = 2 ]"  "concurrent-append: store is valid v2"
chk "[ \"\$(grep -c '^pipeline,step,ref' '$f')\" = 1 ]" "concurrent-append: exactly one header line"
chk "[ \"$got\" = \"$exp\" ]"              "concurrent-append: all $((K*M)) rows present exactly once (got $(echo $got|wc -w))"
chk "! awk -F, 'NR>1 && NF!=17{bad=1} END{exit !bad}' '$f'" "concurrent-append: every row is well-formed 17-col (no interleave)"
chk "[ ! -d '$f'.lock ]"                   "concurrent-append: lock released"

echo "== chaos: concurrent appends + migrations + prunes on ONE store =="
# Mixed workload hammering a single store: appenders (all COMPLETED, so a recency
# prune must keep them), a migrator, and pruners (keep recent-or-COMPLETED). The
# invariant: every appended COMPLETED row is still present, store stays valid v2.
f="$TMP/chaos.txt"; { v1r21 500; v1r20 501; } > "$f"   # start from a legacy store
tnow=$(date -d '7 days ago' +%s 2>/dev/null || echo 0)
for w in $(seq 1 8); do
  ( for j in $(seq 1 10); do
      jrAppendRow "$f" pipeline=p step=bwa ref=hg38 jobID=$((w*1000+j)) inputSize=$j memUsed=$j status=COMPLETED flag=1.0.bwa.c 2>/dev/null
    done ) &
done
for _ in 1 2 3; do ( jrRewrite "$f" 'if(($(col["epoch"])+0>t)||($(col["status"])=="COMPLETED"))print' -v t="$tnow" 2>/dev/null ) & done
( jrMigrateStore "$f" 2>/dev/null ) &
wait
expc=$(norm "$(for w in $(seq 1 8); do for j in $(seq 1 10); do echo $((w*1000+j)); done; done)")
gotc=$(v2ids "$f" | tr ' ' '\n' | awk '$1>=1000' | xargs)   # only the appended ids (legacy 500/501 may age out)
chk "[ \"\$(jrSchemaVersion '$f')\" = 2 ]"  "chaos: store ends valid v2"
chk "no_v1_rows '$f'"                       "chaos: no raw v1 row survives the mixed workload"
chk "! awk -F, 'NR>1 && NF!=17{b=1} END{exit !b}' '$f'" "chaos: no corrupt/interleaved row"
chk "[ \"$gotc\" = \"$expc\" ]"            "chaos: all 80 appended COMPLETED rows survive (got $(echo $gotc|wc -w))"
chk "[ ! -d '$f'.lock ]"                   "chaos: lock released"

echo "== all-junk file is NOT silently emptied =="
f="$TMP/junk.txt"; { echo "total,garbage,here"; echo "another junk line"; } > "$f"
cp "$f" "$TMP/junk.orig"
err=$(jrMigrateStore "$f" 2>&1 1>/dev/null)
# an unrecognized store must be LEFT UNTOUCHED and WARNED ABOUT (not silent).
chk "diff -q '$f' '$TMP/junk.orig' >/dev/null"  "all-junk: left untouched (no conversion)"
chk "echo \"$err\" | grep -qi 'no known schema'" "all-junk: warned on stderr (not silent)"

echo "== non-silent: malformed rows inside a v1 store are reported, kept in backup =="
f="$TMP/mal.txt"
{ v1r21 7101; printf '7102,5000,0,0,4096,120,3800,95,COMPLETED,u,,bwa,hg38,1.0.bwa.s,1,0,0\n'; v2r 7103; } > "$f"  # middle row = 19 col
nbad=$(jrCountUnrecognized "$f")
chk "[ '$nbad' = 1 ]"                            "malformed: jrCountUnrecognized sees the 1 bad row"
err=$(jrMigrateStore "$f" 2>&1 1>/dev/null)
chk "[ \"\$(jrSchemaVersion '$f')\" = 2 ]"       "malformed: store still migrates to v2"
chk "echo \"$err\" | grep -qiE 'WARNING.*1 row.*no known schema'" "malformed: migration WARNS about the dropped row"
got=$(v2ids "$f"); chk "[ \"$got\" = \"$(norm '7101 7103')\" ]" "malformed: the two GOOD rows survive (bad one dropped)"
bak=$(ls "$f".back.*.txt 2>/dev/null | head -1)
chk "[ -n '$bak' ] && awk -F, '\$1==7102' '$bak' | grep -q ." "malformed: bad row recoverable in the backup"
# the dropped row is ALSO saved to a .unrecognized file (easier than diffing backup)
unrecf=$(ls "$f".unrecognized.*.txt 2>/dev/null | head -1)
chk "[ -n '$unrecf' ] && grep -q 7102 '$unrecf' && [ \"\$(grep -c . '$unrecf')\" = 1 ]" "malformed: dropped row saved to a .unrecognized file (just the bad line)"
# jrLint reports it and returns non-zero
jrLint "$bak" >/dev/null; lrc=$?; chk "[ $lrc -ne 0 ]" "jrLint flags a store with unrecognized rows (non-zero)"

echo "== post-conversion integrity gate: a broken conversion ABORTS, store untouched =="
# Force jrEmitV2 to drop a recognizable row; the gate must refuse to commit, leave
# the live store byte-identical, remove the redundant backup, and error loudly.
f="$TMP/gate.txt"
{ v1r21 8101; v1r21 8102; } > "$f"
before=$(md5sum "$f" | cut -d' ' -f1)
( jrEmitV2() { printf '%s\n' "$JR_HEADER_V2"; echo "rnaseq,bwa,hg38,,8101,5000,3800,4096,95,120,1,COMPLETED,u,n,1,2026-01-01,1.0.bwa.s"; }   # drops 8102
  gerr=$(jrMigrateStore "$f" 2>&1 1>/dev/null)
  [ "$(md5sum "$f" | cut -d' ' -f1)" = "$before" ] || { echo "LEAK"; exit 7; }
  echo "$gerr" | grep -qiE 'ERROR.*integrity.*UNTOUCHED' || { echo "NOWARN"; exit 8; }
  [ -z "$(ls "$f".back.*.txt 2>/dev/null)" ] || { echo "BAKLEFT"; exit 9; } )
grc=$?
chk "[ $grc -eq 0 ]" "integrity gate: broken conversion aborts, store untouched, backup removed, errors loudly (rc=$grc)"
chk "[ \"\$(jrSchemaVersion '$f')\" != 0 ]" "integrity gate: store is still its original valid self (not corrupted)"

echo "== item B: all-digit pipeline name + header-word collision =="
# A 17-col row with an all-digit field 1 is ambiguous: a v2 row whose pipeline name
# is all-digit, OR a corrupt/truncated v1 (jobID-first) row. We disambiguate by
# HEADER CONTEXT, never by guessing: under a v2 header it is a v2 record (and must
# be preserved -- it used to vanish silently); headerless it stays UNRECOGNIZED.
adr="2024,bwa,hg38,sA,55555,1000,512,1024,60,120,4,COMPLETED,u,n,1700000000,2026-01-01,1.0.bwa.c"
f="$TMP/bdigit_hdr.txt"; { echo "$H_V2"; echo "$adr"; } > "$f"
chk "[ \"\$(jrSchemaVersion '$f')\" = 2 ]"        "itemB: headed all-digit-pipeline store is v2"
chk "[ \"\$(jrCountRecognized '$f')\" = 1 ]"      "itemB: headed all-digit-pipeline row is RECOGNIZED (not dropped)"
chk "[ \"\$(jrCountUnrecognized '$f')\" = 0 ]"    "itemB: headed all-digit-pipeline row is not flagged unrecognized"
chk "jrEmitV2 '$f' | grep -qxF '$adr'"            "itemB: jrEmitV2 passes the all-digit-pipeline v2 row through verbatim"
jrLint "$f" >/dev/null; chk "[ \$? -eq 0 ]"       "itemB: jrLint calls the headed all-digit-pipeline store clean"
f="$TMP/bdigit_bare.txt"; echo "$adr" > "$f"
chk "[ \"\$(jrCountUnrecognized '$f')\" = 1 ]"    "itemB: headerLESS all-digit \$1 row stays UNRECOGNIZED (no silent guess)"
# A header-shaped EXTRA line (stray duplicate header, or a data row whose
# pipeline+step literally equal the header words "pipeline"/"step") is byte-identical
# to the header and its data is silently excluded -- jrLint must call that out loudly.
f="$TMP/bdup.txt"; { echo "$H_V2"; echo "$H_V2"; echo "rnaseq,bwa,hg38,,100,1,1,1,1,1,1,COMPLETED,u,n,1,2026-01-01,1.0.bwa.c"; } > "$f"
chk "jrLint '$f' | grep -qi 'DUPLICATE HEADER'"   "itemB: jrLint reports a duplicate/header-collision line"
jrLint "$f" >/dev/null; chk "[ \$? -ne 0 ]"       "itemB: jrLint returns non-zero on duplicate header"

echo "== curation pointer: migration directs users to curate legacy history =="
# A v1->v2 migration tags old rows pipeline=legacy (v1 had no pipeline column). The
# migration must POINT the user at the (optional) per-pipeline curation step so the v2
# UX is not a dead end -- but only when legacy rows were actually produced, and it must
# name the jrCurate.sh command only when that tool is actually installed.
f="$TMP/curptr_v1.txt"; v1r21 9001 > "$f"; v1r20 9002 >> "$f"
err=$(jrMigrateStore "$f" 2>&1 1>/dev/null)
chk "echo \"\$err\" | grep -qi 'per-pipeline performance after migrating'" "curptr: v1 migration emits the curation pointer"
chk "echo \"\$err\" | grep -q 'pipeline=legacy'"                           "curptr: pointer explains the pipeline=legacy tagging"
# tool ABSENT -> generic wording, no bare command promise. Scrub PATH to coreutils so
# this is deterministic whether or not SMARTSLURM_ROOT/bin (which may now ship
# jrCurate.sh and be on PATH via ss_activate) is present -- we test the message branch,
# not the ambient environment.
f="$TMP/curptr_abs.txt"; v1r21 9050 > "$f"
err_abs=$(PATH=/usr/bin:/bin; jrMigrateStore "$f" 2>&1 1>/dev/null)
chk "! echo \"\$err_abs\" | grep -q 'jrCurate.sh export'"                  "curptr: generic wording (no command) when jrCurate is not on PATH"
# a migration that produces NO legacy rows (pure v2 passthrough) stays silent on curation
f="$TMP/curptr_none.txt"; v2r 9100 > "$f"   # headerless single v2 row -> migrates, no legacy
err=$(jrMigrateStore "$f" 2>&1 1>/dev/null)
chk "! echo \"\$err\" | grep -qi 'per-pipeline performance after migrating'" "curptr: no curation pointer when migration produced no legacy rows"
# tool PRESENT -> the pointer names the actual command
stubd="$TMP/curptr_bin"; mkdir -p "$stubd"; printf '#!/bin/sh\n' > "$stubd/jrCurate.sh"; chmod +x "$stubd/jrCurate.sh"
f="$TMP/curptr_installed.txt"; v1r21 9201 > "$f"
err=$(PATH="$stubd:$PATH" jrMigrateStore "$f" 2>&1 1>/dev/null)
chk "echo \"\$err\" | grep -q 'jrCurate.sh export'"                        "curptr: names 'jrCurate.sh export' when the tool is installed"

echo "== jrRewrite return contract: 0 iff applied, non-zero if it did not write =="
# A caller that reports success to a user (e.g. curation import) must be able to tell
# an APPLIED rewrite from a silent skip. Previously jrRewrite returned the trailing
# _jrUnlock (rmdir) status = 0 on every path, so a lock-timeout or awk failure looked
# like success (adversary Finding 1). Contract: 0 only when the atomic rename happened
# (or the store was empty); non-zero otherwise.
f="$TMP/rw_ok.txt"; { echo "$H_V2"; v2r 7701 > /dev/null; echo "p,bwa,hg38,,7701,1,1,1,1,1,1,COMPLETED,u,n,1,2026-01-01,1.0.bwa.c"; } > "$f"
jrRewrite "$f" 'print' >/dev/null 2>&1; chk "[ \$? -eq 0 ]"          "rwret: applied rewrite returns 0"
# lock HELD by someone else -> jrRewrite must NOT write and must return non-zero
cp "$f" "$TMP/rw_held.orig"; mkdir "$f.lock"
JR_LOCK_TRIES=2 JR_LOCK_STALE=99999 jrRewrite "$f" 'print' >/dev/null 2>&1; rc=$?
chk "[ $rc -ne 0 ]"                                                   "rwret: lock-timeout returns non-zero (did not write)"
chk "diff -q '$f' '$TMP/rw_held.orig' >/dev/null"                    "rwret: lock-timeout left the store byte-identical"
rmdir "$f.lock" 2>/dev/null
# empty store -> honest no-op success
: > "$TMP/rw_empty.txt"; jrRewrite "$TMP/rw_empty.txt" 'print' >/dev/null 2>&1; chk "[ \$? -eq 0 ]" "rwret: empty store is a no-op success (0)"
# unrecognized store (cannot become v2) -> non-zero (nothing applied)
echo "total,garbage,here" > "$TMP/rw_junk.txt"; jrRewrite "$TMP/rw_junk.txt" 'print' >/dev/null 2>&1; chk "[ \$? -ne 0 ]" "rwret: unrecognized store returns non-zero (not applied)"

echo "== deterministic fuzzer (seed=$SEED) =="
for n in $(seq 1 120); do
  f="$TMP/fz_$n.txt"; : > "$f"; ids=""
  # maybe a header
  case $((RANDOM%4)) in 0) echo "$H_V2" >>"$f";; 1) echo "$H_V1_21" >>"$f";; 2) echo "$H_V1_18" >>"$f";; esac
  rows=$((RANDOM%8+1))
  for r in $(seq 1 $rows); do
    jid=$(( n*1000 + r ))
    case $((RANDOM%5)) in
      0) v1r21 $jid >>"$f";;
      1) v1r20 $jid >>"$f";;
      2) v2r  $jid >>"$f";;
      3) v1r21 $jid missingInputFile >>"$f";;
      4) echo "" >>"$f"; continue;;  # blank line, not a record
    esac
    ids="$ids $jid"
  done
  # sprinkle junk sometimes
  [ $((RANDOM%3)) = 0 ] && echo "JUNK,row,$n" >>"$f"
  want=$(norm "$ids")
  jrMigrateStore "$f" 2>/dev/null
  v=$(jrSchemaVersion "$f"); got=$(v2ids "$f")
  if [ -z "$want" ]; then
    # no records seeded (blank/junk only): nothing to preserve, just assert safety
    if no_v1_rows "$f"; then ok; else no "fuzz#$n: raw v1 left in a record-less file"; fi
  elif [ "$v" != 2 ] || ! no_v1_rows "$f"; then no "fuzz#$n: not v2/raw-v1-left (v=$v)";
  elif [ "$got" != "$want" ]; then no "fuzz#$n: id set changed (want '$want' got '$got')";
  else ok; fi
done

echo
echo "================  RESULT  ================"
echo "PASS=$PASS  FAIL=$FAIL"
[ "$FAIL" -eq 0 ] && echo "ALL GREEN" || { printf '%s\n' "${FAILS[@]}" | sort | uniq -c | sort -rn | head -40; exit 1; }
