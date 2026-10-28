#!/bin/bash
# ---------------------------------------------------------------------------
# harness.sh  --  activate a Slurm-free environment for running SmartSlurm.
#
# Source this from a test or an interactive shell:
#     source test/lib/harness.sh
#     ss_activate                 # set up PATH, sh->bash, scratch state
#     ss_run_pipeline <script> <sbatchOpts> [--dryrun]
#
# It provides:
#   ss_activate            - prepare environment (idempotent)
#   ss_new_workspace <dir> - fresh working dir with clean state
#   ss_run_pipeline ...    - run runAsPipeline under the mocks
#   assertion helpers      - assert_ok / assert_file / assert_empty / etc.
#
# WHY sh->bash: SmartSlurm scripts declare `#!/bin/sh` but use bash features.
# On the production RHEL cluster /bin/sh IS bash, so this is faithful emulation,
# not a hack. On Debian/Ubuntu sandboxes /bin/sh is dash and the scripts break;
# ss_activate repoints /bin/sh at bash when it can (root), else warns.
# ---------------------------------------------------------------------------
set -uo pipefail

# This file lives at <harness>/test/lib/harness.sh, so TEST_DIR is its grandparent.
_HARNESS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export TEST_DIR="$(cd "$_HARNESS_DIR/.." && pwd)"
export MOCK_BIN="$TEST_DIR/mockSlurm"

# SMARTSLURM_ROOT = the SmartSlurm checkout whose bin/ holds the tools under test.
# Resolve it robustly so the harness works at ANY nesting depth:
#   - in-tree  : SmartSlurm/slurmless/test/  -> root is two levels up
#   - legacy   : SmartSlurm/test/            -> root is one level up
#   - external : an out-of-tree SmartSlurm clone -> set SMARTSLURM_ROOT yourself
# Honor an explicit SMARTSLURM_ROOT; otherwise walk up from TEST_DIR looking for
# the dir that actually contains bin/runAsPipeline (the definitive marker), so the
# path is correct no matter how deep slurmless sits. Fall back to the
# legacy one-level-up if no marker is found (e.g. a dry run with no tools present).
if [ -z "${SMARTSLURM_ROOT:-}" ]; then
  _d="$TEST_DIR"
  while [ "$_d" != "/" ]; do
    if [ -e "$_d/bin/runAsPipeline" ]; then SMARTSLURM_ROOT="$_d"; break; fi
    _d="$(dirname "$_d")"
  done
  SMARTSLURM_ROOT="${SMARTSLURM_ROOT:-$(cd "$TEST_DIR/.." && pwd)}"
fi
export SMARTSLURM_ROOT

# ss_activate: prepare a Slurm-free environment. It has two environment-mutating
# side effects appropriate for a disposable sandbox or CI runner but surprising on
# a shared workstation; each can be opted out of (the harness then expects the
# host to already satisfy it):
#   SS_NO_INSTALL=1     do not apt-get install gnuplot (curve-fit blocks skip)
#   SS_NO_SH_REPOINT=1  do not repoint /bin/sh -> bash (assume sh is already bash)
# A one-line summary of what it changed is printed at the end.
ss_activate() {
  local _gp="present" _sh="already-bash"

  # 0) gnuplot is required by estimateResource.sh for curve fitting. Install it
  #    if missing (sandbox images lose apt packages on reset). Non-fatal: the
  #    curve-fit blocks guard on it, functional blocks don't need it.
  if ! command -v gnuplot >/dev/null 2>&1; then
    if [ -n "${SS_NO_INSTALL:-}" ]; then
      _gp="missing(SS_NO_INSTALL)"
    else
      echo "[harness] gnuplot missing; attempting install (needed for curve fitting; set SS_NO_INSTALL=1 to skip)..."
      if { sudo apt-get update -qq && sudo apt-get install -y --no-install-recommends gnuplot-nox; } >/dev/null 2>&1 \
         || { apt-get update -qq && apt-get install -y --no-install-recommends gnuplot-nox; } >/dev/null 2>&1; then
        _gp="installed"
      else
        _gp="missing"
        echo "[harness] WARNING: could not install gnuplot; curve-fit blocks will skip." >&2
      fi
    fi
  fi

  # 1) sh -> bash (match production). Only if currently dash and we can fix it.
  if [ "$(readlink -f /bin/sh 2>/dev/null)" != "$(readlink -f /bin/bash)" ]; then
    if [ -n "${SS_NO_SH_REPOINT:-}" ]; then
      _sh="dash(SS_NO_SH_REPOINT)"
      echo "[harness] WARNING: /bin/sh is not bash and SS_NO_SH_REPOINT is set;" >&2
      echo "[harness]          SmartSlurm's #!/bin/sh scripts may fail under dash." >&2
    elif ln -sf bash /bin/sh 2>/dev/null; then
      ln -sf bash /usr/bin/sh 2>/dev/null || true
      _sh="repointed"
      echo "[harness] repointed /bin/sh -> bash (matches production RHEL)"
    else
      _sh="dash(no-root)"
      echo "[harness] WARNING: /bin/sh is not bash and cannot repoint (need root)." >&2
      echo "[harness]          SmartSlurm's #!/bin/sh scripts may fail under dash." >&2
    fi
  fi

  # 2) Stage the mock commands into a NEUTRAL path (no "SmartSlurm" in it).
  #    ssbatch finds the real sbatch with `which -a sbatch | grep -v SmartSlurm`,
  #    so a mock living under the repo's test/mockSlurm/ (which contains the
  #    substring "SmartSlurm") would be filtered out and realSbatch would be
  #    empty. Copying the shims to $SS_SCRATCH/slurm-shims avoids that.
  export SHIM_BIN="${SS_SCRATCH:-/tmp/ss-scratch}/slurm-shims"
  mkdir -p "$SHIM_BIN"
  cp -f "$MOCK_BIN"/* "$SHIM_BIN"/ 2>/dev/null
  chmod +x "$SHIM_BIN"/* 2>/dev/null
  # Keep only real Slurm-command shims on PATH; drop docs/helpers that share the dir.
  rm -f "$SHIM_BIN/README.md" "$SHIM_BIN/check_sacct_format.sh"

  # 3) PATH: staged mocks first, then SmartSlurm/bin, then the rest
  export PATH="$SHIM_BIN:$SMARTSLURM_ROOT/bin:$PATH"

  # 3) Redirect HOME into the scratch tree. config.txt sets
  #       export smartSlurmJobRecordDir=$HOME/.smartSlurm
  #    UNCONDITIONALLY, and ssbatch/cleanUp.sh re-source config in subprocesses,
  #    so exporting smartSlurmJobRecordDir from here does NOT survive -- it gets
  #    clobbered back to $HOME/.smartSlurm. Pointing HOME at the scratch dir is
  #    the faithful fix (it is exactly the path production resolves) and keeps the
  #    learning DB, stats, and any user config off a real ~/.smartSlurm. We also
  #    set USER (workflowPlot.py and cleanUp.sh read $USER).
  export HOME="${SS_SCRATCH:-/tmp/ss-scratch}/home"
  export USER="${USER:-sstester}"
  mkdir -p "$HOME"
  export smartSlurmJobRecordDir="$HOME/.smartSlurm"
  mkdir -p "$smartSlurmJobRecordDir/stats"
  # shellcheck disable=SC1090
  source "$SMARTSLURM_ROOT/config/config.txt"
  export smartSlurmJobRecordDir="$HOME/.smartSlurm"   # match config's resolution

  # 4) mock scheduler state
  export MOCK_SLURM_STATE="${SS_SCRATCH:-/tmp/ss-scratch}/mockstate"
  mkdir -p "$MOCK_SLURM_STATE/jobs"

  # Concise summary of what was set up / changed, so a human running this
  # interactively can see exactly what the harness touched.
  echo "[harness] activated: SMARTSLURM_ROOT=$SMARTSLURM_ROOT | /bin/sh=$_sh | gnuplot=$_gp | HOME=$HOME (scratch)"
}

# create (or reset) a clean workspace and cd into it
ss_new_workspace() {
  local dir="$1"
  rm -rf "$dir"; mkdir -p "$dir"; cd "$dir" || return 1
  # fresh mock state per workspace so job ids/inject don't leak between tests
  export MOCK_SLURM_STATE="$dir/.mockstate"; mkdir -p "$MOCK_SLURM_STATE/jobs"
  rm -rf "$smartSlurmJobRecordDir"; mkdir -p "$smartSlurmJobRecordDir/stats"
  # Pipelines call their step tools by bare name (e.g. `findNumber.sh`), relying
  # on PATH. Put the workspace itself FIRST on PATH so fixture-local tools win
  # over any same-named tool in the repo's bin, keeping fixtures self-contained.
  case ":$PATH:" in *":$dir:"*) ;; *) export PATH="$dir:$PATH";; esac
}

# --- resource model selection -------------------------------------------------
# ss_use_model NAME|PATH : activate a resource model so the mock sbatch/sacct/
#   memCpuMonitor synthesize mem/time from input size. NAME resolves to
#   $TEST_DIR/models/NAME.model; an explicit path is used as-is.
# ss_clear_model         : deactivate (mocks revert to fixed/real behavior).
ss_use_model() {
  local m="$1"
  [ -f "$m" ] || m="$TEST_DIR/models/$1.model"
  [ -f "$m" ] || { echo "ss_use_model: no such model: $1" >&2; return 1; }
  export SS_MODEL="$(cd "$(dirname "$m")" && pwd)/$(basename "$m")"
}
ss_clear_model() { unset SS_MODEL; }

# ss_run_sweep MODEL "size1 size2 ..." "sbatchOpts"
#   Submit the SAME one-step pipeline once per input size (KB), accumulating rows
#   in the shared jobRecord.txt so a curve can be fit. Each size runs in its own
#   workspace subdir (fresh smartSlurmLog + mock state) -- otherwise the step's
#   .success checkpoint from the previous size would make rAP skip resubmission --
#   but the record DB ($smartSlurmJobRecordDir) is NOT wiped between points.
#   The step name (= program key) is taken from the model's MODEL_PROGRAM.
#   Requires that the caller has already created a base workspace (ss_new_workspace),
#   which establishes a clean record DB for this sweep.
ss_run_sweep() {
  local model="$1" sizes="$2" opts="$3"
  ss_use_model "$model" || return 1
  local prog; prog="$(sed -n 's/^MODEL_PROGRAM=//p' "$SS_MODEL" | head -1)"; : "${prog:=sweepJob}"
  local base="${PWD}/_sweep"; mkdir -p "$base"
  local sz d
  for sz in $sizes; do
    d="$base/size_$sz"; rm -rf "$d"; mkdir -p "$d"
    ( cd "$d" || exit 1
      export MOCK_SLURM_STATE="$d/.mockstate"; mkdir -p "$MOCK_SLURM_STATE/jobs"
      truncate -s "${sz}K" in.dat
      printf '#!/bin/sh\ninput=$1\n#@1,0,%s,,input,%s\ncat $input > /dev/null\n\n' "$prog" "$opts" > pipeline.sh
      chmod +x pipeline.sh
      bash "$SMARTSLURM_ROOT/bin/runAsPipeline" \
        --script "./pipeline.sh in.dat" --sbatch-options "$opts" --tmp noTmp \
        >"$d/rap.log" 2>&1
    )
  done
}

# ss_fit PROGRAM REF  : force estimateResource.sh to (re)build the fit from the
#   current jobRecord.txt and echo the fitted mem.stat keys. Returns non-zero if
#   gnuplot is unavailable. Use a size well within the observed range so the
#   builder runs rather than returning outOfRange.
ss_fit() {
  local prog="$1" ref="${2:-none}"
  command -v gnuplot >/dev/null 2>&1 || { echo "gnuplot-missing"; return 2; }
  # pick a mid-range size from the existing records to avoid the outOfRange guard
  local mid
  mid="$(awk -F, -v p="$prog" -v r="$ref" '$12==p&&$13==r&&$2+0>0{print $2}' \
        "$smartSlurmJobRecordDir/jobRecord.txt" | sort -n | awk '{a[NR]=$1} END{print a[int(NR/2)+1]}')"
  : "${mid:=1000}"
  # Force a fresh fit from ALL current records: estimateResource.sh otherwise
  # reuses a <2-min-old stat file (staleness guard), which during a fast sweep
  # can be an early, few-point fit. Removing it makes the builder re-run.
  rm -f "$smartSlurmJobRecordDir/stats/$prog.$ref.mem.stat" \
        "$smartSlurmJobRecordDir/stats/$prog.$ref.time.stat"
  estimateResource.sh "$prog" "$ref" ",jobSize:$mid" fitflag 4096 120 >/dev/null 2>&1
  cat "$smartSlurmJobRecordDir/stats/$prog.$ref.mem.stat" 2>/dev/null
}

# ss_firstpass_rate PROGRAM REF MODEL PCT DIM PATH "sizes" "seeds"
#   Integer % of (size,seed) trials where the estimate at target PCT covers the
#   model's realized need -- i.e. the job would complete on the FIRST try.
#     DIM  = mem | time        which dimension to score
#     PATH = cmt | full        cmt  -> calculateMemTime.sh (the B1 margin alone)
#                              full -> estimateResource.sh (adds extraMem + floor)
#   A non-numeric estimate (e.g. "outOfRange") counts as a miss. Requires the
#   gnuplot stats to be built already (call ss_fit first).
ss_firstpass_rate() {
  local prog="$1" ref="$2" model="$3" pct="$4" dim="$5" path="$6" sizes="$7" seeds="$8"
  local n=0 ok=0 size seed true est tcol ecol
  [ "$dim" = time ] && tcol=2 || tcol=1          # resourceModel prints "mem time rc"
  for size in $sizes; do
    for seed in $seeds; do
      true=$("$HERE/mockSlurm/resourceModel.sh" "$model" "$size" "$seed" | cut -d' ' -f"$tcol")
      if [ "$path" = full ]; then                # estimateResource: "size mem min extra"
        [ "$dim" = time ] && ecol=3 || ecol=2
        est=$(targetFirstPassPct="$pct" estimateResource.sh "$prog" "$ref" ",jobSize:$size" fp 4096 120 2>/dev/null | cut -d' ' -f"$ecol")
      else                                        # calculateMemTime: "mem time ..."
        [ "$dim" = time ] && ecol=2 || ecol=1
        est=$(targetFirstPassPct="$pct" calculateMemTime.sh "$prog" "$ref" "$size" 2>/dev/null | cut -d' ' -f"$ecol")
      fi
      case "$est" in ''|*[!0-9]*) est=0;; esac    # outOfRange / blank -> miss
      case "$true" in ''|*[!0-9]*) true=0;; esac
      n=$((n+1)); [ "$est" -ge "$true" ] && ok=$((ok+1))
    done
  done
  [ "$n" -gt 0 ] && echo $(( ok * 100 / n )) || echo 0
}

# ss_alloc_ladder PROGRAM : the allocated-mem sequence (col5) across a program's
#   jobRecord rows in time order, with status, as "mem:status" tokens. Lets a
#   convergence block show/verify the monotonic escalation ladder.
ss_alloc_ladder() {
  awk -F, -v p="$1" '$12==p{printf "%s:%s ", $5, $9}' \
      "$smartSlurmJobRecordDir/jobRecord.txt" 2>/dev/null
}

# ss_seed_record PROGRAM REF SIZE MEM TIME STATUS [count] [saccMem]
#   Append crafted row(s) directly to jobRecord.txt, bypassing a real run, for
#   testing the readers / fit / pruning against good OR contaminated data
#   deterministically. Production jobRecord.txt is HEADERLESS (cleanUp.sh just
#   appends rows), so this writes no header either — a seeded file is then
#   byte-indistinguishable from a real one. Column layout matches cleanUp.sh's
#   21-field record:
#     jobID,inputSize,memDef,minDef,allocMem,allocTime,memUsed,timeUsed,status,
#     USER,saccMem,program,ref,flag,core,extraMem,extraTime,ratioM,ratioT,date,epoch
ss_seed_record() {
  local prog="$1" ref="$2" size="$3" mem="$4" tim="$5" status="$6" count="${7:-1}" sacc="${8:-$4}"
  local jr="$smartSlurmJobRecordDir/jobRecord.txt"
  mkdir -p "$smartSlurmJobRecordDir"
  local i jid
  for ((i=0; i<count; i++)); do
    jid=$(( (RANDOM % 90000) + 10000 ))
    printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
      "$jid" "$size" "$mem" "$tim" "$mem" "$tim" "$mem" "$tim" "$status" \
      "${USER:-tester}" "$sacc" "$prog" "$ref" "$prog" 1 500 10 "0.50" "0.50" \
      "$(date)" "$(date +%s)" >> "$jr"
  done
}

# assert_near EXPECTED GOT TOL_FRACTION MSG : pass if |got-expected| <= tol*|expected|
assert_near() {
  local exp="$1" got="$2" tol="$3" msg="$4"
  local ok
  ok="$(awk -v e="$exp" -v g="$got" -v t="$tol" 'BEGIN{
    d=e-g; if(d<0)d=-d; lim=(e<0?-e:e)*t; if(lim<1e-9)lim=1e-9; print (d<=lim)?"1":"0"}')"
  [ "$ok" = "1" ] && _pass "$msg" || _fail "$msg" "expected ~$exp (+/-$(awk -v e="$exp" -v t="$tol" 'BEGIN{print e*t}')), got $got"
}

# stat_val KEY FILE : echo the value of KEY= in a stat file (empty if absent)
stat_val() { sed -n "s/^$1=//p" "$2" 2>/dev/null | head -1; }

# run the pipeline. usage: ss_run_pipeline "<script args>" "<sbatch opts>" [--dryrun]
# writes combined output to $LAST_RUN_LOG, returns rAP's exit code.
ss_run_pipeline() {
  local script="$1" opts="$2" mode="${3:-}"
  export LAST_RUN_LOG="${MOCK_SLURM_STATE%/*}/rap.log"
  local mflag=()
  [ "$mode" = "--dryrun" ] && mflag=(--mode dryrun)
  bash "$SMARTSLURM_ROOT/bin/runAsPipeline" \
       --script "$script" --sbatch-options "$opts" --tmp noTmp "${mflag[@]}" \
       >"$LAST_RUN_LOG" 2>&1
}

# --- assertion helpers + grouped PASS/FAIL report -----------------------------
: "${FAILS:=0}" "${CHECKS:=0}"
# Results are accumulated to a temp file so the end-of-run report can group each
# check under its parent block (survives subshells too). Records are TSV:
#   B <key> <purpose>        (one per block, in run order)
#   C <key> PASS|FAIL <name> (one per check)
: "${SS_REPORT:=}"
[ -n "$SS_REPORT" ] || { SS_REPORT="$(mktemp)"; export SS_REPORT; }
CUR_KEY="-"

# _begin <key> [purpose] : start a block; prints its live header and records it.
# If no purpose is given, it is looked up from the suite's ss_purpose() catalog.
_begin() {
  CUR_KEY="$1"
  local purpose="${2:-}"
  [ -n "$purpose" ] || purpose="$(ss_purpose "$1")"
  printf 'B\t%s\t%s\n' "$1" "$purpose" >> "$SS_REPORT"
  echo "== [$1] $purpose =="
}

_pass() { CHECKS=$((CHECKS+1)); printf 'C\t%s\tPASS\t%s\n' "$CUR_KEY" "$1" >> "$SS_REPORT"; echo "    PASS - $1"; }
_fail() { CHECKS=$((CHECKS+1)); FAILS=$((FAILS+1)); printf 'C\t%s\tFAIL\t%s\n' "$CUR_KEY" "$1" >> "$SS_REPORT"; echo "    FAIL - $1"; [ -n "${2:-}" ] && echo "           $2"; }

# ss_report : print the tidy grouped report (blocks, with checks indented).
ss_report() {
  [ -f "$SS_REPORT" ] || return 0
  awk -F'\t' '
    $1=="B" { if (!(($2) in seen)) { order[++n]=$2; seen[$2]=1 } purpose[$2]=$3 }
    $1=="C" { k=$2; tot[k]++; if ($3=="FAIL") fail[k]++
              line[k]=line[k] sprintf("           %-4s  %s\n", $3, $4) }
    END {
      print  ""
      print  "=========================  TEST REPORT  ========================="
      print  ""
      bpass=0; bfail=0
      for (i=1; i<=n; i++) {
        k=order[i]; f=fail[k]+0; t=tot[k]+0; p=t-f
        st=(f>0)?"FAIL":"PASS"
        if (f>0) bfail++; else bpass++
        printf "[%s]  %s  (%s)   %d/%d\n", st, purpose[k], k, p, t
        printf "%s", line[k]
        print ""
      }
      ctot=0; cfail=0
      for (k in tot) { ctot+=tot[k]; cfail+=fail[k] }
      print  "----------------------------------------------------------------"
      printf "Blocks: %d/%d PASS    Checks: %d/%d PASS    RESULT: %s\n", \
             bpass, bpass+bfail, ctot-cfail, ctot, (cfail>0 || bfail>0)?"FAIL":"PASS"
    }' "$SS_REPORT"
}

# ss_catalog <key...> : print "<key> — <purpose>" for --list, without running.
# Purposes come from the suite-level ss_purpose() defined in run_tests.sh.
ss_catalog() {
  local k
  for k in "$@"; do printf '  %-22s  %s\n' "$k" "$(ss_purpose "$k")"; done
}

assert_rc()        { [ "$1" -eq "$2" ] && _pass "$3" || _fail "$3" "expected rc=$2 got rc=$1"; }
assert_file()      { [ -f "$1" ] && _pass "${2:-file exists: $1}" || _fail "${2:-file exists: $1}" "missing: $1"; }
assert_no_file()   { [ ! -e "$1" ] && _pass "${2:-absent: $1}" || _fail "${2:-absent: $1}" "unexpected: $1"; }
assert_nonempty()  { [ -s "$1" ] && _pass "${2:-non-empty: $1}" || _fail "${2:-non-empty: $1}" "empty/missing: $1"; }
assert_empty()     { [ -e "$1" ] && [ ! -s "$1" ] && _pass "${2:-empty: $1}" || _fail "${2:-empty: $1}" "not empty/missing: $1"; }
assert_contains()  { grep -q -- "$2" "$1" 2>/dev/null && _pass "${3:-'$2' in $1}" || _fail "${3:-'$2' in $1}" "not found in $1"; }
assert_valid_bash(){ bash -n "$1" 2>/dev/null && _pass "${2:-valid bash: $1}" || _fail "${2:-valid bash: $1}" "syntax error in $1"; }
assert_invalid_bash(){ bash -n "$1" 2>/dev/null && _fail "${2:-invalid bash: $1}" "unexpectedly valid" || _pass "${2:-invalid bash (known issue): $1}"; }
assert_lacks()     { grep -q -- "$2" "$1" 2>/dev/null && _fail "${3:-'$2' absent from $1}" "unexpectedly found in $1" || _pass "${3:-'$2' absent from $1}"; }

# ss_scn <homedir> : run a shell scenario in a PRISTINE bash and echo its
# combined output. The scenario script is read from stdin. The child shell gets
# a clean environment (env -i) with a sandbox HOME, SMARTSLURM_HOME pointing at
# the checkout under test, and a minimal PATH that contains no 'sbatch' -- so
# 'command -v sbatch' resolving to the SmartSlurm shim is proof of aliasing, and
# resolving to nothing is proof of a clean deactivate. Used by the activation
# blocks, which must exercise PATH/PS1 changes without touching the harness's
# own shell. Caller redirects 2>&1 to capture warnings too.
ss_scn() {
  local home="$1"; mkdir -p "$home"
  env -i HOME="$home" SMARTSLURM_HOME="$SMARTSLURM_ROOT" \
      PATH="/usr/bin:/bin" TERM=dumb \
      bash --noprofile --norc
}
