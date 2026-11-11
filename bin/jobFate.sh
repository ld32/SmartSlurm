#!/bin/bash
# jobFate.sh -- resolve a vanished job's terminal fate from sacct and persist it,
# so checkRunGrid (marker files) and diagnose.sh (the appended job-report region)
# reflect OOM / timeout / cancelled / node-fail instead of a phantom "running" cell
# or an "UNCONFIRMED" diagnose card.
#
# WHY: all of SmartSlurm's status artifacts (.success / .failed / the sacct report
# in <flag>.out) are written by cleanUp.sh from the job's EXIT trap. A SIGKILL
# (cgroup OOM-killer; the hard kill after a timeout grace period; node failure)
# never runs that trap, so none of those artifacts exist -- the only record of the
# job's fate is in Slurm's accounting DB. This helper reads that record.
#
# COST / INVARIANTS: used ONLY on explicit, cold paths -- checkRun's launch/refresh
# reconcile (batched: one sacct call for many flags) and a diagnose drill-in (single
# flag). Never the grid's hot redraw loop. It is read-mostly: it writes only the
# authoritative .failed marker + the sacct report, and only for a job sacct says
# ended in a failure state with no .success flag. Fully graceful -- no jobid, no
# sacct, or no record yet (post-job lag) -> prints UNRESOLVED and writes nothing, so
# the caller keeps showing the live "?" (unknown) until a later refresh resolves it.
#
# Usage:
#   jobFate.sh <flag> [--log DIR]              # single: prints the state token
#   jobFate.sh --batch <flag>... [--log DIR]   # batch: one sacct call; prints "<flag> <state>" per resolved flag
#
# State tokens: OOM OOT NODE_FAIL CANCELLED FAILED COMPLETED UNRESOLVED
# ---------------------------------------------------------------------------
set -uo pipefail

LOGDIR="${smartSlurmLogDir:-smartSlurmLog}"
BATCH=0; FLAGS=()
while [ $# -gt 0 ]; do
  case "$1" in
    --log)          LOGDIR="$2"; shift 2;;
    --batch)        BATCH=1; shift;;
    -h|--help)      sed -n '2,29p' "$0" | sed 's/^# \{0,1\}//'; exit 0;;
    *)              FLAGS+=("$1"); shift;;
  esac
done
[ "${#FLAGS[@]}" -gt 0 ] || { echo "usage: jobFate.sh <flag>... [--log DIR] [--batch]" >&2; exit 2; }

ALLJOBS="$LOGDIR/allJobs.txt"
# Same columns cleanUp.sh emits, so the appended report reads like a real one and
# diagnose.sh's report_state() recognises it unchanged.
SACCT_FMT="JobID,Submit,Start,End,MaxRSS,State,NodeList%30,Partition,ReqTRES%30,TotalCPU,Elapsed%14,Timelimit%14,ExitCode"

have_sacct() { command -v sacct >/dev/null 2>&1; }

# Last jobid recorded for a flag. allJobs.txt is "job_id depend_on job_flag ...";
# on requeue a flag can appear more than once, so the LAST entry wins. Exact string
# compare on $3 (no regex) -> safe for the dotted/underscored real flags.
jobid_of() {
  [ -f "$ALLJOBS" ] || { printf ''; return; }
  local id; id="$(awk -v f="$1" 'NR>1 && $3==f{v=$1} END{printf "%s", v}' "$ALLJOBS")"
  # A real Slurm jobid is all digits. A non-numeric / garbage value (corrupt
  # allJobs.txt) is not a job to query -- return empty so the caller degrades to
  # UNRESOLVED instead of handing sacct a bogus -j argument.
  [[ "$id" =~ ^[0-9]+$ ]] || id=""
  printf '%s' "$id"
}

# Map a sacct blob (a job's alloc line plus its .batch/.0/.extern substeps) to one
# state token. Precedence mirrors cleanUp.sh and diagnose.sh report_state(): a
# memory kill on any substep outranks a bare FAILED/CANCELLED on the batch line.
map_state() {
  # Match whole, case-sensitive State tokens (sacct states are UPPERCASE), not bare
  # substrings: a NodeList/Partition/QOS name containing "cancel" or "failed" must
  # not be read as a job state. OUT_OF_MEMORY/TIMEOUT may appear on a .batch/.0
  # substep while the alloc line shows FAILED, so scan all of this job's lines.
  if   grep -qw 'OUT_OF_MEMORY'        <<<"$1"; then echo OOM
  elif grep -qw 'TIMEOUT'              <<<"$1"; then echo OOT
  elif grep -qwE 'NODE_FAIL|BOOT_FAIL' <<<"$1"; then echo NODE_FAIL
  elif grep -qw 'CANCELLED'            <<<"$1"; then echo CANCELLED
  elif grep -qw 'FAILED'               <<<"$1"; then echo FAILED
  elif grep -qw 'COMPLETED'            <<<"$1"; then echo COMPLETED
  else echo ""; fi
}

# Persist a resolved failure. Idempotent, and authoritative-success always wins.
persist() {
  local flag="$1" state="$2" blob="$3" out="$LOGDIR/$flag.out"
  [ -e "$LOGDIR/$flag.success" ] && return 0            # a real success landed -> never override
  case "$state" in ""|UNRESOLVED) return 0;; esac
  # Append the sacct table once, under the sentinel diagnose.sh scans for, so a
  # drill-in recovers the specific cause via its existing report_state(). Done for
  # any terminal state (incl. COMPLETED) so diagnose can classify the
  # completed-but-no-success case.
  if [ -f "$out" ] && ! grep -q 'SmartSlurm job report:' "$out" 2>/dev/null; then
    { echo "========== SmartSlurm job report: $flag (reconciled from sacct) =========="
      printf '%s\n' "$blob"
    } >> "$out" 2>/dev/null || true
  fi
  # Write the .failed marker (which the grid reads as error) ONLY for an actual
  # tool/resource/infra failure. NOT for COMPLETED (Slurm says it finished -- leave
  # the grid's "?" and let diagnose classify the missing-.success case) and NOT for
  # CANCELLED: a cancelled job is usually an orphan (upstream failed -> afterok
  # cancel) or a manual scancel, not a tool failure, and stamping .failed would mask
  # checkRunGrid's more informative orphan glyph. The appended report still lets
  # diagnose explain the cancellation.
  case "$state" in
    OOM|OOT|NODE_FAIL|FAILED) touch "$LOGDIR/$flag.failed" 2>/dev/null || true;;
  esac
}

resolve_single() {
  local flag="$1" id blob state
  id="$(jobid_of "$flag")"
  [ -n "$id" ] || { echo UNRESOLVED; return 1; }
  have_sacct     || { echo UNRESOLVED; return 1; }
  blob="$(sacct -n -j "$id" --format="$SACCT_FMT" --units=M 2>/dev/null)"
  [ -n "$blob" ]  || { echo UNRESOLVED; return 1; }    # post-job lag: no record yet
  state="$(map_state "$blob")"
  [ -n "$state" ] || { echo UNRESOLVED; return 1; }
  persist "$flag" "$state" "$blob"
  echo "$state"
}

run_batch() {
  local flag id csv blob jblob state
  local -a ids=(); local -A id_of=()
  for flag in "${FLAGS[@]}"; do
    id="$(jobid_of "$flag")"
    [ -n "$id" ] && { ids+=("$id"); id_of["$flag"]="$id"; }
  done
  [ "${#ids[@]}" -gt 0 ] || return 0
  have_sacct              || return 0
  csv="$(IFS=,; echo "${ids[*]}")"
  blob="$(sacct -n -j "$csv" --format="$SACCT_FMT" --units=M 2>/dev/null)"   # ONE call
  [ -n "$blob" ] || return 0
  for flag in "${FLAGS[@]}"; do
    id="${id_of[$flag]:-}"; [ -n "$id" ] || continue
    # this job's lines: JobID (col 1) is exactly <id> or <id>.<step>
    jblob="$(awk -v id="$id" '$1==id || index($1, id".")==1' <<<"$blob")"
    [ -n "$jblob" ] || continue
    state="$(map_state "$jblob")"
    [ -n "$state" ] || continue
    persist "$flag" "$state" "$jblob"
    echo "$flag $state"
  done
}

if [ "$BATCH" = 1 ]; then run_batch; else resolve_single "${FLAGS[0]}"; fi
