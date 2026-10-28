#!/bin/bash
# ---------------------------------------------------------------------------
# run_tests.sh  --  validate that runAsPipeline still works under mock Slurm.
#
#   ./test/run_tests.sh                    run every block
#   ./test/run_tests.sh clean              run blocks whose name matches "clean"
#   ./test/run_tests.sh efficiency_file checkrun_efficiency
#                                          run several named blocks (OR of matches)
#   ./test/run_tests.sh estimator          run a named group (see --groups)
#   ./test/run_tests.sh activation email_consolidation
#                                          groups and names can be mixed freely
#   ./test/run_tests.sh --groups   list the named groups and their blocks
#   ./test/run_tests.sh --list     list every block and what it verifies (no run)
#   ./test/run_tests.sh --help     usage
#
# Each argument is a substring matched against block names (so one token can fan
# out, e.g. "extrapolate" -> 3 blocks) OR a group alias that expands to a block
# list. A block runs if it matches ANY argument; no arguments runs everything.
#
# Ends with a grouped PASS/FAIL report (each block, with its checks indented).
# Exit status is non-zero if any assertion fails (CI-friendly).
# ---------------------------------------------------------------------------
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Every block, in run order, and a one-line functional description of what it
# verifies. This is the single source of truth for block purposes: _begin looks
# the text up here, and --list prints it.
ALL_BLOCKS="clean_2step nested_fanout ragged_axes fail_exit0_empty known_broken_swallowed_done \
loop_detection oom_injection tool_stderr diagnose rerun end_block verbosity scaling_linear \
noise_fit scaling_exp contaminate_exit0 contaminate_exclude contaminate_prune firstpass_rate \
failure_floor manual_rerun_floor constant_model negslope_robust \
setup_rc_block activate_path activate_prompt deactivate_surgical conda_coexist \
uninstall_precise unexport_deprecated \
firstpass_fullpath extrapolate_linear extrapolate_superlinear extrapolate_lowdata \
resubmit_convergence manual_resubmit_sequence \
efficiency_file checkrun_efficiency email_consolidation"

# Named groups: a convenience alias that expands to a block list, so related
# blocks that do NOT share a common name substring can be run in one go (e.g.
# the activation feature spans activate/deactivate/uninstall/conda/setup). A
# non-group argument is returned unchanged and matched as a substring. Keep the
# membership lists here in sync with ALL_BLOCKS above.
GROUP_NAMES="estimator activation efficiency"
ss_group() {
  case "$1" in
    estimator)  echo "scaling_linear noise_fit scaling_exp contaminate_exit0 contaminate_exclude contaminate_prune firstpass_rate failure_floor manual_rerun_floor constant_model negslope_robust firstpass_fullpath extrapolate_linear extrapolate_superlinear extrapolate_lowdata resubmit_convergence manual_resubmit_sequence";;
    activation) echo "setup_rc_block activate_path activate_prompt deactivate_surgical conda_coexist uninstall_precise unexport_deprecated";;
    efficiency) echo "efficiency_file checkrun_efficiency email_consolidation";;
    *)          echo "$1";;
  esac
}

ss_purpose() {
  case "$1" in
    clean_2step)        echo "Baseline: convert a 2-step pipeline, submit under mock Slurm, and wire afterok dependencies end-to-end";;
    nested_fanout)      echo "Nested per-sample/per-chrom loops — a fan-out step (one grid cell backs multiple jobs) converts to valid bash and every job completes";;
    ragged_axes)        echo "Two steps whose loop axes are disjoint (sample ids vs group ids) — ragged addressing converts and runs end-to-end";;
    loop_detection)     echo "Documents the textual for*/done* loop foot-gun: a one-line for...done and a for-prefixed command corrupt the next step's flag";;
    fail_exit0_empty)   echo "Producer exits 0 but writes nothing and the consumer fails — documents the current symptom-vs-cause misattribution";;
    known_broken_swallowed_done) echo "A blank-line block right before a loop 'done' makes the converter emit invalid bash (the M4 block-closure foot-gun)";;
    oom_injection)      echo "Injected OOM/timeout/nonzero jobs are attributed correctly by cleanUp (no spurious .success)";;
    tool_stderr)        echo "Tool chatter on stderr is bracketed by sentinels, not misread as a failure";;
    diagnose)           echo "diagnose.sh classifies failure state/cause/action and emits the expected JSON contract";;
    rerun)              echo "Rerun menu (failed / all / from-step:N) materializes markers and honors the reachable guard";;
    end_block)          echo "'#@ … #@end' scripts-on-disk: closure validation, --lint, legacy back-compat, and in-block multivar escaping";;
    verbosity)          echo "Stdout verbosity: default A+ per-sample rows vs -v firehose vs -q silence";;
    scaling_linear)     echo "Linear mem/time vs input size — the real gnuplot fit recovers slope/intercept and prediction tracks the line";;
    noise_fit)          echo "~8% stochastic noise — the fit still recovers the slope; STDFIT>0 and RSquare are computed";;
    scaling_exp)        echo "Exponential truth vs SmartSlurm's linear fit — the model mismatch under-provisions at large input sizes";;
    contaminate_exit0)  echo "An exit-0-empty / mis-sized COMPLETED row drags the fitted slope downward (data contamination)";;
    contaminate_exclude) echo "OOM/failed rows are excluded from the fit — only COMPLETED rows are used";;
    contaminate_prune)  echo "An OOM storm prunes the offending program's failed and failure-disproven rows, keeping its validated history and unrelated programs";;
    firstpass_rate)     echo "Measured first-pass success rate tracks the targetFirstPassPct knob (and responds monotonically to it)";;
    failure_floor)      echo "Estimates never drop to/below a proven-insufficient allocation; a genuinely smaller input is unaffected";;
    manual_rerun_floor) echo "After a real emergent-OOM run, a fresh (manual) rerun estimate exceeds the proven-insufficient allocation (cross-run)";;
    constant_model)     echo "Constant memory (size-independent) with scaling time — mem uses the constant model, not a negative-slope collapse to the floor";;
    negslope_robust)    echo "Numeric edge cases (negative/sub-1 results, empty STDFIT) always yield a valid floored integer, never bc '-' garbage";;
    setup_rc_block)     echo "smartSlurmSetup writes exactly one marked initialize block to the rc file, is idempotent on re-run, backs up, and preserves the user's own lines";;
    activate_path)      echo "'smartslurm activate' puts bin/ and sbatchBin/ on PATH, aliases sbatch -> ssbatch, and announces the aliasing (the previously-silent behavior)";;
    activate_prompt)    echo "'smartslurm activate' prefixes the prompt with '{smartSlurm} ' (curly, not conda's parens) exactly once even on repeat activate";;
    deactivate_surgical) echo "'smartslurm deactivate' removes ONLY smartSlurm's own PATH/prompt additions, preserving anything added after activation";;
    conda_coexist)      echo "smartSlurm and a conda-style '(env) ' prompt + PATH entry coexist: activating/deactivating either in any order never clobbers the other";;
    uninstall_precise)  echo "'smartslurm uninstall' deletes exactly the initialize block from the rc file, leaving surrounding user lines byte-identical, with a backup";;
    unexport_deprecated) echo "the deprecated 'unExport' warns and delegates to 'smartslurm deactivate' (no longer strips the system sbatch dir)";;
    firstpass_fullpath) echo "First-pass coverage through the FULL estimateResource path (extraMem + failure floor), for mem AND time (not just the calculateMemTime margin)";;
    extrapolate_linear) echo "[Phase 1 baseline] Linear truth trained on small sizes, predicted at 2x-50x: with SCount>9 it extrapolates the line (not outOfRange); reports first-pass per scale";;
    extrapolate_superlinear) echo "[Phase 1 baseline] Quadratic truth trained on small sizes: the straight-line fit UNDER-provisions at 10x-50x; records the shortfall ratio (the problem to fix)";;
    extrapolate_lowdata) echo "Thin data (<10 records) beyond the trained range uses a trend-scaled fallback (in-range estimate x size/Maximum), not a bare default";;
    resubmit_convergence) echo "Auto-resubmit: an under-provisioned job OOMs and is auto-requeued with monotonically escalating resources until it completes; verifies the ladder never regresses";;
    manual_resubmit_sequence) echo "Multiple sequential MANUAL resubmissions (each with its own auto-requeue), bootstrapping to completion; each manual start exceeds the prior proven-insufficient alloc";;
    efficiency_file)    echo "cleanUp writes one unconditional per-job row to a per-run efficiency.tsv (flag, sample, status, alloc/used mem+time) for checkRun to read";;
    checkrun_efficiency) echo "checkRunGrid --efficiency drill-down views (overview/sample/step): Overall + per-sample bars, per-step and per-sample breakdowns, eff/used/alloc columns, time/mem toggles, color";;
    email_consolidation) echo "cleanUp sends consolidated TEXT-ONLY mail (no pngs): per-sample break/complete alerts + one per-run efficiency email; no per-job spam, no hard-coded dev recipient";;
    *)                  echo "$1";;
  esac
}

source "$HERE/lib/harness.sh"

# --list / --groups / --help short-circuit before any setup
case "${1:-}" in
  --list)
    echo "Test blocks (name — what it verifies):"; echo
    ss_catalog $ALL_BLOCKS
    exit 0;;
  --groups)
    echo "Named groups (alias — blocks it expands to):"; echo
    for g in $GROUP_NAMES; do printf '  %-12s %s\n' "$g" "$(ss_group "$g")"; done
    exit 0;;
  --help|-h)
    sed -n '2,22p' "$0" | sed 's/^# \{0,1\}//'
    exit 0;;
esac

ss_activate

# Build the filter list from every argument: group aliases expand to their block
# lists, any other token is kept as a literal substring. No arguments => empty
# filter list => every block runs.
FILTERS=""
for a in "$@"; do FILTERS="$FILTERS $(ss_group "$a")"; done
WORK="${SS_SCRATCH:-/tmp/ss-scratch}/work"
FAILS=0; CHECKS=0

# A block runs when the filter list is empty, or its name contains ANY filter.
_want() {
  [ -z "${FILTERS// /}" ] && return 0
  local b="$1" f
  for f in $FILTERS; do [[ "$b" == *"$f"* ]] && return 0; done
  return 1
}

# ===========================================================================
_want clean_2step && {
  _begin clean_2step
  ss_new_workspace "$WORK/clean_2step"
  cp "$HERE"/fixtures/clean_2step/* .
  ss_run_pipeline "./pipeline.sh 777" "sbatch -p short -c 1 --mem 4G -t 50:0"; rc=$?
  assert_rc "$rc" 0 "runAsPipeline exits 0"
  driver=$(ls smartSlurmLog/slurmPipeLine.*.sh 2>/dev/null | grep -v dryrun | head -1)
  assert_file    "$driver"                                   "converted driver produced"
  assert_valid_bash "$driver"                                "driver is valid bash"
  assert_file    "smartSlurmLog/1.0.findNumber.1.success"    "step 1a marked success"
  assert_file    "smartSlurmLog/1.0.findNumber.2.success"    "step 1b marked success"
  assert_file    "smartSlurmLog/2.1.mergeNumber.success"     "step 2 marked success"
  assert_no_file "smartSlurmLog/2.1.mergeNumber.failed"      "step 2 not marked failed"
  assert_nonempty "777.1.txt"                                "payload 1a wrote output"
  assert_nonempty "all777.txt"                               "merge produced combined output"
  assert_contains "all777.txt" "777"                         "merged output contains data"
  assert_contains "smartSlurmLog/allJobs.txt" "mergeNumber"  "allJobs.txt records the DAG"
  # dependency wiring: merge job lists both producers as afterok deps
  assert_contains "smartSlurmLog/allJobs.txt" ":"            "merge depends on upstream job ids"
  echo
}

# ===========================================================================
# Nested loops produce a FAN-OUT step: one (sample,step) grid cell backs multiple
# jobs (per-chromosome). Proves rAP converts nested #@ loops to valid bash and
# every fanned-out job runs — the structural case behind checkRunGrid's '@' cells.
_want nested_fanout && {
  _begin nested_fanout
  ss_new_workspace "$WORK/nested_fanout"
  cp "$HERE"/fixtures/nested_fanout/* .; chmod +x ./*.sh 2>/dev/null || true
  ss_run_pipeline "./pipeline.sh 777" "sbatch -p short -c 1 --mem 2G -t 50:0"; rc=$?
  assert_rc "$rc" 0 "runAsPipeline exits 0"
  driver=$(ls smartSlurmLog/slurmPipeLine.*.sh 2>/dev/null | grep -v dryrun | head -1)
  assert_valid_bash "$driver" "nested-loop driver is valid bash"
  assert_file "smartSlurmLog/1.0.prep.1.success"          "prep (sample 1) succeeded"
  assert_file "smartSlurmLog/1.0.prep.2.success"          "prep (sample 2) succeeded"
  # fan-out: step 2 is one grid cell (2.1.callChrom) backing 2 samples x 3 chroms
  n=$(ls smartSlurmLog/2.1.callChrom.*.success 2>/dev/null | wc -l)
  assert_rc 6 "$n" "fan-out produced 6 callChrom jobs (2 samples x 3 chroms)"
  assert_file "smartSlurmLog/3.2.summarize.success"       "summarize (downstream of the fan-out) succeeded"
  assert_contains "smartSlurmLog/allJobs.txt" "callChrom"  "allJobs.txt records the fan-out step"
  assert_nonempty "summary777.txt"                        "summarize combined the fanned-out outputs"
  echo
}

# ===========================================================================
# Two steps whose loop axes are DISJOINT (sample ids {1,2} vs group ids {7,8}).
# Proves rAP addresses ragged axes without cross-contaminating the flags.
_want ragged_axes && {
  _begin ragged_axes
  ss_new_workspace "$WORK/ragged_axes"
  cp "$HERE"/fixtures/ragged_axes/* .; chmod +x ./*.sh 2>/dev/null || true
  ss_run_pipeline "./pipeline.sh 777" "sbatch -p short -c 1 --mem 2G -t 50:0"; rc=$?
  assert_rc "$rc" 0 "runAsPipeline exits 0"
  driver=$(ls smartSlurmLog/slurmPipeLine.*.sh 2>/dev/null | grep -v dryrun | head -1)
  assert_valid_bash "$driver" "ragged-axis driver is valid bash"
  assert_file "smartSlurmLog/1.0.align.1.success"      "align axis: sample 1 succeeded"
  assert_file "smartSlurmLog/1.0.align.2.success"      "align axis: sample 2 succeeded"
  assert_file "smartSlurmLog/2.1.perGroup.7.success"   "group axis: group 7 succeeded (disjoint id)"
  assert_file "smartSlurmLog/2.1.perGroup.8.success"   "group axis: group 8 succeeded (disjoint id)"
  assert_nonempty "group7.txt"                         "group-axis step produced output"
  echo
}

# ===========================================================================
_want fail_exit0_empty && {
  _begin fail_exit0_empty
  ss_new_workspace "$WORK/fail_exit0_empty"
  cp "$HERE"/fixtures/fail_exit0_empty/* .
  ss_run_pipeline "./pipeline.sh 42" "sbatch -p short -c 1 --mem 4G -t 50:0"; rc=$?
  assert_rc "$rc" 0 "runAsPipeline exits 0 (submission succeeds)"
  assert_empty   "all42.txt"                              "producer wrote an EMPTY output"
  # CURRENT (pre-overhaul) behavior: producer exited 0 -> marked success (the bug)
  assert_file    "smartSlurmLog/1.0.makeTable.success"    "[baseline] producer WRONGLY marked success"
  # the symptom step is the one that actually fails
  assert_file    "smartSlurmLog/2.1.consume.failed"       "[baseline] consumer (symptom) marked failed"
  assert_contains "smartSlurmLog/2.1.consume.out" "empty" "consumer error mentions empty input"
  echo "    note: after M1+M2 this fixture flips — producer FAILED, consumer a symptom."
  echo
}

# ===========================================================================
_want known_broken_swallowed_done && {
  _begin known_broken_swallowed_done
  ss_new_workspace "$WORK/known_broken"
  cp "$HERE"/fixtures/known_broken_swallowed_done/* .
  chmod +x ./*.sh 2>/dev/null || true
  ss_run_pipeline "./pipeline.sh 5" "sbatch -p short -c 1 --mem 4G -t 50:0" --dryrun || true
  driver=$(ls smartSlurmLog/slurmPipeLine.*.dryrun.sh 2>/dev/null | head -1)
  assert_file "$driver" "driver was emitted"
  # the outer `done` is swallowed into the step-2 --wrap; driver is invalid bash
  assert_invalid_bash "$driver" "driver is INVALID bash (swallowed 'done' — known issue)"
  assert_contains "$driver" 'done;' "the 'done' was absorbed into a wrapped command"
  echo "    note: M4's #@ … #@end closure fixes this; fixture is the regression anchor."
  echo
}

_want loop_detection && {
  _begin loop_detection
  ss_new_workspace "$WORK/loop_detection"
  cp "$HERE"/fixtures/loop_detection/* .
  chmod +x ./*.sh 2>/dev/null || true
  ss_run_pipeline "./pipeline.sh input" "sbatch -p short -c 1 --mem 100M -t 5:0" --dryrun || true
  driver=$(ls smartSlurmLog/slurmPipeLine.*.dryrun.sh 2>/dev/null | head -1)
  assert_file "$driver" "driver was emitted"
  # control: the step before any trap carries a clean flag
  assert_contains "$driver" 'flag: 1\.0\.first' "step before the loop trap is clean (1.0.first)"
  # Once a one-line `for x ...; done` used to push .$x and never pop it, corrupting
  # the NEXT step's flag with a stray .$x loop suffix. rAP loop detection is now
  # hardened, so the next step's flag must stay clean.
  assert_contains "$driver" 'flag: 2\.1\.second'  "one-line for...done no longer corrupts the next flag (2.1.second)"
  assert_lacks    "$driver" 'flag: 2\.1\.second\.\$x' "no stray .\$x loop suffix leaks onto 2.1.second"
  # `format_report ...` (first word starts with "for") used to be misread as a loop;
  # its 2nd token (--nope) became a bogus loop var. It must now read clean too.
  assert_contains "$driver" 'flag: 3\.2\.third'   "for-prefixed command no longer corrupts the flag (3.2.third)"
  assert_lacks    "$driver" 'flag: 3\.2\.third\.\$' "no bogus loop-var suffix leaks onto 3.2.third"
  echo "    note: regression guard — if these fail, rAP loop detection (word-boundary"
  echo "    match + loop-var identifier validation + one-line guard + balance check) broke."
  echo
}

# ===========================================================================
_want oom_injection && {
  _begin oom_injection
  ss_new_workspace "$WORK/oom_injection"
  cp "$HERE"/fixtures/clean_2step/* .
  chmod +x ./*.sh 2>/dev/null || true
  MOCK_FAIL_JOB=mergeNumber MOCK_FAIL_MODE=oom \
    ss_run_pipeline "./pipeline.sh 777" "sbatch -p short -c 1 --mem 4G -t 50:0"; rc=$?
  assert_rc "$rc" 0 "runAsPipeline exits 0 (submission still succeeds)"
  assert_file    "smartSlurmLog/1.0.findNumber.1.success"      "unrelated producer still succeeds"
  assert_no_file "smartSlurmLog/2.1.mergeNumber.success"       "OOM job did NOT touch .success"
  assert_contains "smartSlurmLog/2.1.mergeNumber.out" "Out Of Memory" "OOM message present in job log"
  assert_contains "smartSlurmLog/2.1.mergeNumber.out" "jobStatus: OOM" "cleanUp classifies job as OOM via sacct"
  echo
}

# ===========================================================================
_want tool_stderr && {
  _begin tool_stderr
  ss_new_workspace "$WORK/tool_stderr"
  cp "$HERE"/fixtures/tool_stderr/* .
  chmod +x ./*.sh 2>/dev/null || true
  ss_run_pipeline "./pipeline.sh 777" "sbatch -p short -c 1 --mem 4G -t 50:0"; rc=$?
  out=smartSlurmLog/1.0.noisyTool.out
  assert_file    "$out"                                    "job .out exists"
  assert_contains "$out" "SmartSlurm setup:"               "setup sentinel present"
  assert_contains "$out" "SmartSlurm step output:"         "step-output sentinel present"
  assert_contains "$out" "SmartSlurm job report:"          "job-report sentinel present"
  # the tool region = lines between the step-output and job-report sentinels
  awk '/SmartSlurm step output:/{f=1;next} /SmartSlurm job report:/{f=0} f' "$out" > region.txt
  assert_contains region.txt "ERROR: no usable records"   "killing error sits inside the tool region"
  assert_contains region.txt "WARNING: 3 records skipped"  "tool stderr sits inside the tool region"
  if grep -q "memCpuMonitor" region.txt; then _fail "monitor noise leaked into tool region"; else _pass "monitor noise kept out of tool region"; fi
  if grep -qE "Cleanup counter|gnuplot|Job summary" region.txt; then _fail "teardown noise leaked into tool region"; else _pass "teardown noise kept out of tool region"; fi
  echo
}

# ===========================================================================
_want diagnose && {
  _begin diagnose
  # helper: assert a JSON field value from diagnose.sh
  dfield() { diagnose.sh "$2" --log "$1" --json | grep -o "\"$3\":\"\?[^\",}]*" | head -1 | sed 's/.*"'"$3"'":"\?//'; }
  adiag() { local got; got="$(dfield "$1" "$2" "$3")"; [ "$got" = "$4" ] && _pass "diag $2 → $3=$4" || _fail "diag $2 → $3=$4" "got '$got'"; }

  ss_new_workspace "$WORK/d_clean"; cp "$HERE"/fixtures/clean_2step/* .; chmod +x ./*.sh
  ss_run_pipeline "./pipeline.sh 777" "sbatch -p short -c 1 --mem 4G -t 50:0"
  adiag smartSlurmLog 1.0.findNumber.1 failed false
  adiag smartSlurmLog 1.0.findNumber.1 action NONE

  ss_new_workspace "$WORK/d_tool"; cp "$HERE"/fixtures/tool_stderr/* .; chmod +x ./*.sh
  ss_run_pipeline "./pipeline.sh 777" "sbatch -p short -c 1 --mem 4G -t 50:0"
  adiag smartSlurmLog 1.0.noisyTool state COMPLETED_NOSUCCESS
  if diagnose.sh 1.0.noisyTool --log smartSlurmLog --json | grep -q 'ERROR: no usable records'; then _pass "diag noisyTool evidence = tool error line"; else _fail "diag noisyTool evidence = tool error line"; fi

  ss_new_workspace "$WORK/d_cmd"; cp "$HERE"/fixtures/tool_cmd_notfound/* .; chmod +x ./*.sh 2>/dev/null
  ss_run_pipeline "./pipeline.sh 777" "sbatch -p short -c 1 --mem 4G -t 50:0"
  adiag smartSlurmLog 1.0.badcmd category environment:missing-command
  adiag smartSlurmLog 1.0.badcmd action DOOMED

  ss_new_workspace "$WORK/d_segf"; cp "$HERE"/fixtures/tool_segfault/* .; chmod +x ./*.sh
  ss_run_pipeline "./pipeline.sh 777" "sbatch -p short -c 1 --mem 4G -t 50:0"
  adiag smartSlurmLog 1.0.crasher category tool:crash
  adiag smartSlurmLog 1.0.crasher action DOOMED

  ss_new_workspace "$WORK/d_oom"; cp "$HERE"/fixtures/clean_2step/* .; chmod +x ./*.sh
  MOCK_FAIL_JOB=findNumber.1 MOCK_FAIL_MODE=oom ss_run_pipeline "./pipeline.sh 777" "sbatch -p short -c 1 --mem 4G -t 50:0"
  adiag smartSlurmLog 1.0.findNumber.1 state OOM
  adiag smartSlurmLog 1.0.findNumber.1 retryable true
  adiag smartSlurmLog 2.1.mergeNumber category orphaned:upstream-failed
  adiag smartSlurmLog 2.1.mergeNumber action CHECK_UPSTREAM
  echo
}

# ===========================================================================
_want rerun && {
  _begin rerun
  reissue() { smartSlurmRerun="$1" bash "$SMARTSLURM_ROOT/bin/runAsPipeline" --script "./pipeline.sh 9" \
      --sbatch-options "sbatch -p short -c 1 --mem 2G -t 50:0" --tmp noTmp </dev/null >"$2" 2>&1; }

  # --- failed: keep successful, rerun failed + orphaned downstream ---
  ss_new_workspace "$WORK/rerun_failed"
  cp "$HERE"/fixtures/three_step/* .; chmod +x ./*.sh
  MOCK_FAIL_JOB=count.1 MOCK_FAIL_MODE=oom ss_run_pipeline "./pipeline.sh 9" "sbatch -p short -c 1 --mem 2G -t 50:0"
  assert_no_file "smartSlurmLog/2.1.count.1.success" "run1: count.1 failed (OOM)"
  assert_file    "smartSlurmLog/1.0.trim.1.success"  "run1: trim.1 succeeded"
  t_keep=$(stat -c '%Y' smartSlurmLog/1.0.trim.1.success); sleep 1.1
  reissue failed reissue.log
  assert_rc $? 0 "failed: re-issue exits 0 (non-interactive)"
  assert_file "smartSlurmLog/2.1.count.1.success" "failed: count.1 reran and succeeded"
  assert_file "smartSlurmLog/3.2.merge.success"   "failed: orphaned merge reran and succeeded"
  [ "$(stat -c '%Y' smartSlurmLog/1.0.trim.1.success)" = "$t_keep" ] && _pass "failed: successful trim.1 was KEPT (not rerun)" || _fail "failed: successful trim.1 was KEPT (not rerun)"

  # --- from-step: keep step1, rerun step2 onward ---
  ss_new_workspace "$WORK/rerun_fromstep"
  cp "$HERE"/fixtures/three_step/* .; chmod +x ./*.sh
  ss_run_pipeline "./pipeline.sh 9" "sbatch -p short -c 1 --mem 2G -t 50:0"
  s_trim=$(stat -c '%Y' smartSlurmLog/1.0.trim.1.success); s_cnt=$(stat -c '%Y' smartSlurmLog/2.1.count.1.success); sleep 1.1
  reissue fromStep:2 fs.log
  [ "$(stat -c '%Y' smartSlurmLog/1.0.trim.1.success)" = "$s_trim" ] && _pass "from-step2: step1 trim.1 KEPT" || _fail "from-step2: step1 trim.1 KEPT"
  [ "$(stat -c '%Y' smartSlurmLog/2.1.count.1.success)" != "$s_cnt" ] && _pass "from-step2: step2 count.1 RERAN" || _fail "from-step2: step2 count.1 RERAN"
  assert_file "smartSlurmLog/reRunAllSuccessJobs2.txt" "from-step2: step2 marked rerun"
  assert_file "smartSlurmLog/skipAllSuccessJobs1.txt"  "from-step2: step1 marked keep"

  # --- reachable guard: cannot pick a step the run never reached ---
  reissue fromStep:9 guard.log
  assert_contains "guard.log" "not in the completed set" "guard: unreached step 9 rejected"
  echo
}

# ===========================================================================
_want end_block && {
  _begin end_block
  # (1) a real #@end pipeline: converts to valid bash, runs, writes per-job scripts
  ss_new_workspace "$WORK/m4_run"
  cp "$HERE"/fixtures/end_block/* .; chmod +x ./*.sh
  ss_run_pipeline "./pipeline.sh 777" "sbatch -p short -c 1 --mem 4G -t 50:0"
  assert_file "smartSlurmLog/1.0.analyze.1.success" "analyze.1 succeeded"
  assert_file "smartSlurmLog/2.1.merge.success"     "merge succeeded"
  assert_file "smartSlurmLog/jobs/1.0.analyze.1/cmd.sh" "per-job script persisted on disk"
  # submit-scope var baked ($s->1, $number->777), in-block var left literal ($label)
  assert_contains "smartSlurmLog/jobs/1.0.analyze.1/cmd.sh" 'label="sample 1"' "submit-scope \$s baked per job"
  assert_contains "smartSlurmLog/jobs/1.0.analyze.1/cmd.sh" 'processed $label'  "in-block \$label left for the node"
  assert_contains "777.1.log" "processed sample 1" "in-block var expanded on the node at run time"

  # (1b) regression guard: MULTIPLE in-block vars must each be escaped (the
  # conversion runs with IFS='', which once broke the var-splitting so only a
  # single-var block escaped correctly). report step assigns a,b,c then uses them.
  ss_new_workspace "$WORK/m4_multivar"
  printf '%s\n' '#!/bin/sh' 'for s in 1 2; do' '  in=in$s.txt' '' \
    '  #@1,0,multi,,in,sbatch -p short -c 1 --mem 100M -t 5:00' \
    '  a=$(wc -w < $in)' '  b=$(head -1 $in)' '  c=$(wc -c < $in)' \
    '  printf "%s w=%s h=%s c=%s\n" "$s" "$a" "$b" "$c" > out$s.txt' \
    '  #@end' '' 'done' > pipeline.sh
  printf 'hello world\n' > in1.txt; printf 'foo bar baz\n' > in2.txt; chmod +x pipeline.sh
  ss_run_pipeline "./pipeline.sh" "sbatch -p short -c 1 --mem 100M -t 5:00"
  assert_file "smartSlurmLog/1.0.multi.1.success" "multivar: job succeeded"
  assert_contains "smartSlurmLog/jobs/1.0.multi.1/cmd.sh" 'w=%s h=%s c=%s' "multivar: printf format intact"
  # all three in-block vars must survive as literal $a/$b/$c in the per-job script
  if grep -q 'w=\$a h=\$b c=\$c' "smartSlurmLog/jobs/1.0.multi.1/cmd.sh" || grep -qE '"\$a" "\$b" "\$c"' "smartSlurmLog/jobs/1.0.multi.1/cmd.sh"; then _pass "multivar: all in-block vars left for the node"; else _fail "multivar: all in-block vars left for the node"; fi
  assert_contains "out1.txt" "w=2 h=hello world c=12" "multivar: node-side expansion correct"


  ss_new_workspace "$WORK/m4_legacy"
  cp "$HERE"/fixtures/clean_2step/* .; chmod +x ./*.sh
  ss_run_pipeline "./pipeline.sh 777" "sbatch -p short -c 1 --mem 4G -t 50:0"
  assert_file "smartSlurmLog/1.0.findNumber.1.success" "legacy blank-line block still runs"
  [ ! -d smartSlurmLog/jobs ] && _pass "legacy block used --wrap (no jobs/ dir)" || _fail "legacy block used --wrap (no jobs/ dir)"

  # (3) closure validations (dryrun; must emit Error and not build a runnable driver)
  _m4bad() { ss_new_workspace "$WORK/m4_$1" >/dev/null 2>&1; printf '%s\n' "$2" > pipeline.sh
    printf '777\n' > numbers1.txt; chmod +x pipeline.sh
    ss_run_pipeline "./pipeline.sh 777" "sbatch -p short -c 1 --mem 2G -t 50:0" --dryrun
    if grep -qi 'Error:' "$LAST_RUN_LOG"; then _pass "validation: $1 rejected"; else _fail "validation: $1 rejected"; fi; }
  _m4bad orphan_end  $'#!/bin/sh\n#@1,0,a,,,sbatch -p short\necho hi\n#@end\n#@end'
  _m4bad nested_open $'#!/bin/sh\n#@1,0,a,,,sbatch -p short\necho hi\n#@2,0,b,,,sbatch -p short\necho bye\n#@end'
  _m4bad unclosed    $'#!/bin/sh\n#@1,0,a,,,sbatch -p short\necho hi\n#@end\n#@2,0,b,,,sbatch -p short\necho bye'

  # (4) --lint exit codes
  ss_new_workspace "$WORK/m4_lintok"; cp "$HERE"/fixtures/end_block/* .; chmod +x ./*.sh
  bash "$SMARTSLURM_ROOT/bin/runAsPipeline" --script "./pipeline.sh 777" --sbatch-options "sbatch -p short -c 1 --mem 4G -t 50:0" --tmp noTmp --lint </dev/null >/dev/null 2>&1
  assert_rc $? 0 "--lint passes a good #@end pipeline (exit 0)"
  ss_new_workspace "$WORK/m4_lintbad"; printf '%s\n' $'#!/bin/sh\n#@1,0,a,,,sbatch -p short\necho hi\n#@2,0,b,,,sbatch -p short\necho bye\n#@end' > pipeline.sh; printf '777\n' > numbers1.txt; chmod +x pipeline.sh
  bash "$SMARTSLURM_ROOT/bin/runAsPipeline" --script "./pipeline.sh 777" --sbatch-options "sbatch -p short -c 1 --mem 4G -t 50:0" --tmp noTmp --lint </dev/null >/dev/null 2>&1
  assert_rc $? 1 "--lint fails a malformed pipeline (exit 1)"
  echo
}

_want verbosity && {
  _begin verbosity
  vrun() { ss_new_workspace "$WORK/vb_$1" >/dev/null 2>&1; cp "$HERE"/fixtures/three_step/* .; chmod +x ./*.sh
    bash "$SMARTSLURM_ROOT/bin/runAsPipeline" --script "./pipeline.sh 9" --sbatch-options "sbatch -p short -c 1 --mem 2G -t 50:0" --tmp noTmp $2 </dev/null 2>&1 | sed 's/\x1b\[[0-9;]*m//g' > "$3"; }

  vrun normal ""   n.txt
  assert_contains  n.txt "sample .*steps"       "normal: shows sample/steps header"
  assert_contains  n.txt "merge .*3✓"           "normal: shows a per-sample row"
  assert_contains  n.txt "1✓ 2✓"               "normal: appends per-step glyphs"
  assert_contains  n.txt "Submitted 5, kept 0" "normal: one-line summary"
  if grep -q "Working on inputs" n.txt; then _fail "normal: firehose suppressed"; else _pass "normal: firehose suppressed"; fi

  vrun verbose "-v" v.txt
  assert_contains  v.txt "Working on inputs"      "verbose: per-job firehose restored"
  assert_contains  v.txt "All submitted jobs"     "verbose: full job table on stdout"
  if grep -q "1✓ 2✓" v.txt; then _fail "verbose: no A+ rows"; else _pass "verbose: no A+ rows"; fi

  vrun quiet "-q" q.txt
  if grep -qE "Submitting jobs|1✓|Working on inputs" q.txt; then _fail "quiet: silent on success"; else _pass "quiet: silent on success"; fi
  echo
}

# ===========================================================================
# M1: resource capture / curve fitting / jobRecords (resource-model layer)
# ===========================================================================

_want scaling_linear && {
  _begin scaling_linear
  if ! command -v gnuplot >/dev/null 2>&1; then
    echo "    (skipped: gnuplot unavailable)"; echo
  else
  ss_new_workspace "$WORK/scaling_linear" >/dev/null
  ss_run_sweep scaling_linear "1000 2000 3000 4000 6000 8000" "sbatch -p short -c 1 --mem 16G -t 600:0"
  jr="$smartSlurmJobRecordDir/jobRecord.txt"
  n=$(awk -F, '$9=="COMPLETED" && $12=="scaleLin"{c++} END{print c+0}' "$jr")
  assert_rc 6 "$n" "six COMPLETED records, one per size"
  # size 4000 recorded exactly per model: mem=0.5*4000+200=2200, time=0.02*4000+10=90
  assert_contains "$jr" "4000,.*,2200,90,COMPLETED,.*,scaleLin,none" "size 4000 row matches model exactly"
  ss_fit scaleLin none >/dev/null
  mstat="$smartSlurmJobRecordDir/stats/scaleLin.none.mem.stat"
  assert_near 0.5 "$(stat_val Finala "$mstat")" 0.05 "real gnuplot fit: mem slope ~0.5"
  assert_near 200 "$(stat_val Finalb "$mstat")" 0.10 "real gnuplot fit: mem intercept ~200"
  assert_rc 6 "$(stat_val SCount "$mstat")" "fit used all six points"
  # prediction at a new in-range size tracks truth + margin: 0.5*5000+200 +500extra = 3200
  pm=$(estimateResource.sh scaleLin none ,jobSize:5000 p 4096 120 2>/dev/null | cut -d' ' -f2)
  assert_near 3200 "$pm" 0.10 "prediction at size 5000 tracks the fitted line (+margin)"
  echo
  fi
}

# ===========================================================================
_want noise_fit && {
  _begin noise_fit
  if ! command -v gnuplot >/dev/null 2>&1; then
    echo "    (skipped: gnuplot unavailable)"; echo
  else
  ss_new_workspace "$WORK/noise_fit" >/dev/null
  ss_run_sweep scaling_noise "1000 2000 3000 4000 5000 6000 7000 8000 9000 10000" \
               "sbatch -p short -c 1 --mem 32G -t 900:0"
  jr="$smartSlurmJobRecordDir/jobRecord.txt"
  # noise really perturbs the data: the size-1000 row is NOT the exact 700
  exact1000=$(awk -F, '$2==1000 && $12=="scaleLin"{print $7; exit}' "$jr")
  if [ "$exact1000" != "700" ]; then _pass "noise perturbs recorded mem (size 1000: $exact1000 != exact 700)"
  else _fail "noise perturbs recorded mem" "size-1000 mem is exactly 700 (no noise?)"; fi
  ss_fit scaleLin none >/dev/null
  mstat="$smartSlurmJobRecordDir/stats/scaleLin.none.mem.stat"
  assert_near 0.5 "$(stat_val Finala "$mstat")" 0.20 "noisy fit recovers slope ~0.5 (within 20%)"
  assert_rc 10 "$(stat_val SCount "$mstat")" "fit used all ten points"
  # STDFIT (residual spread) is populated and positive once there IS noise -- it
  # is EMPTY for perfectly-collinear (noiseless) history (see scaling_linear).
  sd=$(stat_val STDFIT "$mstat")
  if [ -n "$sd" ] && awk -v s="$sd" 'BEGIN{exit !(s>0)}'; then _pass "noisy fit reports STDFIT>0 (residual spread detected): $sd"
  else _fail "noisy fit reports STDFIT>0" "STDFIT='$sd'"; fi
  # RSquare is computed correctly (gnuplot-6 'Correlation: r = X' parsed): strong
  # linear trend here -> R^2 well above 0.5. (On unpatched bin/ this reads 0.)
  rsq=$(stat_val RSquare "$mstat")
  if [ -n "$rsq" ] && awk -v r="$rsq" 'BEGIN{exit !(r>0.5)}'; then _pass "RSquare computed (>0.5): $rsq"
  else _fail "RSquare computed (>0.5)" "RSquare='$rsq' (unpatched bin/ records 0 under gnuplot 6)"; fi
  echo
  fi
}

# ===========================================================================
_want scaling_exp && {
  _begin scaling_exp
  if ! command -v gnuplot >/dev/null 2>&1; then
    echo "    (skipped: gnuplot unavailable)"; echo
  else
  ss_new_workspace "$WORK/scaling_exp" >/dev/null
  # generous allocation so every point COMPLETES (no OOM -> no fail-pruning -> fast)
  ss_run_sweep scaling_exp "1000 2000 3000 4000 6000 8000" "sbatch -p short -c 1 --mem 64G -t 900:0"
  jr="$smartSlurmJobRecordDir/jobRecord.txt"
  n=$(awk -F, '$9=="COMPLETED" && $12=="scaleExp"{c++} END{print c+0}' "$jr")
  assert_rc 6 "$n" "six COMPLETED exp records"
  # recorded curve is convex (exp), not linear: the step 6000->8000 jump far
  # exceeds the step 1000->3000 jump (a straight line would make them equal-ish).
  read -r m1 m3 m6 m8 <<<"$(awk -F, '$12=="scaleExp"{a[$2]=$7} END{print a[1000],a[3000],a[6000],a[8000]}' "$jr")"
  if [ $(( m8 - m6 )) -gt $(( (m3 - m1) * 3 )) ]; then
    _pass "recorded mem is convex/exponential (d(6000->8000)=$((m8-m6)) >> d(1000->3000)=$((m3-m1)))"
  else _fail "recorded mem is convex/exponential" "jumps: $((m8-m6)) vs $((m3-m1))"; fi
  # the LINEAR fit under-predicts the exponential truth at a large size
  ss_fit scaleExp none >/dev/null
  big=16000
  truth=$("$HERE/mockSlurm/resourceModel.sh" "$HERE/models/scaling_exp.model" "$big" | cut -d' ' -f1)
  pred=$(estimateResource.sh scaleExp none ,jobSize:$big p 4096 120 2>/dev/null | cut -d' ' -f2)
  echo "    (info) size=$big  linear-prediction=$pred  exponential-truth=$truth"
  if [ "${pred:-0}" -lt "$(( truth / 2 ))" ]; then
    _pass "linear fit under-provisions exp truth at size $big (pred $pred < half of truth $truth)"
  else _fail "linear fit under-provisions exp truth at size $big" "pred=$pred truth=$truth"; fi
  echo
  fi
}

# ===========================================================================
# Contamination of jobRecords: the headline question -- under what conditions do
# bad/error data enter and corrupt the learning DB? Each block constructs ONE
# condition deterministically and asserts the consequence.
# ===========================================================================

_want contaminate_exit0 && {
  _begin contaminate_exit0
  if ! command -v gnuplot >/dev/null 2>&1; then
    echo "    (skipped: gnuplot unavailable)"; echo
  else
  ss_new_workspace "$WORK/contaminate_exit0" >/dev/null
  # clean linear history (mem = 0.5*size + 200)
  for s in 1000 2000 4000 6000 8000; do ss_seed_record cjob none "$s" $((s/2+200)) 30 COMPLETED 1; done
  ss_fit cjob none >/dev/null
  clean=$(stat_val Finala "$smartSlurmJobRecordDir/stats/cjob.none.mem.stat")
  assert_near 0.5 "$clean" 0.05 "clean history fits slope ~0.5"
  # a large-input job that exits 0 having done nothing (or whose inputSize was
  # mis-measured) records a near-zero mem at a big size -> drags the slope down
  ss_seed_record cjob none 10000 60 30 COMPLETED 1
  ss_fit cjob none >/dev/null
  dirty=$(stat_val Finala "$smartSlurmJobRecordDir/stats/cjob.none.mem.stat")
  echo "    (info) slope clean=$clean  contaminated=$dirty"
  if awk -v c="$clean" -v d="$dirty" 'BEGIN{exit !(d < c*0.8)}'; then
    _pass "one exit-0-empty row corrupts the slope downward ($clean -> $dirty)"
  else _fail "exit-0-empty row corrupts the slope" "clean=$clean dirty=$dirty"; fi
  echo
  fi
}

# ===========================================================================
_want contaminate_exclude && {
  _begin contaminate_exclude
  if ! command -v gnuplot >/dev/null 2>&1; then
    echo "    (skipped: gnuplot unavailable)"; echo
  else
  ss_new_workspace "$WORK/contaminate_exclude" >/dev/null
  for s in 1000 2000 4000 6000 8000; do ss_seed_record xjob none "$s" $((s/2+200)) 30 COMPLETED 1; done
  ss_fit xjob none >/dev/null
  before=$(stat_val Finala "$smartSlurmJobRecordDir/stats/xjob.none.mem.stat")
  # add OOM rows with absurd mem at a big size: if they leaked into the fit the
  # slope would explode. They must be ignored (estimateResource greps COMPLETED).
  ss_seed_record xjob none 10000 999999 120 OOM 4
  ss_fit xjob none >/dev/null
  after=$(stat_val Finala "$smartSlurmJobRecordDir/stats/xjob.none.mem.stat")
  if awk -v b="$before" -v a="$after" 'BEGIN{exit !(a==b || (a<b*1.05 && a>b*0.95))}'; then
    _pass "OOM rows excluded from the fit (slope unchanged: $before -> $after)"
  else _fail "OOM rows excluded from the fit" "slope moved $before -> $after"; fi
  echo "    (note) the danger of OOM rows is not direct fit entry but fail-pruning -- see contaminate_prune"
  echo
  fi
}

# ===========================================================================
_want contaminate_prune && {
  _begin contaminate_prune
  if ! command -v gnuplot >/dev/null 2>&1; then
    echo "    (skipped: gnuplot unavailable)"; echo
  else
  ss_new_workspace "$WORK/contaminate_prune" >/dev/null
  jr="$smartSlurmJobRecordDir/jobRecord.txt"
  # progA history: small legit records (survive), one BOGUS under-reporter at a
  # large size (tiny mem -- the exit-0-empty kind), one legit HIGH-mem large
  # record (survives). Plus an UNRELATED progB, and 5 recent OOM fails for progA.
  ss_seed_record progA none 1000 700    15  COMPLETED 3     # legit small
  ss_seed_record progA none 12000 50    30  COMPLETED 1     # BOGUS: huge input, ~0 mem
  ss_seed_record progA none 15000 100000 200 COMPLETED 1    # legit: huge input, huge mem
  ss_seed_record progB none 1000 900    20  COMPLETED 4     # unrelated program
  ss_seed_record progA none 9000 8000   120 OOM 6           # proven-insufficient allocs (>threshold)
  assert_rc 4 "$(awk -F, '$12=="progB"{c++}END{print c+0}' "$jr")" "progB starts with 4 unrelated records"
  # one more real emergent-OOM progA job (size 9000) tips recentFails over threshold
  # True need (80000) sits ABOVE what the smarter estimator will predict for size
  # 9000 from this history (~50k) but BELOW the legit-high record (100000), so the
  # job OOMs on first attempt (triggering the prune) and the OOM alloc stays under
  # the legit-high record (which must be kept).
  cat > "$HERE/models/_prune.model" <<EOF
MODEL_PROGRAM=progA
MODEL_REF=none
MEM_FORM=const
MEM_B=80000
TIME_FORM=const
TIME_B=30
SEED=1
EOF
  pd="$WORK/contaminate_prune/oomjob"; mkdir -p "$pd/.mockstate/jobs"
  ( cd "$pd"
    export MOCK_SLURM_STATE="$pd/.mockstate"
    truncate -s 9000K in.dat
    printf '#!/bin/sh\ninput=$1\n#@1,0,progA,,input,sbatch -p short -c 1 --mem 4G -t 60:0\ncat $input >/dev/null\n\n' > pipeline.sh
    MOCK_FAIL_MODE=oom timeout 115 bash "$SMARTSLURM_ROOT/bin/runAsPipeline" \
      --script "./pipeline.sh in.dat" --sbatch-options "sbatch -p short -c 1 --mem 4G -t 60:0" \
      --tmp noTmp >rap.log 2>&1 )
  rm -f "$HERE/models/_prune.model"
  qa() { awk -F, -v p="$1" -v s="$2" '$12==p && $2==s && $9=="COMPLETED"{c++}END{print c+0}' "$jr"; }
  assert_file "$(ls "$smartSlurmJobRecordDir"/jobRecord.back.*.txt 2>/dev/null | head -1)" "fail-pruning backed up jobRecord.txt"
  # Conservative suspect-removal: keep-filter drops the offending program's failed
  # rows AND any COMPLETED row disproven by a failure (size>=failing size and
  # mem<=proven-insufficient allocation). Unrelated programs and validated
  # (consistent) COMPLETED history are kept.
  assert_rc 4 "$(awk -F, '$12=="progB"{c++}END{print c+0}' "$jr")"  "unrelated progB preserved (4)"
  assert_rc 3 "$(qa progA 1000)"   "progA legit small history preserved (size 1000 x3)"
  assert_rc 1 "$(qa progA 15000)"  "progA legit HIGH-mem large record preserved (size 15000)"
  assert_rc 0 "$(qa progA 12000)"  "progA BOGUS under-reporter removed (size 12000, ~0 mem)"
  echo
  fi
}

# ===========================================================================
_want firstpass_rate && {
  _begin firstpass_rate
  if ! command -v gnuplot >/dev/null 2>&1; then
    echo "    (skipped: gnuplot unavailable)"; echo
  else
  ss_new_workspace "$WORK/firstpass_rate" >/dev/null
  mkdir -p smartSlurmLog
  mdl="$HERE/models/scaling_bignoise.model"
  # Seed ~100 realistic COMPLETED records (model + 8% noise) so the fit's residual
  # percentile is well-resolved. (A 10-point sweep makes p93 collapse to the max.)
  i=1; while [ "$i" -le 100 ]; do
    size=$(( 1000 + (i * 97) % 9000 ))
    read tm tt tr <<<"$("$HERE/mockSlurm/resourceModel.sh" "$mdl" "$size" "$i")"
    ss_seed_record bigJob none "$size" "$tm" "$tt" COMPLETED 1
    i=$((i+1))
  done
  ss_fit bigJob none >/dev/null
  # Simulate "future" jobs (fresh noise draws = unseen seeds). A job passes first
  # try when the estimate at targetPct >= its realized need.
  measure() { # <targetPct> -> integer success percentage over the sample
    local pct="$1" n=0 ok=0 size seed true est
    for size in 1500 3500 5500 7500 9500; do
      for seed in $(seq 201 224); do
        true=$("$HERE/mockSlurm/resourceModel.sh" "$mdl" "$size" "$seed" | cut -d' ' -f1)
        # Measure the B1 MARGIN itself (calculateMemTime), isolated from the separate
        # extraMem safety buffer that estimateResource adds on top.
        est=$(targetFirstPassPct="$pct" calculateMemTime.sh bigJob none "$size" 2>/dev/null | cut -d' ' -f1)
        n=$((n+1)); [ "${est:-0}" -ge "${true:-0}" ] && ok=$((ok+1))
      done
    done
    echo $(( ok * 100 / n ))
  }
  r93=$(measure 93); r70=$(measure 70); r99=$(measure 99)
  echo "    (info) first-pass rate: p70=$r70%  p93=$r93%  p99=$r99%"
  if [ "$r93" -ge 85 ] && [ "$r93" -le 98 ]; then _pass "default (p93) first-pass rate near the ~90-95% target band ($r93%)"
  else _fail "p93 first-pass rate in band" "got $r93%"; fi
  if [ "$r99" -ge "$r93" ] && [ "$r93" -ge "$r70" ]; then _pass "rate increases with target (p70 $r70% <= p93 $r93% <= p99 $r99%)"
  else _fail "rate increases with target" "p70=$r70 p93=$r93 p99=$r99"; fi
  if [ "$r70" -lt 100 ]; then _pass "a lower target leaves more jobs to auto-rerun (p70 $r70% < 100%)"
  else _fail "lower target reduces first-pass coverage" "p70=$r70%"; fi
  echo
  fi
}

# ===========================================================================
_want failure_floor && {
  _begin failure_floor
  ss_new_workspace "$WORK/failure_floor" >/dev/null
  mkdir -p smartSlurmLog
  jr="$smartSlurmJobRecordDir/jobRecord.txt"
  # small legit history, an OOM at size 5000 (alloc mem 8000), an OOT at size 5000
  # (alloc time 120). Both are proven-insufficient at size 5000.
  ss_seed_record ff none 1000 700  15  COMPLETED 3
  ss_seed_record ff none 5000 8000 120 OOM 1
  ss_seed_record ff none 5000 2000 120 OOT 1
  rf=$(awk -F, -v a=ff -v b=none -v t="$(date -d '1 hour ago' +%s)" '$9!="COMPLETED"&&$21>t&&$12==a&&$13==b' "$jr" | wc -l)
  expMem=$(nextMemAfterFail 8000 "$rf"); expTime=$(nextTimeAfterFail 120 "$rf")
  # size >= failing size: must be floored strictly above the proven-insufficient values
  o=$(estimateResource.sh ff none ,jobSize:6000 f 4096 120 2>/dev/null)
  mem=$(echo "$o" | cut -d' ' -f2); tim=$(echo "$o" | cut -d' ' -f3)
  echo "    (info) size 6000 -> mem=$mem time=$tim ; proven-insufficient mem=8000 time=120"
  if [ "${mem:-0}" -gt 8000 ]; then _pass "size 6000: mem ($mem) > proven-insufficient 8000"
  else _fail "size 6000: mem > 8000" "got $mem"; fi
  if [ "${tim:-0}" -gt 120 ]; then _pass "size 6000: time ($tim) > proven-insufficient 120"
  else _fail "size 6000: time > 120" "got $tim"; fi
  assert_rc "$expMem" "$mem" "mem floor == shared auto-rerun escalation of 8000"
  assert_rc "$expTime" "$tim" "time floor == shared auto-rerun escalation of 120"
  # size < failing size: that failure must NOT floor a genuinely smaller input
  o=$(estimateResource.sh ff none ,jobSize:3000 f 4096 120 2>/dev/null); mem=$(echo "$o" | cut -d' ' -f2)
  if [ "${mem:-0}" -lt 8000 ]; then _pass "size 3000 (< failing size): not floored by the OOM ($mem)"
  else _fail "size 3000 not floored" "got $mem"; fi
  echo
}

# ===========================================================================
_want manual_rerun_floor && {
  _begin manual_rerun_floor
  if ! command -v gnuplot >/dev/null 2>&1; then
    echo "    (skipped: gnuplot unavailable)"; echo
  else
  ss_new_workspace "$WORK/manual_rerun_floor" >/dev/null
  jr="$smartSlurmJobRecordDir/jobRecord.txt"
  # A job that really needs ~9000 M; the first run submits 4G (4096) and OOMs.
  cat > "$HERE/models/_mrf.model" <<EOF
MODEL_PROGRAM=mrf
MODEL_REF=none
MEM_FORM=const
MEM_B=9000
TIME_FORM=const
TIME_B=30
SEED=3
EOF
  d="$WORK/manual_rerun_floor/run1"; mkdir -p "$d/.mockstate/jobs"
  ( cd "$d"
    export MOCK_SLURM_STATE="$d/.mockstate"; export SS_MODEL="$HERE/models/_mrf.model"
    truncate -s 5000K in.dat
    printf '#!/bin/sh\ninput=$1\n#@1,0,mrf,,input,sbatch -p short -c 1 --mem 4G -t 60:0\ncat $input >/dev/null\n\n' > pipeline.sh
    timeout 115 bash "$SMARTSLURM_ROOT/bin/runAsPipeline" \
      --script "./pipeline.sh in.dat" --sbatch-options "sbatch -p short -c 1 --mem 4G -t 60:0" \
      --tmp noTmp >rap.log 2>&1 )
  rm -f "$HERE/models/_mrf.model"
  mkdir -p smartSlurmLog
  oomAllocs=$(awk -F, '$12=="mrf" && $9 ~ /OOM/{print $5}' "$jr" | tr '\n' ' ')
  echo "    (info) run1 OOM alloc(s): $oomAllocs"
  if [ -n "$(echo $oomAllocs | tr -d ' ')" ]; then _pass "run1 left at least one OOM record (4096 was insufficient)"
  else _fail "run1 left an OOM record" "none found"; fi
  # The estimate a manual rerun would use (exactly what ssbatch calls) for the same
  # step/size must exceed the proven-insufficient 4096 -- the escalation persists
  # across runs via the failure floor, not just within the auto-requeue.
  est=$(estimateResource.sh mrf none ,jobSize:5000 f 4096 120 2>/dev/null | cut -d' ' -f2)
  echo "    (info) manual-rerun estimate for size 5000 = $est M (initial 4096 OOM'd)"
  if [ "${est:-0}" -gt 4096 ]; then _pass "manual rerun estimate ($est) > proven-insufficient 4096"
  else _fail "manual rerun estimate > 4096" "got $est"; fi
  echo
  fi
}

# ===========================================================================
_want constant_model && {
  _begin constant_model
  if ! command -v gnuplot >/dev/null 2>&1; then
    echo "    (skipped: gnuplot unavailable)"; echo
  else
  ss_new_workspace "$WORK/constant_model" >/dev/null
  ss_run_sweep scaling_const_mem "1000 2000 3000 4000 6000 8000 10000" "sbatch -p short -c 1 --mem 16G -t 1200:0"
  ss_fit constMem none >/dev/null
  # mem prediction must be ~constant across sizes (constant model), NOT negative-
  # sloped and NOT floored to 100. time prediction must scale up with size.
  ms=$(estimateResource.sh constMem none ,jobSize:1000 p 4096 120 2>/dev/null)
  ml=$(estimateResource.sh constMem none ,jobSize:9000 p 4096 120 2>/dev/null)
  mem_s=$(echo "$ms" | cut -d' ' -f2); mem_l=$(echo "$ml" | cut -d' ' -f2)
  time_s=$(echo "$ms" | cut -d' ' -f3); time_l=$(echo "$ml" | cut -d' ' -f3)
  echo "    (info) mem: size1000=$mem_s size9000=$mem_l | time: size1000=$time_s size9000=$time_l"
  if [ "${mem_s:-0}" -gt 1000 ]; then _pass "constant mem not collapsed to floor (size 1000 -> $mem_s M)"
  else _fail "constant mem not collapsed to floor" "size 1000 -> $mem_s"; fi
  assert_near "$mem_s" "$mem_l" 0.15 "mem is ~constant across sizes ($mem_s vs $mem_l)"
  if [ "${time_l:-0}" -gt "${time_s:-0}" ]; then _pass "time scales up with size ($time_s -> $time_l)"
  else _fail "time scales up with size" "$time_s -> $time_l"; fi
  echo
  fi
}

# ===========================================================================
_want negslope_robust && {
  _begin negslope_robust
  ss_new_workspace "$WORK/negslope_robust" >/dev/null
  sdir="$smartSlurmJobRecordDir/stats"; mkdir -p "$sdir"
  mkstat() { # prog ref slope intercept max stdfit scount
    printf 'Finala=%s\nFinalb=%s\nMaximum=%s\nSTDFIT=%s\nSCount=%s\n' "$3" "$4" "$5" "$6" "$7" > "$sdir/$1.$2.mem.stat"
    printf 'Finala=%s\nFinalb=%s\nMaximum=%s\nSTDFIT=%s\nSCount=%s\nRSquare=%s\n' "$3" "$4" "$5" "$6" "$7" 0.5 > "$sdir/$1.$2.time.stat"
  }
  isint() { printf '%s' "$1" | grep -qE '^[0-9]+$'; }
  # LINE path (positive slope) with a result in (-1,0): old code produced mem="-"
  # (bc "-.3" -> ${mem%.*} -> "-") and passed garbage downstream. Must floor to 100.
  mkstat P none 0.5 -0.7 100000 0 10
  m=$(calculateMemTime.sh P none 1 2>/dev/null | awk '{print $1}')
  if isint "$m" && [ "$m" -eq 100 ]; then _pass "line path, sub-1 negative result floors to 100 (was bc '-' bug): $m"
  else _fail "line path sub-1 negative floors to 100" "got '$m'"; fi
  # LINE path, large negative result -> floor 100 (not a negative number).
  mkstat P none 0.5 -100000 100000 0 10
  m=$(calculateMemTime.sh P none 1 2>/dev/null | awk '{print $1}')
  if isint "$m" && [ "$m" -eq 100 ]; then _pass "line path, large negative result floors to 100: $m"
  else _fail "line path large negative floors to 100" "got '$m'"; fi
  # Empty STDFIT (perfectly collinear/noiseless history) must not break arithmetic.
  mkstat P none 0.5 200 100000 '' 10
  m=$(calculateMemTime.sh P none 4000 2>/dev/null | awk '{print $1}')
  assert_rc 2200 "$m" "empty STDFIT handled (line mem=2200, no error)"
  # CONSTANT path (slope <= 0) with raw data present -> valid integer percentile.
  mkstat P none -0.3 50 100000 0 10
  printf '%s\n' "1000 1500" "2000 1600" "4000 1550" "8000 1700" > "$sdir/P.none.mem"
  m=$(calculateMemTime.sh P none 5000 2>/dev/null | awk '{print $1}')
  if isint "$m" && [ "$m" -ge 100 ]; then _pass "constant path (neg slope) returns valid integer from data ($m)"
  else _fail "constant path returns valid integer" "got '$m'"; fi
  echo
}

# ===========================================================================
# Activation / setup feature (smartSlurmSetup + the `smartslurm` shell hook).
# These blocks need no Slurm and no gnuplot: they exercise the rc-block writer
# and the activate/deactivate/uninstall shell functions in a pristine sandbox
# shell (ss_scn), asserting on PATH, PS1, and the rc file.
# ===========================================================================
_want setup_rc_block && {
  _begin setup_rc_block
  h="$WORK/setup_rc_block_home"; rm -rf "$h"; mkdir -p "$h"; o="$h/out.txt"
  # Seed the rc with the user's own content, to prove setup preserves it.
  printf 'export USER_SENTINEL=1\nalias ll="ls -l"\n' > "$h/.bashrc"
  ss_scn "$h" > "$o" 2>&1 <<'SCN'
bash "$SMARTSLURM_HOME/smartSlurmSetup.sh" >/dev/null 2>&1
echo "BLOCKS1=$(grep -c '>>> smartSlurm initialize >>>' "$HOME/.bashrc")"
bash "$SMARTSLURM_HOME/smartSlurmSetup.sh" >/dev/null 2>&1   # re-run: must not duplicate
echo "BLOCKS2=$(grep -c '>>> smartSlurm initialize >>>' "$HOME/.bashrc")"
echo "USERLINE=$(grep -c 'USER_SENTINEL' "$HOME/.bashrc")"
echo "ALIASLINE=$(grep -c 'alias ll' "$HOME/.bashrc")"
echo "HOMEREF=$(grep -c "SMARTSLURM_HOME=\"$SMARTSLURM_HOME\"" "$HOME/.bashrc")"
ls "$HOME"/.bashrc.smartslurm.bak.* >/dev/null 2>&1 && echo "BACKUP=yes" || echo "BACKUP=no"
SCN
  assert_contains "$o" "BLOCKS1=1"  "setup writes exactly one initialize block"
  assert_contains "$o" "BLOCKS2=1"  "re-running setup is idempotent (no duplicate block)"
  assert_contains "$o" "USERLINE=1" "user's existing export line is preserved"
  assert_contains "$o" "ALIASLINE=1" "user's existing alias line is preserved"
  assert_contains "$o" "HOMEREF=1"  "block points SMARTSLURM_HOME at this checkout"
  assert_contains "$o" "BACKUP=yes" "re-run backs up the rc before editing"
  echo
}

# ===========================================================================
_want activate_path && {
  _begin activate_path
  h="$WORK/activate_path_home"; rm -rf "$h"; mkdir -p "$h"; o="$h/out.txt"
  ss_scn "$h" > "$o" 2>&1 <<'SCN'
. "$SMARTSLURM_HOME/bin/smartslurm.sh"
smartslurm activate
echo "SBATCH=$(command -v sbatch)"
echo "SSBATCH=$(command -v ssbatch)"
echo "PATHVAL=$PATH"
SCN
  assert_contains "$o" "now runs smart ssbatch"              "activate announces the sbatch aliasing"
  assert_contains "$o" "SBATCH=$SMARTSLURM_ROOT/sbatchBin/sbatch" "plain sbatch resolves into sbatchBin (aliased)"
  assert_contains "$o" "SSBATCH=$SMARTSLURM_ROOT/bin/ssbatch"     "ssbatch resolves into bin/"
  assert_contains "$o" "PATHVAL=$SMARTSLURM_ROOT/bin:$SMARTSLURM_ROOT/sbatchBin:" "both dirs are prepended, bin first"
  echo
}

# ===========================================================================
_want activate_prompt && {
  _begin activate_prompt
  h="$WORK/activate_prompt_home"; rm -rf "$h"; mkdir -p "$h"; o="$h/out.txt"
  ss_scn "$h" > "$o" 2>&1 <<'SCN'
. "$SMARTSLURM_HOME/bin/smartslurm.sh"
PS1='base$ '
smartslurm activate >/dev/null
echo "PS1A=<<<$PS1>>>"
smartslurm activate >/dev/null 2>&1        # second activate must not re-prefix
echo "PS1B=<<<$PS1>>>"
echo "COUNT=$(printf '%s' "$PS1" | grep -oF '{smartSlurm} ' | wc -l | tr -d ' ')"
SCN
  assert_contains "$o" "PS1A=<<<{smartSlurm} base"  "activate prefixes prompt with {smartSlurm}"
  assert_lacks    "$o" "(env)"                      "prompt token uses braces, not conda-style parens"
  assert_contains "$o" "COUNT=1"                    "double activate leaves exactly one prefix (idempotent)"
  echo
}

# ===========================================================================
_want deactivate_surgical && {
  _begin deactivate_surgical
  h="$WORK/deactivate_surgical_home"; rm -rf "$h"; mkdir -p "$h"; o="$h/out.txt"
  ss_scn "$h" > "$o" 2>&1 <<'SCN'
. "$SMARTSLURM_HOME/bin/smartslurm.sh"
PS1='base$ '
smartslurm activate >/dev/null
# Things added AFTER activation (a module load, a conda activate, a user export):
export PATH="/opt/added/bin:$PATH"
PS1="(env) $PS1"
smartslurm deactivate >/dev/null
echo "PATHVAL=$PATH"
echo "SBATCH=$(command -v sbatch)"
echo "PS1VAL=<<<$PS1>>>"
SCN
  assert_contains "$o" "/opt/added/bin"                     "post-activation PATH addition is preserved"
  assert_lacks    "$o" "$SMARTSLURM_ROOT/bin"               "smartSlurm bin/ removed from PATH"
  assert_lacks    "$o" "$SMARTSLURM_ROOT/sbatchBin"         "smartSlurm sbatchBin/ removed from PATH"
  assert_lacks    "$o" "/sbatchBin/sbatch"                  "sbatch alias is gone (system sbatch restored)"
  assert_contains "$o" "PS1VAL=<<<(env) base"               "post-activation prompt prefix '(env) ' is preserved"
  assert_lacks    "$o" "{smartSlurm}"                       "our prompt token removed, leaving others intact"
  echo
}

# ===========================================================================
_want conda_coexist && {
  _begin conda_coexist
  h="$WORK/conda_coexist_home"; rm -rf "$h"; mkdir -p "$h"; o="$h/out.txt"
  ss_scn "$h" > "$o" 2>&1 <<'SCN'
. "$SMARTSLURM_HOME/bin/smartslurm.sh"
PS1='base$ '
# conda activates FIRST (prepends its dir + "(env) " prompt), then smartSlurm.
export PATH="/opt/conda/envs/x/bin:$PATH"
PS1="(env) $PS1"
smartslurm activate >/dev/null
echo "BOTH_PATH=$PATH"
echo "BOTH_PS1=<<<$PS1>>>"
smartslurm deactivate >/dev/null            # ss off; conda must be untouched
echo "AFTER_PATH=$PATH"
echo "AFTER_PS1=<<<$PS1>>>"
SCN
  assert_contains "$o" "BOTH_PATH=$SMARTSLURM_ROOT/bin"   "with both active, smartSlurm is on PATH"
  assert_contains "$o" "/opt/conda/envs/x/bin"            "with both active, conda is still on PATH"
  assert_contains "$o" "BOTH_PS1=<<<{smartSlurm} (env)"   "both prompt prefixes present, neither clobbered"
  assert_contains "$o" "AFTER_PATH=/opt/conda/envs/x/bin" "after ss deactivate, conda PATH entry remains"
  assert_lacks    "$o" "AFTER_PATH=$SMARTSLURM_ROOT/bin"  "after ss deactivate, smartSlurm is off PATH"
  assert_contains "$o" "AFTER_PS1=<<<(env) base"          "after ss deactivate, conda prompt prefix remains"
  echo
}

# ===========================================================================
_want uninstall_precise && {
  _begin uninstall_precise
  h="$WORK/uninstall_precise_home"; rm -rf "$h"; mkdir -p "$h"; o="$h/out.txt"
  printf 'line-before\nexport USER_SENTINEL=1\nline-after\n' > "$h/.bashrc"
  ss_scn "$h" > "$o" 2>&1 <<'SCN'
bash "$SMARTSLURM_HOME/smartSlurmSetup.sh" >/dev/null 2>&1
. "$SMARTSLURM_HOME/bin/smartslurm.sh"
smartslurm uninstall >/dev/null 2>&1
echo "BLOCKLEFT=$(grep -c 'smartSlurm initialize' "$HOME/.bashrc")"
echo "BEFORE=$(grep -c '^line-before$' "$HOME/.bashrc")"
echo "MID=$(grep -c 'USER_SENTINEL' "$HOME/.bashrc")"
echo "AFTER=$(grep -c '^line-after$' "$HOME/.bashrc")"
ls "$HOME"/.bashrc.smartslurm.bak.* >/dev/null 2>&1 && echo "BACKUP=yes" || echo "BACKUP=no"
smartslurm uninstall 2>&1 | grep -qi 'no initialize block' && echo "ABSENT_OK=yes" || echo "ABSENT_OK=no"
SCN
  assert_contains "$o" "BLOCKLEFT=0"  "uninstall removes the initialize block entirely"
  assert_contains "$o" "BEFORE=1"     "line before the block is preserved"
  assert_contains "$o" "MID=1"        "user line is preserved"
  assert_contains "$o" "AFTER=1"      "line after the block is preserved"
  assert_contains "$o" "BACKUP=yes"   "uninstall backs up the rc first"
  assert_contains "$o" "ABSENT_OK=yes" "uninstall on an already-clean rc reports nothing-to-remove"
  echo
}

# ===========================================================================
_want unexport_deprecated && {
  _begin unexport_deprecated
  h="$WORK/unexport_deprecated_home"; rm -rf "$h"; mkdir -p "$h"; o="$h/out.txt"
  ss_scn "$h" > "$o" 2>&1 <<'SCN'
. "$SMARTSLURM_HOME/bin/smartslurm.sh"
. "$SMARTSLURM_HOME/bin/unExport"
PS1='base$ '
smartslurm activate >/dev/null
unExport                                    # deprecated path: warn + delegate
echo "SBATCH=$(command -v sbatch)"
SCN
  assert_contains "$o" "deprecated"               "unExport prints a deprecation warning"
  assert_contains "$o" "smartSlurm deactivated"   "unExport delegates to smartslurm deactivate"
  assert_lacks    "$o" "/sbatchBin/sbatch"        "after unExport, the smart sbatch alias is gone"
  echo
}

# ===========================================================================
# Estimator first-pass / extrapolation / auto-resubmit measurement (Phase 1).
# These quantify how the patched estimator behaves at the 90-95% first-pass and
# 100%-on-resubmit targets, in-range and when extrapolating far beyond training.
# Several assertions pin the CURRENT (imperfect) behavior as a baseline; the
# purpose text marks those "[Phase 1 baseline]".
# ===========================================================================
_want firstpass_fullpath && {
  _begin firstpass_fullpath
  if ! command -v gnuplot >/dev/null 2>&1; then
    echo "    (skipped: gnuplot unavailable)"; echo
  else
  ss_new_workspace "$WORK/firstpass_fullpath" >/dev/null; mkdir -p smartSlurmLog
  mdl="$HERE/models/scaling_bignoise.model"
  i=1; while [ "$i" -le 100 ]; do
    size=$(( 1000 + (i*97) % 9000 ))
    read tm tt tr <<<"$("$HERE/mockSlurm/resourceModel.sh" "$mdl" "$size" "$i")"
    ss_seed_record bigJob none "$size" "$tm" "$tt" COMPLETED 1
    i=$((i+1))
  done
  ss_fit bigJob none >/dev/null
  szs="1500 3500 5500 7500 9500"; sds="$(seq 201 216)"
  fm=$(ss_firstpass_rate bigJob none "$mdl" 93 mem  full "$szs" "$sds")
  ft=$(ss_firstpass_rate bigJob none "$mdl" 93 time full "$szs" "$sds")
  cm=$(ss_firstpass_rate bigJob none "$mdl" 93 mem  cmt  "$szs" "$sds")
  echo "    (info) full-path first-pass: mem=$fm%  time=$ft%   (margin-only mem=$cm%)"
  if [ "$fm" -ge 85 ]; then _pass "full-path (estimateResource) mem first-pass in target band ($fm%)"
  else _fail "full-path mem first-pass in band" "got $fm%"; fi
  if [ "$ft" -ge 85 ]; then _pass "full-path TIME first-pass in target band ($ft%)"
  else _fail "full-path time first-pass in band" "got $ft%"; fi
  if [ "$fm" -ge "$cm" ]; then _pass "full path (+extraMem/floor) covers at least as well as the margin alone ($fm% >= $cm%)"
  else _fail "full path >= margin-only" "full=$fm cmt=$cm"; fi
  echo
  fi
}

# ===========================================================================
_want extrapolate_linear && {
  _begin extrapolate_linear
  if ! command -v gnuplot >/dev/null 2>&1; then
    echo "    (skipped: gnuplot unavailable)"; echo
  else
  ss_new_workspace "$WORK/extrapolate_linear" >/dev/null; mkdir -p smartSlurmLog
  mdl="$HERE/models/extrap_lin.model"
  # Train on SMALL sizes 100..1000 (dozens of points => SCount>9), linear+noise.
  i=1; while [ "$i" -le 40 ]; do
    size=$(( 100 + (i*23) % 900 ))
    read tm tt tr <<<"$("$HERE/mockSlurm/resourceModel.sh" "$mdl" "$size" "$i")"
    ss_seed_record extraLin none "$size" "$tm" "$tt" COMPLETED 1
    i=$((i+1))
  done
  ss_fit extraLin none >/dev/null
  o5=$(calculateMemTime.sh extraLin none 5000 2>/dev/null | cut -d' ' -f1)
  case "$o5" in ''|*[!0-9]*) r5=outOfRange;; *) r5=$o5;; esac
  if [ "$r5" != outOfRange ]; then _pass "SCount>9 extrapolates beyond training max (size 5000 -> $r5 M, not outOfRange)"
  else _fail "extrapolates beyond max with SCount>9" "got outOfRange"; fi
  e2=$(calculateMemTime.sh extraLin none 2000 2>/dev/null | cut -d' ' -f1)
  e50=$(calculateMemTime.sh extraLin none 50000 2>/dev/null | cut -d' ' -f1)
  case "$e2" in ''|*[!0-9]*) e2=0;; esac; case "$e50" in ''|*[!0-9]*) e50=0;; esac
  if [ "$e50" -gt "$e2" ]; then _pass "estimate grows with size (2k=$e2 -> 50k=$e50 M)"
  else _fail "estimate monotonic in size" "2k=$e2 50k=$e50"; fi
  for sc in 2000 5000 10000 25000 50000; do
    r=$(ss_firstpass_rate extraLin none "$mdl" 93 mem cmt "$sc" "$(seq 201 212)")
    echo "    (info) linear-truth first-pass at size $sc: $r%"
  done
  t10=$("$HERE/mockSlurm/resourceModel.sh" "$mdl" 10000 207 | cut -d' ' -f1)
  p10=$(calculateMemTime.sh extraLin none 10000 2>/dev/null | cut -d' ' -f1)
  case "$p10" in ''|*[!0-9]*) p10=0;; esac
  if [ "$p10" -ge $(( t10 / 2 )) ]; then _pass "linear 10x estimate not grossly short (pred $p10 >= half of truth $t10)"
  else _fail "linear 10x coverage" "pred=$p10 truth=$t10"; fi
  echo
  fi
}

# ===========================================================================
_want extrapolate_superlinear && {
  _begin extrapolate_superlinear
  if ! command -v gnuplot >/dev/null 2>&1; then
    echo "    (skipped: gnuplot unavailable)"; echo
  else
  ss_new_workspace "$WORK/extrapolate_superlinear" >/dev/null; mkdir -p smartSlurmLog
  mdl="$HERE/models/extrap_quad.model"
  i=1; while [ "$i" -le 40 ]; do
    size=$(( 100 + (i*23) % 900 ))
    read tm tt tr <<<"$("$HERE/mockSlurm/resourceModel.sh" "$mdl" "$size" "$i")"
    ss_seed_record extraQuad none "$size" "$tm" "$tt" COMPLETED 1
    i=$((i+1))
  done
  ss_fit extraQuad none >/dev/null
  for sc in 10000 50000; do
    tr=$("$HERE/mockSlurm/resourceModel.sh" "$mdl" "$sc" 207 | cut -d' ' -f1)
    pr=$(calculateMemTime.sh extraQuad none "$sc" 2>/dev/null | cut -d' ' -f1)
    case "$pr" in ''|*[!0-9]*) pr=0;; esac
    echo "    (info) quad truth@${sc}=${tr}M  linear-fit pred=${pr}M  shortfall=$(( tr / (pr+1) ))x"
  done
  tr50=$("$HERE/mockSlurm/resourceModel.sh" "$mdl" 50000 207 | cut -d' ' -f1)
  pr50=$(calculateMemTime.sh extraQuad none 50000 2>/dev/null | cut -d' ' -f1)
  case "$pr50" in ''|*[!0-9]*) pr50=0;; esac
  if [ "$pr50" -lt "$tr50" ]; then _pass "[baseline] straight-line fit UNDER-provisions quadratic truth at 50x (pred $pr50 < truth $tr50) -- the gap Phase 2 must close"
  else _fail "baseline: under-provision at 50x expected" "pred=$pr50 truth=$tr50"; fi
  r=$(ss_firstpass_rate extraQuad none "$mdl" 93 mem cmt "10000 50000" "$(seq 201 212)")
  echo "    (info) quad-truth first-pass at 10x+50x: $r%"
  echo
  fi
}

# ===========================================================================
_want extrapolate_lowdata && {
  _begin extrapolate_lowdata
  if ! command -v gnuplot >/dev/null 2>&1; then
    echo "    (skipped: gnuplot unavailable)"; echo
  else
  ss_new_workspace "$WORK/extrapolate_lowdata" >/dev/null; mkdir -p smartSlurmLog
  mdl="$HERE/models/extrap_lin.model"
  for size in 200 400 600 800 1000; do           # only 5 points => SCount<=9
    read tm tt tr <<<"$("$HERE/mockSlurm/resourceModel.sh" "$mdl" "$size" 5)"
    ss_seed_record extraLow none "$size" "$tm" "$tt" COMPLETED 1
  done
  ss_fit extraLow none >/dev/null
  # Fixed behavior: thin data beyond the range no longer returns outOfRange/bare
  # default -- it returns a trend-scaled estimate (in-range value x size/Maximum).
  m1=$(calculateMemTime.sh extraLow none 10000 2>/dev/null | awk '{print $1}')
  case "$m1" in
    ''|*[!0-9]*) _fail "thin data returns a numeric estimate (not outOfRange)" "got '$m1'";;
    *) _pass "thin data beyond range returns a numeric trend-scaled estimate, not outOfRange (${m1}M)";;
  esac
  inrange=$(calculateMemTime.sh extraLow none 1000 2>/dev/null | cut -d' ' -f1)
  m10=$(calculateMemTime.sh extraLow none 10000 2>/dev/null | cut -d' ' -f1)
  m50=$(calculateMemTime.sh extraLow none 50000 2>/dev/null | cut -d' ' -f1)
  case "$m10" in ''|*[!0-9]*) m10=0;; esac; case "$m50" in ''|*[!0-9]*) m50=0;; esac
  full=$(estimateResource.sh extraLow none ,jobSize:10000 f 4096 120 2>/dev/null | cut -d' ' -f2)
  case "$full" in ''|*[!0-9]*) full=0;; esac
  echo "    (info) thin-data trend-scaled: in-range(1x@1000)=${inrange}M  10x=${m10}M  50x=${m50}M  full-path(10x)=${full}M"
  if [ "$m10" -gt 4096 ]; then _pass "thin-data 10x estimate beats the bare default (${m10}M > 4096M)"
  else _fail "thin-data 10x beats default" "got ${m10}M"; fi
  if [ "$m50" -gt "$m10" ]; then _pass "estimate grows with size on thin data (10x=${m10} -> 50x=${m50}M)"
  else _fail "thin-data estimate monotonic in size" "10x=${m10} 50x=${m50}"; fi
  # trend-scaled ~= in-range * ratio (anchor at max, linear scale); ratio 10x here
  lo=$(( inrange * 10 * 7 / 10 )); hi=$(( inrange * 10 * 13 / 10 ))
  if [ "$m10" -ge "$lo" ] && [ "$m10" -le "$hi" ]; then _pass "10x estimate tracks in-range x (size/Maximum) within +-30% (${m10}M vs ~$(( inrange * 10 ))M)"
  else _fail "trend-scaled tracks in-range*ratio" "m10=${m10} expected ~$(( inrange*10 )) [$lo..$hi]"; fi
  echo
  fi
}

# ===========================================================================
_want resubmit_convergence && {
  _begin resubmit_convergence
  ss_new_workspace "$WORK/resubmit_convergence" >/dev/null
  jr="$smartSlurmJobRecordDir/jobRecord.txt"
  # A job that really needs ~7000 M; the first submit uses 4096 (4G) and OOMs,
  # then auto-requeue escalates until it fits -- within the 3-requeue cap.
  cat > "$HERE/models/_conv.model" <<EOF
MODEL_PROGRAM=conv
MODEL_REF=none
MEM_FORM=const
MEM_B=7000
TIME_FORM=const
TIME_B=30
SEED=5
EOF
  d="$WORK/resubmit_convergence/run"; mkdir -p "$d/.mockstate/jobs"
  ( cd "$d"
    export MOCK_SLURM_STATE="$d/.mockstate"; export SS_MODEL="$HERE/models/_conv.model"
    truncate -s 5000K in.dat
    printf '#!/bin/sh\ninput=$1\n#@1,0,conv,,input,sbatch -p short -c 1 --mem 4G -t 60:0\ncat $input >/dev/null\n\n' > pipeline.sh
    timeout 115 bash "$SMARTSLURM_ROOT/bin/runAsPipeline" \
      --script "./pipeline.sh in.dat" --sbatch-options "sbatch -p short -c 1 --mem 4G -t 60:0" \
      --tmp noTmp >rap.log 2>&1 )
  rm -f "$HERE/models/_conv.model"
  ladder=$(ss_alloc_ladder conv)
  echo "    (info) auto-requeue alloc ladder (mem:status): $ladder"
  succ=$(find "$d" -name '*conv*.success' 2>/dev/null | head -1)
  if [ -n "$succ" ]; then _pass "under-provisioned job converged to success via auto-requeue"
  else _fail "auto-requeue converges to success" "no .success; ladder=$ladder"; fi
  firstoom=$(awk -F, '$12=="conv" && $9 ~ /OOM/{print $5; exit}' "$jr")
  if [ -n "$firstoom" ]; then _pass "first attempt recorded OOM (initial 4096 was insufficient: ${firstoom}M)"
  else _fail "first attempt OOMed" "no OOM record"; fi
  mono=$(echo "$ladder" | awk '{p=0; for(i=1;i<=NF;i++){split($i,a,":"); if(a[1]+0<p){print "no"; exit} p=a[1]+0} print "yes"}')
  if [ "$mono" = yes ]; then _pass "escalation ladder never regresses (monotonically non-decreasing allocs)"
  else _fail "ladder monotonic" "$ladder"; fi
  echo
}

# ===========================================================================
_want manual_resubmit_sequence && {
  _begin manual_resubmit_sequence
  ss_new_workspace "$WORK/manual_resubmit_sequence" >/dev/null
  jr="$smartSlurmJobRecordDir/jobRecord.txt"
  # A need far above what a single run's 3-requeue auto-chain reaches from 4096,
  # so completion requires several sequential MANUAL resubmissions -- each one
  # re-estimating from the accumulated OOM history (failure floor), so it starts
  # above the previous proven-insufficient allocation.
  cat > "$HERE/models/_mseq.model" <<EOF
MODEL_PROGRAM=mseq
MODEL_REF=none
MEM_FORM=const
MEM_B=150000
TIME_FORM=const
TIME_B=30
SEED=5
EOF
  base="$WORK/manual_resubmit_sequence"; rounds=0; converged=0
  for r in 1 2 3 4 5 6 7 8; do
    rounds=$r
    d="$base/m$r"; mkdir -p "$d/.mockstate/jobs"
    ( cd "$d"
      export MOCK_SLURM_STATE="$d/.mockstate"; export SS_MODEL="$HERE/models/_mseq.model"
      truncate -s 5000K in.dat
      printf '#!/bin/sh\ninput=$1\n#@1,0,mseq,,input,sbatch -p short -c 1 --mem 4G -t 60:0\ncat $input >/dev/null\n\n' > pipeline.sh
      timeout 115 bash "$SMARTSLURM_ROOT/bin/runAsPipeline" \
        --script "./pipeline.sh in.dat" --sbatch-options "sbatch -p short -c 1 --mem 4G -t 60:0" \
        --tmp noTmp >rap.log 2>&1 )
    find "$d" -name '*mseq*.success' 2>/dev/null | grep -q . && { converged=1; break; }
  done
  rm -f "$HERE/models/_mseq.model"
  ladder=$(ss_alloc_ladder mseq)
  echo "    (info) manual rounds used: $rounds; full alloc ladder (mem:status): $ladder"
  if [ "$converged" = 1 ]; then _pass "bootstrapped to success across $rounds sequential manual submission(s)"
  else _fail "converges within the manual budget" "ladder=$ladder"; fi
  if [ "$rounds" -ge 2 ]; then _pass "auto-requeue alone did not suffice; multiple manual rounds were needed ($rounds)"
  else _fail "scenario needs >1 manual round" "converged in $rounds (raise MEM_B)"; fi
  mono=$(echo "$ladder" | awk '{p=0; for(i=1;i<=NF;i++){split($i,a,":"); if(a[1]+0<p){print "no"; exit} p=a[1]+0} print "yes"}')
  if [ "$mono" = yes ]; then _pass "resources never regress across the whole manual+auto sequence (monotonic)"
  else _fail "monotonic non-regression across manual rounds" "$ladder"; fi
  echo
}

# ===========================================================================
# Efficiency reporting (feat/efficiency-reporting).
# ===========================================================================
_want efficiency_file && {
  _begin efficiency_file
  ss_new_workspace "$WORK/efficiency_file"
  cp "$HERE"/fixtures/clean_2step/* .
  chmod +x ./*.sh 2>/dev/null || true
  ss_run_pipeline "./pipeline.sh 777" "sbatch -p short -c 1 --mem 4G -t 50:0"; rc=$?
  eff="smartSlurmLog/efficiency.tsv"
  assert_file     "$eff"                "per-run efficiency.tsv written"
  assert_contains "$eff" "allocMemMB"   "header names carry units (allocMemMB)"
  assert_contains "$eff" "usedTimeMin"  "header names carry units (usedTimeMin)"
  # trailing human-readable ISO date column (alongside the machine epoch)
  if awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)if($i=="date")d=i} NR>1 && d && $d ~ /^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}$/{c++} END{exit !(c>=1)}' "$eff"; then
    _pass "rows carry a human-readable ISO date column"
  else _fail "ISO date column" "$(head -2 "$eff")"; fi
  # one row per job (2-step fixture = 2 producers + 1 merge = 3 jobs)
  rows=$(grep -cvE '^epoch' "$eff" 2>/dev/null)
  if [ "${rows:-0}" -ge 3 ]; then _pass "one row per job (unconditional, not jobRecord-gated): $rows rows"
  else _fail "row per job" "got ${rows:-0}"; fi
  # resource columns are numeric on COMPLETED rows (alloc/used mem, alloc/used time)
  ok=$(awk -F'\t' '$6=="COMPLETED" && $7 ~ /^[0-9]+$/ && $8 ~ /^[0-9]+$/ && $9 ~ /^[0-9]+$/ && $10 ~ /^[0-9]+$/{c++}END{print c+0}' "$eff")
  if [ "${ok:-0}" -ge 1 ]; then _pass "COMPLETED rows carry numeric alloc/used mem+time ($ok rows)"
  else _fail "numeric resource columns on completed rows" "none"; fi
  # sample key populated for every data row
  blank=$(awk -F'\t' 'NR>1 && $3==""{c++}END{print c+0}' "$eff")
  if [ "${blank:-1}" -eq 0 ]; then _pass "sample key populated on every row"
  else _fail "sample key populated" "$blank blank"; fi
  # used <= alloc sanity for completed rows (efficiency <= 100%)
  bad=$(awk -F'\t' '$6=="COMPLETED" && ($8+0>$7+0 || $10+0>$9+0){c++}END{print c+0}' "$eff")
  if [ "${bad:-1}" -eq 0 ]; then _pass "completed rows have used <= alloc (mem and time)"
  else _fail "used <= alloc on completed rows" "$bad violating rows"; fi
  echo
}

# ===========================================================================
_want checkrun_efficiency && {
  _begin checkrun_efficiency
  # Pin a wide terminal so the label width is deterministic (it now adapts to the
  # terminal); tests that need a narrow width set CHECKRUN_SAMPLE_W inline.
  export COLUMNS=200
  ss_new_workspace "$WORK/checkrun_efficiency" >/dev/null
  mkdir -p smartSlurmLog
  eff=smartSlurmLog/efficiency.tsv
  printf 'epoch\tflag\tsample\tprogram\tref\tstatus\tallocMem\tusedMem\tallocTime\tusedTime\tinputs\n' > "$eff"
  # sampleA: 2 steps. time SUM=(420+60)/(600+120)=67%; mem PEAK=step2 22077/22528=98%
  # (the peak must ignore step1's smaller 10240 alloc). sampleB: single step at 50%.
  printf '1\t1.0.trim.sampleA\tsampleA\ttrim\tnone\tCOMPLETED\t10240\t9216\t600\t420\t,a.fq\n'   >> "$eff"
  printf '2\t2.1.align.sampleA\tsampleA\talign\tnone\tCOMPLETED\t22528\t22077\t120\t60\t,a.bam\n' >> "$eff"
  printf '1\t1.0.trim.sampleB\tsampleB\ttrim\tnone\tCOMPLETED\t10240\t5120\t600\t300\t,b.fq\n'   >> "$eff"
  printf '2\t2.1.align.sampleB\tsampleB\talign\tnone\tCOMPLETED\t22528\t11264\t120\t90\t,b.bam\n' >> "$eff"

  # ---- OVERVIEW: Overall total + per-sample rows + Samples/Steps lists --------
  checkRunGrid --efficiency --no-color --log smartSlurmLog > ov.txt 2>&1
  assert_contains ov.txt "Overall" "overview shows an Overall total row"
  if grep -qE '\(time\) .* 67%' ov.txt; then _pass "sampleA time% is the SUM across steps (67%)"
  else _fail "time% = sum of steps" "$(grep -E '\(time\)' ov.txt | head -3)"; fi
  if grep -qE '\(mem\) .* 98%' ov.txt; then _pass "sampleA mem% is the PEAK step, not the low step (98%)"
  else _fail "mem% = peak step" "$(grep -E '\(mem\)' ov.txt | head -3)"; fi
  # (the Overall mem row carries a blank label by design -- identify it by its
  # unique 74% + 44.0G aggregate alloc)
  if grep -qE '\(mem\) .* 74% .* 44.0G' ov.txt; then _pass "Overall mem% = sum of per-sample peaks (74%)"
  else _fail "Overall mem aggregate" "$(grep -E '\(mem\)' ov.txt | head -2)"; fi
  # eff/used/alloc columns come AFTER the bar, with real numbers
  assert_contains ov.txt "12.0h" "used/alloc column shows time in hours (12.0h)"
  assert_contains ov.txt "22.0G" "used/alloc column shows mem in GB (22.0G)"
  # sample names are shortened (common prefix stripped -> "sampleA" shows as "A"),
  # each tagged with a base-26 UPPERCASE letter, and the footer gives the count.
  if grep -q 'Samples:  2 total' ov.txt && grep -qE '\(A\).*\(time\)' ov.txt && grep -qE '\(B\).*\(time\)' ov.txt; then _pass "overview lists samples with base-26 letters (A,B) + count footer"
  else _fail "samples lettered list" "$(grep -E 'Samples:|\(time\)' ov.txt)"; fi
  if grep -q 'Steps:' ov.txt && grep -q 'trim (1)' ov.txt && grep -q 'align (2)' ov.txt; then _pass "overview lists steps with numbers (1,2)"
  else _fail "steps numbered list" "$(grep Steps: ov.txt)"; fi
  barA=$(grep -E '\(time\) .* 67%' ov.txt | grep -oE '\|+' | head -1)
  if [ "${#barA}" -eq 20 ]; then _pass "bar fill is percent-scaled (67% -> 20 of 30 cells)"
  else _fail "bar fill percent-scaled" "67% bar = ${#barA} cells"; fi
  # Overall (and all total rows) now carry a bar, not just numbers
  if grep -E 'Overall +\(time\)' ov.txt | grep -q '\['; then _pass "Overall total row has a bar (not numbers-only)"
  else _fail "Overall row has a bar" "$(grep -E 'Overall' ov.txt | head -1)"; fi

  # ---- low-efficiency RED override + effLowPct knob ---------------------------
  # No row here is <=25% (min is sampleB mem 50%), so the default threshold paints
  # NO red bars; raising effLowPct to 60 turns sampleB's 50-54% bars red.
  checkRunGrid --efficiency --color --log smartSlurmLog > ovcol.txt 2>&1
  redDefault=$(grep -c $'\033\[38;5;196m' ovcol.txt)
  if [ "${redDefault:-0}" -eq 0 ]; then _pass "no red bars when all rows are above the low-eff threshold (default 25%)"
  else _fail "red only below threshold" "$redDefault red at default"; fi
  effLowPct=60 checkRunGrid --efficiency --color --log smartSlurmLog > ovcol60.txt 2>&1
  red60=$(grep -c $'\033\[38;5;196m' ovcol60.txt)
  if [ "${red60:-0}" -ge 2 ]; then _pass "effLowPct knob raises the cutoff -> low rows (<=60%) turn red ($red60 bars)"
  else _fail "effLowPct threshold knob" "got $red60 red bars at effLowPct=60"; fi

  # ---- SAMPLE view (letter) : that sample's per-step breakdown ----------------
  checkRunGrid --efficiency --eff-view sample --eff-select A --no-color --log smartSlurmLog > sv.txt 2>&1
  assert_contains sv.txt "Sample view: showing sampleA" "sample view header"
  if grep -q 'trim (1)' sv.txt && grep -q 'align (2)' sv.txt; then _pass "sample view breaks down into per-step rows (trim, align)"
  else _fail "sample view per-step rows" "$(cat sv.txt)"; fi

  # ---- STEP view (number) : that step across all samples ----------------------
  checkRunGrid --efficiency --eff-view step --eff-select 1 --no-color --log smartSlurmLog > stv.txt 2>&1
  assert_contains stv.txt "Step view: showing trim (1)" "step view header"
  if grep -qE '\(A\).*\(time\)' stv.txt && grep -qE '\(B\).*\(time\)' stv.txt; then _pass "step view shows the step across all samples (A,B)"
  else _fail "step view per-sample rows" "$(cat stv.txt)"; fi
  if grep -qE 'trim +\(time\) .* 60%' stv.txt; then _pass "step total time% across samples (trim 60%)"
  else _fail "step total aggregate" "$(grep -E 'trim +\(time\)' stv.txt)"; fi

  # ---- toggles + colors ------------------------------------------------------
  checkRunGrid --efficiency --eff-show time --no-color --log smartSlurmLog > tonly.txt 2>&1
  if grep -q '(mem)' tonly.txt; then _fail "--eff-show time hides mem rows" "mem row present"
  else _pass "--eff-show time shows only time rows"; fi
  checkRunGrid --efficiency --color --log smartSlurmLog > ovc.txt 2>&1
  # time = 256-cyan (38;5;45), mem = 256-magenta (38;5;201) -- deliberately NOT
  # green, so the time/mem colors never read as a red/green pass-fail scale next to
  # the 256-red (38;5;196) low-eff bars.
  if grep -q $'\033\[38;5;45m' ovc.txt && grep -q $'\033\[38;5;201m' ovc.txt && ! grep -q $'\033\[32m' ovc.txt; then _pass "time bars cyan + mem bars magenta (not green) under --color"
  else _fail "distinct bar colors under --color" "$(grep -oaE $'\033\[3[0-9]m' ovc.txt | sort -u | tr '\n' ' ')"; fi
  # the SAME color scheme must reach the submenus (one scheme at all levels)
  checkRunGrid --efficiency --eff-view sample --eff-select A --color --log smartSlurmLog > svc.txt 2>&1
  if grep -q $'\033\[38;5;45m' svc.txt && grep -q $'\033\[38;5;201m' svc.txt; then _pass "cyan/magenta scheme reaches the sample submenu too"
  else _fail "submenu colors" "$(grep -oaE $'\033\[3[0-9]m' svc.txt | sort -u | tr '\n' ' ')"; fi
  if grep -q $'\033\[' ov.txt; then _fail "--no-color strips ANSI" "escapes present"
  else _pass "--no-color renders plain (no ANSI)"; fi

  # ---- terabyte-scale memory rolls M->G->T (column never overflows) -----------
  # Own log dir so the big numbers don't perturb the aggregates asserted above.
  # alloc 2 TiB = 2097152 M, used 1.5 TiB = 1572864 M -> 1.5T / 2.0T, 75%.
  mkdir -p smartSlurmLogTB
  printf 'epoch\tflag\tsample\tprogram\tref\tstatus\tallocMem\tusedMem\tallocTime\tusedTime\tinputs\n' > smartSlurmLogTB/efficiency.tsv
  printf '1\t1.0.huge.big\tbig\thuge\tnone\tCOMPLETED\t2097152\t1572864\t600\t420\t,big.fq\n' >> smartSlurmLogTB/efficiency.tsv
  checkRunGrid --efficiency --no-color --log smartSlurmLogTB > tb.txt 2>&1
  if grep -qE '2\.0T' tb.txt && grep -qE '1\.5T' tb.txt; then _pass "memory >=1 TiB rolls G->T (1.5T / 2.0T)"
  else _fail "memory rolls to T past 1024G" "$(grep -E '\(mem\)' tb.txt | head -2)"; fi
  colW=$(grep -E '\(mem\)' tb.txt | head -1 | awk '{print length($NF)}')
  if [ "${colW:-99}" -le 7 ]; then _pass "T-scaled mem value stays within the %7s column (${colW} chars)"
  else _fail "T value fits the column" "alloc field = ${colW} chars"; fi

  # ---- long step/sample labels keep the bars aligned (adaptive label width) ---
  # A step name wider than the 14-col default must widen the whole column, not
  # shove its own bar right. Seed one long-named step beside a short one under a
  # single sample and assert every bar's '[' lands in the same byte column.
  mkdir -p smartSlurmLogLW
  printf 'epoch\tflag\tsample\tprogram\tref\tstatus\tallocMemMB\tusedMemMB\tallocTimeMin\tusedTimeMin\tinputs\tdate\n' > smartSlurmLogLW/efficiency.tsv
  printf '1\t2.0.umiCutAdaptSeqtk.liver1\tliver1\tumiCutAdaptSeqtk\tnone\tCOMPLETED\t16384\t3000\t600\t100\t,x\n' >> smartSlurmLogLW/efficiency.tsv
  printf '2\t3.0.bowtie2.liver1\tliver1\tbowtie2\tnone\tCOMPLETED\t16384\t8000\t600\t310\t,y\n'                   >> smartSlurmLogLW/efficiency.tsv
  checkRunGrid --efficiency --eff-view sample --eff-select A --no-color --log smartSlurmLogLW > lw.txt 2>&1
  bars=$(awk '/\[/{print index($0,"[")}' lw.txt | sort -u | wc -l)
  if [ "${bars:-0}" -eq 1 ]; then _pass "long labels keep every bar aligned (single '[' column)"
  else _fail "bar alignment with long labels" "$(grep -n '\[' lw.txt)"; fi

  # ---- long SAMPLE names are shortened (…) but stay recoverable on drill-in ----
  # Real fastq paths made the overview a wall of text. The view must strip the
  # common path + Illumina _S.._L.. suffix, window to an … label, hide full paths
  # in the overview, and still show the full name in the sample-view header.
  mkdir -p smartSlurmLogSN
  printf 'epoch\tflag\tsample\tprogram\tref\tstatus\tallocMemMB\tusedMemMB\tallocTimeMin\tusedTimeMin\tinputs\tdate\n' > smartSlurmLogSN/efficiency.tsv
  # the sample is derived from the flag loper (rAP mangles the path, / -> _), so put
  # the long name there; $3 mirrors it (the view reads the flag, not $3).
  sn1='_n_data1_cores_ntc_FC_09468_SRN00333229_Durbin001-PR1800-5H7-DMSO-1_S17_L005_R1_001'
  sn2='_n_data1_cores_ntc_FC_09468_SRN00333230_Durbin001-PR1801-5H7-DMSO-2_S18_L005_R1_001'
  printf '1\t2.0.trim.%s\t%s\ttrim\tnone\tCOMPLETED\t16384\t8000\t600\t300\t,x\n' "$sn1" "$sn1" >> smartSlurmLogSN/efficiency.tsv
  printf '2\t2.0.trim.%s\t%s\ttrim\tnone\tCOMPLETED\t16384\t8000\t600\t300\t,y\n' "$sn2" "$sn2" >> smartSlurmLogSN/efficiency.tsv
  CHECKRUN_SAMPLE_W=24 checkRunGrid --efficiency --no-color --log smartSlurmLogSN > sn.txt 2>&1
  if grep -q '…' sn.txt; then _pass "long sample names are shortened with an … ellipsis"
  else _fail "long sample name shortening" "$(cat sn.txt)"; fi
  if grep -q '_n_data1_cores' sn.txt; then _fail "overview hides the shared path prefix" "prefix present in overview"
  else _pass "overview no longer dumps the full name"; fi
  if grep -q 'S17_L005' sn.txt; then _fail "sequencing suffix stripped" "found _S17_L005 in label"
  else _pass "Illumina _S.._L.. suffix stripped from labels"; fi
  # right-anchor (default): the discriminating TAIL is kept, not the shared head noise
  if grep -q 'DMSO-1' sn.txt && grep -q 'DMSO-2' sn.txt; then _pass "right-anchor keeps the discriminating tail (DMSO-1/DMSO-2)"
  else _fail "right-anchor keeps tail" "$(grep -E '…' sn.txt)"; fi
  # the anchor is overridable: left keeps the head instead
  CHECKRUN_SAMPLE_W=24 CHECKRUN_LABEL_ANCHOR=left checkRunGrid --efficiency --no-color --log smartSlurmLogSN > snl.txt 2>&1
  if grep -qE 'SRN0033322[0-9].*…' snl.txt; then _pass "CHECKRUN_LABEL_ANCHOR=left keeps the head instead"
  else _fail "left-anchor override" "$(grep -E '…' snl.txt)"; fi
  checkRunGrid --efficiency --eff-view sample --eff-select A --no-color --log smartSlurmLogSN > sna.txt 2>&1
  if grep -q 'Sample view: showing _n_data1' sna.txt && grep -q 'Durbin001-PR1800-5H7-DMSO-1' sna.txt; then _pass "drill-in header carries the full sample name"
  else _fail "full name on drill-in" "$(head -1 sna.txt)"; fi

  # ---- common-prefix drop snaps to a token boundary (no mid-token slice) -------
  # char-level LCP of these is "run_093" (mid-token); token-aware trimming must
  # keep the whole flowcell token (09328/09390), never "328"/"390".
  mkdir -p smartSlurmLogTK
  printf 'epoch\tflag\tsample\tprogram\tref\tstatus\tallocMemMB\tusedMemMB\tallocTimeMin\tusedTimeMin\tinputs\tdate\n' > smartSlurmLogTK/efficiency.tsv
  printf '1\t2.0.trim.run_09328_alpha\trun_09328_alpha\ttrim\tnone\tCOMPLETED\t1024\t512\t60\t30\t,x\n' >> smartSlurmLogTK/efficiency.tsv
  printf '2\t2.0.trim.run_09390_beta\trun_09390_beta\ttrim\tnone\tCOMPLETED\t1024\t512\t60\t30\t,y\n'  >> smartSlurmLogTK/efficiency.tsv
  checkRunGrid --efficiency --no-color --log smartSlurmLogTK > tk.txt 2>&1
  if grep -q '09328_alpha' tk.txt && ! grep -qE '(^|[^0-9])328_alpha' tk.txt; then _pass "common-prefix drop snaps to a token boundary (keeps 09328, not 328)"
  else _fail "token-boundary prefix drop" "$(grep -E '0?9?328' tk.txt)"; fi

  # ---- base-26 row letters roll A..Z, AA.. (not into punctuation past 'Z') ------
  mkdir -p smartSlurmLogZZ
  printf 'epoch\tflag\tsample\tprogram\tref\tstatus\tallocMemMB\tusedMemMB\tallocTimeMin\tusedTimeMin\tinputs\tdate\n' > smartSlurmLogZZ/efficiency.tsv
  for i in $(seq -w 1 27); do printf '%s\t2.0.trim.s%s\tsmp%s\ttrim\tnone\tCOMPLETED\t1024\t512\t60\t30\t,x\n' "$i" "$i" "$i" >> smartSlurmLogZZ/efficiency.tsv; done
  checkRunGrid --efficiency --no-color --log smartSlurmLogZZ > zz.txt 2>&1
  if grep -qE '\(AA\)' zz.txt; then _pass "27th sample gets base-26 letter (AA)"
  else _fail "base-26 letter rollover" "$(grep -oE '\([^)]+\)' zz.txt | tail -3 | tr '\n' ' ')"; fi
  if grep -qF '([' zz.txt; then _fail "no punctuation letters past Z" "found '([' (old chr(64+n) overflow)"
  else _pass "no punctuation row letters past Z (chr(64+n) overflow fixed)"; fi

  # ---- experiment-level steps (no per-sample loper) are NOT listed as samples ----
  # A step whose flag has no sample field (e.g. 9.0.multiqc) spans all samples; it
  # must appear in the STEP list, never as a per-sample row (regression: such steps
  # leaked into the sample axis via cleanUp's program-name fallback in $3).
  mkdir -p smartSlurmLogEX
  printf 'epoch\tflag\tsample\tprogram\tref\tstatus\tallocMemMB\tusedMemMB\tallocTimeMin\tusedTimeMin\tinputs\tdate\n' > smartSlurmLogEX/efficiency.tsv
  printf '1\t2.0.trim.s1\ts1\ttrim\tnone\tCOMPLETED\t1024\t512\t60\t30\t,x\n'          >> smartSlurmLogEX/efficiency.tsv
  printf '2\t2.0.trim.s2\ts2\ttrim\tnone\tCOMPLETED\t1024\t512\t60\t30\t,y\n'          >> smartSlurmLogEX/efficiency.tsv
  printf '3\t9.0.multiqc\tmultiqc\tmultiqc\tnone\tCOMPLETED\t2048\t1024\t120\t60\t,none\n' >> smartSlurmLogEX/efficiency.tsv
  checkRunGrid --efficiency --no-color --log smartSlurmLogEX > ex.txt 2>&1
  nsmp=$(grep -cE '\([A-Z]+\).*\(time\)' ex.txt)
  if [ "${nsmp:-0}" -eq 2 ] && ! grep -qE 'multiqc +\([A-Z]+\)' ex.txt; then _pass "experiment-level step (multiqc) is not listed as a sample (2 real samples)"
  else _fail "experiment-level excluded from samples" "$(grep -E '\(time\)|Samples:|Steps:' ex.txt)"; fi
  if grep -qE 'Steps:.*multiqc \(9\)' ex.txt; then _pass "experiment-level step still shows in the step list"
  else _fail "experiment-level step in step list" "$(grep Steps: ex.txt)"; fi
  # experiment-level steps get their own overview section so they are still viewable
  if grep -q 'Experiment-level (across all samples):' ex.txt && grep -qE 'multiqc \(9\).*\(time\)' ex.txt; then _pass "experiment-level steps shown in their own overview section"
  else _fail "experiment-level overview section" "$(cat ex.txt)"; fi

  # ---- a step that ran with NO recorded allocation shows n/a, not a red 0% -------
  # (e.g. a summary/multistats job that never went through resource estimation:
  # alloc 0, used 15m). It must appear, with n/a efficiency and a grey (not red) bar.
  mkdir -p smartSlurmLogNA
  printf 'epoch\tflag\tsample\tprogram\tref\tstatus\tallocMemMB\tusedMemMB\tallocTimeMin\tusedTimeMin\tinputs\tdate\n' > smartSlurmLogNA/efficiency.tsv
  printf '1\t2.0.trim.s1\ts1\ttrim\tnone\tCOMPLETED\t1024\t512\t60\t30\t,x\n'              >> smartSlurmLogNA/efficiency.tsv
  printf '2\t11.0.multistats\tmultistats\tmultistats\tnone\tCOMPLETED\t0\t0\t0\t15\t,none\n' >> smartSlurmLogNA/efficiency.tsv
  checkRunGrid --efficiency --no-color --log smartSlurmLogNA > na.txt 2>&1
  if grep -qE 'multistats \(11\).*n/a' na.txt; then _pass "unallocated step (multistats) shows n/a efficiency, and is still listed"
  else _fail "n/a for unallocated step" "$(grep -E 'multistats|n/a' na.txt)"; fi
  # its used time is still shown (it ran 15m), so the row is informative
  if grep -qE 'multistats \(11\).*15m' na.txt; then _pass "unallocated step still shows its actual usage (15m)"
  else _fail "usage shown for unallocated step" "$(grep multistats na.txt)"; fi
  # the n/a bar must NOT be red (grey); red would read as a 0% waste alert
  checkRunGrid --efficiency --color --log smartSlurmLogNA > nac.txt 2>&1
  if grep -E 'multistats' nac.txt | grep -q $'\033\[38;5;196m'; then _fail "n/a bar is not red" "red escape on multistats row"
  else _pass "n/a bar is grey, not the red low-eff alert"; fi

  # ---- step list is in ascending NUMERIC order (10 after 2, not after 1) --------
  mkdir -p smartSlurmLogORD
  printf 'epoch\tflag\tsample\tprogram\tref\tstatus\tallocMemMB\tusedMemMB\tallocTimeMin\tusedTimeMin\tinputs\tdate\n' > smartSlurmLogORD/efficiency.tsv
  printf '1\t10.0.heatmap10.s1\ts1\theatmap10\tnone\tCOMPLETED\t1024\t512\t60\t30\t,x\n' >> smartSlurmLogORD/efficiency.tsv
  printf '2\t2.0.trim.s1\ts1\ttrim\tnone\tCOMPLETED\t1024\t512\t60\t30\t,x\n'           >> smartSlurmLogORD/efficiency.tsv
  printf '3\t3.0.bowtie2.s1\ts1\tbowtie2\tnone\tCOMPLETED\t1024\t512\t60\t30\t,x\n'     >> smartSlurmLogORD/efficiency.tsv
  checkRunGrid --efficiency --no-color --log smartSlurmLogORD > ord.txt 2>&1
  # the Steps legend must read (2) ... (3) ... (10), despite being seeded 10,2,3
  if grep -qE 'Steps:.*\(2\).*\(3\).*\(10\)' ord.txt; then _pass "step list sorted numerically (2,3,10 not 10,2,3)"
  else _fail "numeric step order" "$(grep Steps: ord.txt)"; fi
  echo
}

# ===========================================================================
_want email_consolidation && {
  _begin email_consolidation
  # --- fully successful run: consolidated, text-only emails, no pngs ----------
  ss_new_workspace "$WORK/email_consolidation"
  cp "$HERE"/fixtures/clean_2step/* .; chmod +x ./*.sh 2>/dev/null || true
  ss_run_pipeline "./pipeline.sh 777" "sbatch -p short -c 1 --mem 4G -t 50:0"
  mlog="$MOCK_SLURM_STATE/mail.log"
  assert_file "$mlog" "mail mock recorded emails"
  if ls smartSlurmLog/*.png >/dev/null 2>&1; then _fail "no per-job/barchart pngs written" "$(ls smartSlurmLog/*.png 2>/dev/null)"
  else _pass "no per-job/barchart pngs written (plot generation removed)"; fi
  # the per-job/barchart gnuplot CSVs (<jobid>.dataMem.csv, .dataTime.csv,
  # .output.csv) fed the removed pngs; they must no longer be written either.
  pc=$(ls smartSlurmLog/*.dataMem.csv smartSlurmLog/*.dataTime.csv smartSlurmLog/*.output.csv 2>/dev/null)
  if [ -n "$pc" ]; then _fail "no plot CSVs written" "$pc"
  else _pass "no plot CSVs (dataMem/dataTime/output) written"; fi
  # 0-byte "<flag>.out.<jobid>" job-id marker files are no longer created: they
  # were only read cosmetically by checkRun; the job id lives in allJobs.txt and
  # the .out record. Real job logs are "<flag>.out" (no trailing .<jobid>).
  jm=$(ls smartSlurmLog/*.out.[0-9]* 2>/dev/null)
  if [ -n "$jm" ]; then _fail "no 0-byte .out.<jobid> job-id markers written" "$jm"
  else _pass "no .out.<jobid> job-id marker files written"; fi
  if grep -q 'att=1' "$mlog" 2>/dev/null; then _fail "all emails text-only (no attachments)" "found att=1"
  else _pass "all emails are text-only (no -a attachments)"; fi
  if grep -q ' to=ld32' "$mlog" 2>/dev/null; then _fail "no hard-coded dev recipient" "found ld32"
  else _pass "no hard-coded dev (ld32) recipient"; fi
  effN=$(grep -c 'run complete: resource efficiency' "$mlog" 2>/dev/null)
  if [ "${effN:-0}" -eq 1 ]; then _pass "exactly one per-run efficiency email"
  else _fail "one per-run efficiency email" "got ${effN:-0}"; fi
  assert_contains "$mlog" "Resource efficiency -- overview" "efficiency email carries the per-sample overview table"
  assert_contains "$mlog" "press e"                          "efficiency email points to checkRun"
  cmpN=$(grep -c 'done: sample' "$mlog" 2>/dev/null)
  if [ "${cmpN:-0}" -ge 1 ]; then _pass "per-sample completion email(s) sent ($cmpN)"
  else _fail "per-sample completion email sent" "got ${cmpN:-0}"; fi

  # --- a sample breaks: one warning per broken sample, no "run complete" email -
  ss_new_workspace "$WORK/email_consolidation_fail"
  cp "$HERE"/fixtures/clean_2step/* .; chmod +x ./*.sh 2>/dev/null || true
  MOCK_FAIL_JOB=mergeNumber MOCK_FAIL_MODE=oom \
    ss_run_pipeline "./pipeline.sh 777" "sbatch -p short -c 1 --mem 4G -t 50:0"
  mlog2="$MOCK_SLURM_STATE/mail.log"
  failN=$(grep -c 'FAIL: sample' "$mlog2" 2>/dev/null)
  if [ "${failN:-0}" -ge 1 ]; then _pass "broken sample triggers a text FAIL warning ($failN)"
  else _fail "break warning email sent" "got ${failN:-0}"; fi
  assert_contains "$mlog2" "Re-run the same command" "break warning tells the user to re-run"
  if grep -q 'run complete: resource efficiency' "$mlog2" 2>/dev/null; then _fail "no run-complete email when a sample broke" "found one"
  else _pass "no per-run efficiency email when the run did not fully succeed"; fi
  # the per-run "SmartSlurm Saved <$>" figure is computed from the plot CSV, which
  # is gated off with the pngs; with plots off it must be skipped cleanly, never
  # emitted with a blank amount ("Saved  by") or leaking awk/bc errors into mail.
  if grep -qE 'Saved +(by|$)|Saved  ' "$mlog2" 2>/dev/null; then _fail "no blank savings line when plots are off" "$(grep -n 'Saved' "$mlog2")"
  else _pass "no blank/garbled per-run savings line when plots are off"; fi
  echo
}

# ===========================================================================
ss_report
[ "$FAILS" -eq 0 ] && exit 0 || exit 1
