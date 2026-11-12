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
ALL_BLOCKS="clean_2step wf_afterok wf_scriptargs wf_depgrammar wf_arrays wf_wait wf_scancel wf_inject wf_tiered wf_pipeline wf_isolation wf_jobseq wf_exit0_empty ma_extract_dag ma_emit_rap ma_parse_dag ma_compare ma_differential nested_fanout ragged_axes fail_exit0_empty known_broken_swallowed_done \
loop_detection oom_injection tool_stderr diagnose fate_reconcile hardkill oomoot_escalate rerun end_block verbosity scaling_linear \
noise_fit scaling_exp contaminate_exit0 contaminate_exclude contaminate_prune firstpass_rate \
failure_floor manual_rerun_floor constant_model negslope_robust \
setup_rc_block activate_path activate_prompt deactivate_surgical conda_coexist \
uninstall_precise unexport_deprecated config_overlay config_guided config_migrate \
firstpass_fullpath extrapolate_linear extrapolate_superlinear extrapolate_lowdata \
resubmit_convergence manual_resubmit_sequence \
efficiency_file checkrun_efficiency email_consolidation fitqc migrate_v1 jr_drift_guard cross_pipeline legacy_pooling"

# Named groups: a convenience alias that expands to a block list, so related
# blocks that do NOT share a common name substring can be run in one go (e.g.
# the activation feature spans activate/deactivate/uninstall/conda/setup). A
# non-group argument is returned unchanged and matched as a substring. Keep the
# membership lists here in sync with ALL_BLOCKS above.
GROUP_NAMES="estimator activation efficiency workflow p6"
ss_group() {
  case "$1" in
    estimator)  echo "scaling_linear noise_fit scaling_exp contaminate_exit0 contaminate_exclude contaminate_prune firstpass_rate failure_floor manual_rerun_floor constant_model negslope_robust firstpass_fullpath extrapolate_linear extrapolate_superlinear extrapolate_lowdata resubmit_convergence manual_resubmit_sequence fitqc migrate_v1 jr_drift_guard cross_pipeline legacy_pooling";;
    activation) echo "setup_rc_block activate_path activate_prompt deactivate_surgical conda_coexist uninstall_precise unexport_deprecated";;
    efficiency) echo "efficiency_file checkrun_efficiency email_consolidation";;
    workflow)   echo "wf_afterok wf_scriptargs wf_depgrammar wf_arrays wf_wait wf_scancel wf_inject wf_tiered wf_pipeline wf_isolation wf_jobseq wf_exit0_empty ma_extract_dag ma_emit_rap ma_parse_dag ma_compare ma_differential";;
    p6)         echo "p6_signal_sequence p6_determinism cleanup_reentrancy";;
    *)          echo "$1";;
  esac
}

ss_purpose() {
  case "$1" in
    clean_2step)        echo "Baseline: convert a 2-step pipeline, submit under mock Slurm, and wire afterok dependencies end-to-end";;
    wf_afterok)         echo "Migration Assistant: an UNMODIFIED hand-rolled SLURM afterok script runs invisibly end-to-end via register-don't-run + topological drain";;
    wf_scriptargs)      echo "Migration Assistant: a worker invoked as 'sbatch [opts] worker.sh --flag value …' keeps its --flag args verbatim (real sbatch passes everything after the script to the script) — registration + drained worker both receive them";;
    wf_depgrammar)      echo "Migration Assistant: the full sbatch --dependency grammar — after/afterany/afterok/afternotok/singleton, AND(,)/OR(?), +time, and kill-vs-hold on an unsatisfiable dep";;
    wf_arrays)          echo "Migration Assistant: job arrays (--array) — per-task execution, whole-array afterok, and aftercorr element-wise correspondence (incl. a failing task)";;
    wf_wait)            echo "Migration Assistant: sbatch --wait runs the job (and its deps) inline and returns the job's exit code";;
    wf_scancel)         echo "Migration Assistant: scancel cancels a pending job (orphaning its afterok downstream); squeue reports live pending state";;
    wf_inject)          echo "Migration Assistant: resource-sim fault injection (MOCK_FAIL_JOB/MODE) forces a workflow step to fail, exercising afterok orphaning + afternotok cleanup";;
    wf_tiered)          echo "Migration Assistant: tiered execution — run real payloads, or stub a step whose tool/input is unavailable (explicit SS_WF_STUB, or auto SS_WF_STUB_MISSING)";;
    wf_pipeline)        echo "Migration Assistant: realistic end-to-end hand-rolled SLURM pipeline (per-sample afterok chains + multi-id afterok fan-in + report) with verified outputs and fan-in failure orphaning";;
    wf_isolation)       echo "Migration Assistant: no crosstalk between the sbatch (wf-mode) and rAP arms — an rAP pipeline still executes while SS_WF_MODE is set, in a separate workspace (differential-test isolation)";;
    wf_jobseq)          echo "Migration Assistant: the mock sbatch job-id counter is atomic — many CONCURRENT submissions against one mock state all get unique ids (no read-modify-write collisions)";;
    wf_exit0_empty)     echo "Migration Assistant: silent data failure — a step exits 0 without its output, so the afterok consumer RUNS (not orphaned) and then fails on the missing input (contrast with a hard failure)";;
    ma_extract_dag)     echo "Migration Assistant skill: extract_dag.sh recovers a hand-rolled SLURM workflow's DAG (nodes + name-resolved dependency edges) by reading the register-don't-run state";;
    ma_emit_rap)        echo "Migration Assistant skill: emit_rap.sh renders a logical plan into an annotated rAP script; round-trip — the emitted pipeline runs under runAsPipeline and reproduces the original's outputs";;
    ma_parse_dag)       echo "Migration Assistant skill: parse_rap_dag.sh reads an rAP driver's #@ step DAG back (dependency step-numbers resolved to names), round-tripping emit_rap and resolving dot-joined multi-deps";;
    ma_compare)         echo "Migration Assistant skill: compare_outputs.sh diffs two runs' authoritative outputs under an identity contract — matches, normalizes volatile fields, flags real diffs and present-in-one-only";;
    ma_differential)    echo "Migration Assistant skill: full Scenario-2 — run the ORIGINAL (wf-mode) and its rAP PORT in isolated workspaces and confirm identical authoritative outputs under a contract";;
    nested_fanout)      echo "Nested per-sample/per-chrom loops — a fan-out step (one grid cell backs multiple jobs) converts to valid bash and every job completes";;
    ragged_axes)        echo "Two steps whose loop axes are disjoint (sample ids vs group ids) — ragged addressing converts and runs end-to-end";;
    loop_detection)     echo "Documents the textual for*/done* loop foot-gun: a one-line for...done and a for-prefixed command corrupt the next step's flag";;
    fail_exit0_empty)   echo "Producer exits 0 but writes nothing and the consumer fails — documents the current symptom-vs-cause misattribution";;
    known_broken_swallowed_done) echo "A blank-line block right before a loop 'done' makes the converter emit invalid bash (the M4 block-closure foot-gun)";;
    oom_injection)      echo "Injected OOM/timeout/nonzero jobs are attributed correctly by cleanUp (no spurious .success)";;
    tool_stderr)        echo "Tool chatter on stderr is bracketed by sentinels, not misread as a failure";;
    diagnose)           echo "diagnose.sh classifies failure state/cause/action and emits the expected JSON contract";;
    fate_reconcile)     echo "Hard-killed (OOM/timeout) jobs that bypass cleanUp are flagged unknown by the grid, resolved from sacct via jobFate (marker + report) only on an explicit reconcile, kept read-only by default, left unmarked when CANCELLED (orphan-safe), and degrade gracefully on sacct lag";;
    hardkill)           echo "Hard-kill an arbitrary pipeline step mid-run (MOCK_FAIL_MODE=sigkill/sigterm): the whole job process group dies, so cleanUp's EXIT trap never runs -- a faithful cgroup-OOM/node-kill model (sacct OUT_OF_MEMORY, no markers), the substrate for cleanUp signal/trap dev";;
    oomoot_escalate)    echo "OOMOOT: a job that hit BOTH the memory and time limit (sacct TIMEOUT + OUT_OF_MEMORY) -> cleanUp classifies OOMOOT and escalates mem AND time together (the dual-escalation branch)";;
    p6_signal_sequence) echo "P6 (opt-in via _optin; never in the default suite): deterministic signal/trap SEQUENCE matrix over a self-contained MODEL of the proposed rAP job body (tool as '& wait', run-once-guarded cleanup trapped on EXIT/TERM/USR1). Barrier-controlled (no wall-clock race): S0/S0b back-compat (clean success + nonzero exit identical to today's EXIT-only path), S1 TERM mid-tool, S2 USR1 grace-signal mid-tool, S6 post-completion no-op -- each asserting cleanup body runs EXACTLY once. Executable spec for the parked cleanUp change; touches no bin/";;
    p6_determinism)     echo "P6 (opt-in via _optin; never in the default suite): the anti-flake guarantee -- re-run the barrier-controlled TERM/USR1/no-signal cases 20x each and assert ZERO variance (one distinct outcome per case). Proves the timing harness is reproducible before any cleanUp/ssbatch edit; touches no bin/";;
    cleanup_reentrancy) echo "P6 SCAFFOLD (opt-in via _optin; never in the default suite; expected to FAIL until cleanUp gains a re-entrancy guard): invoking the REAL cleanUp twice on one job must not double-append jobRecords or double-requeue";;
    rerun)              echo "Rerun menu (failed / all / from-step:N) materializes markers and honors the reachable guard";;
    end_block)          echo "'#@ … #@end' scripts-on-disk: closure validation, --lint, legacy back-compat, and in-block multivar escaping";;
    verbosity)          echo "Stdout verbosity: default A+ per-sample rows vs -v firehose vs -q silence";;
    scaling_linear)     echo "Linear mem/time vs input size — the real gnuplot fit recovers slope/intercept and prediction tracks the line";;
    noise_fit)          echo "~8% stochastic noise — the fit still recovers the slope; STDFIT>0 and RSquare are computed";;
    scaling_exp)        echo "Exponential truth vs SmartSlurm's linear fit — the model mismatch under-provisions at large input sizes";;
    contaminate_exit0)  echo "An exit-0-empty / mis-sized COMPLETED row drags the fitted slope downward (data contamination)";;
    contaminate_exclude) echo "OOM/failed rows are excluded from the fit — only COMPLETED rows are used";;
    contaminate_prune)  echo "An OOM storm prunes the offending program's failed and failure-disproven rows, keeping its validated history and unrelated programs";;
    fitqc)              echo "fitQC: triage flags the planted outlier group, detail ranks it worst, --remove/--auto clean it (R2 jumps), groups stay per pipeline, high-leverage extremes are never auto-removed";;
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
    config_overlay)     echo "smartslurm config: overlay model — set/get/unset write only ~/.smartSlurm/config/user.conf (never the repo file, Q1); non-overridden repo defaults still load (Q2); schema validation rejects bad types; unknown keys flagged custom";;
    config_guided)      echo "smartslurm config: guided editors drive real behavior — refgroup add changes what jrRefClass pools; partition sets the partition<N>TimeLimit, resolving a partition by name";;
    config_migrate)     echo "smartslurm config migrate: a legacy full ~/.smartSlurm/config/config.txt is converted to the overlay — only keys differing from repo defaults are kept in user.conf, the legacy file is backed up";;
    firstpass_fullpath) echo "First-pass coverage through the FULL estimateResource path (extraMem + failure floor), for mem AND time (not just the calculateMemTime margin)";;
    extrapolate_linear) echo "[Phase 1 baseline] Linear truth trained on small sizes, predicted at 2x-50x: with SCount>9 it extrapolates the line (not outOfRange); reports first-pass per scale";;
    extrapolate_superlinear) echo "[Phase 1 baseline] Quadratic truth trained on small sizes: the straight-line fit UNDER-provisions at 10x-50x; records the shortfall ratio (the problem to fix)";;
    extrapolate_lowdata) echo "Thin data (<10 records) beyond the trained range uses a trend-scaled fallback (in-range estimate x size/Maximum), not a bare default";;
    resubmit_convergence) echo "Auto-resubmit: an under-provisioned job OOMs and is auto-requeued with monotonically escalating resources until it completes; verifies the ladder never regresses";;
    manual_resubmit_sequence) echo "Multiple sequential MANUAL resubmissions (each with its own auto-requeue), bootstrapping to completion; each manual start exceeds the prior proven-insufficient alloc";;
    efficiency_file)    echo "cleanUp writes one unconditional per-job row to a per-run efficiency.tsv (flag, sample, status, alloc/used mem+time) for checkRun to read";;
    checkrun_efficiency) echo "checkRunGrid --efficiency drill-down views (overview/sample/step): Overall + per-sample bars, per-step and per-sample breakdowns, eff/used/alloc columns, time/mem toggles, color";;
    email_consolidation) echo "cleanUp sends consolidated TEXT-ONLY mail (no pngs): per-sample break/complete alerts + one per-run efficiency email; no per-job spam, no hard-coded dev recipient";;
    migrate_v1)         echo "jobRecords v1->v2 auto-migration: headerless/20-col/21-col legacy, MIXED, both HYBRID orderings (v1 appended under a v2 header) and headerless-v2 stores migrate with zero loss + backup, survive the real cleanUp prune (no self-destruct), jrRepair diverts only junk; plus a property/fuzz + concurrency harness over randomized shapes (the ld32 data-loss regression)";;
    jr_drift_guard)     echo "Schema-drift guard: no SmartSlurm tool (outside the jobRecordSchema.sh API + jrRepair) writes a jobRecord header literal, redirects to a store, or reads it by raw column number -- so the 18-vs-21-col drift and raw-append hybrid cannot reappear";;
    cross_pipeline)     echo "Cross-pipeline isolation (grouping key = pipeline/step/ref): a colliding (inputSize,memUsed) row in another pipeline is NOT deduped away; a fail-prune in one pipeline does NOT delete another's COMPLETED history (and still prunes its own); standalone unknown never pools into/broadens a pipeline fit -- the ld32 cross-pipeline data-loss regression";;
    legacy_pooling)     echo "Standalone ssbatch estimation uses migrated jR1 (legacy) history: a pipeline=unknown job pools the {unknown, legacy} bucket (legacy grouped by step/ref as before the pipeline column), so old history still sizes jobs; unknown still never pools a NAMED pipeline; a fresh named pipeline still bootstraps from legacy -- the ld32 standalone-legacy regression";;
    *)                  echo "$1";;
  esac
}

source "$HERE/lib/harness.sh"

# Extract -jN / -j N (parallel block count; default 1 = serial) from the args;
# everything else is left as block-name / group filters.
PARJOBS=1; _rtargs=()
while [ $# -gt 0 ]; do
  case "$1" in
    -j)  if [[ "${2:-}" =~ ^[0-9]+$ ]]; then PARJOBS="$2"; shift 2; else PARJOBS=4; shift; fi;;
    -j*) _v="${1#-j}"; [[ "$_v" =~ ^[0-9]+$ ]] && PARJOBS="$_v"; shift;;
    *)   _rtargs+=("$1"); shift;;
  esac
done
set -- ${_rtargs[@]+"${_rtargs[@]}"}

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

# _optin -- stricter than _want for OPT-IN scaffolds (the P6 signal/trap blocks).
# These model or probe an UNBUILT cleanUp change, so some are expected-red by
# design; they must NEVER run in the default suite. Unlike _want, an EMPTY filter
# (a bare `./run_tests.sh` or bootstrap) selects NOTHING here -- a block runs only
# when its name is given explicitly on the command line -- so the default suite
# (serial OR -jN, since opt-in blocks are also absent from ALL_BLOCKS) stays green.
# Run them during P6 dev with `./run_tests.sh p6` or a single block name.
_optin() {
  [ -z "${FILTERS// /}" ] && return 1
  local b="$1" f
  for f in $FILTERS; do [[ "$b" == *"$f"* ]] && return 0; done
  return 1
}

# --- parallel dispatch (-jN) ------------------------------------------------
# With -jN (N>1), run each SELECTED block as its own isolated process and exit
# before the inline blocks below. Block-level only; each block stays internally
# sequential. ss_activate above did the one-time global setup (gnuplot,
# /bin/sh->bash); the workers skip it via SS_NO_INSTALL/SS_NO_SH_REPOINT.
if [ "${PARJOBS:-1}" -gt 1 ]; then
  _sel=""
  for _b in $ALL_BLOCKS; do _want "$_b" && _sel="$_sel $_b"; done
  [ -n "${_sel// /}" ] || { echo "no blocks match the given filter"; exit 0; }
  ss_parallel_run "$PARJOBS" "$HERE/$(basename "$0")" $_sel
  exit $?
fi

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
# Migration Assistant — thin end-to-end slice. Prove the invisible user-workflow
# harness: a hand-rolled SLURM script (naked sbatch + afterok, knowing nothing
# about SmartSlurm) runs UNMODIFIED via register-don't-run + topological drain.
_want wf_afterok && {
  _begin wf_afterok
  ss_new_workspace "$WORK/wf_afterok"
  cp "$HERE"/fixtures/wf_afterok/* .; chmod +x ./*.sh 2>/dev/null || true
  ss_wf_activate
  # Phase 1: run the user's submit script. In workflow mode each sbatch only
  # REGISTERS a job and returns — no payload runs yet.
  bash ./submit.sh > wf.log 2>&1; rc=$?
  assert_rc "$rc" 0 "submit script exits 0 (job-id capture + afterok submit work)"
  assert_contains "wf.log" "submitted stepA as"            "stepA job id captured from mock sbatch"
  wfdir="$MOCK_SLURM_STATE/wf"
  n=$(ls "$wfdir"/*.run 2>/dev/null | wc -l)
  assert_rc 2 "$n" "two jobs registered (stepA, stepB)"
  assert_no_file "marker.txt"  "register-don't-run: no payload executed at submit time"
  assert_no_file "fromB.txt"   "register-don't-run: downstream payload also deferred"
  # Phase 2: drain the DAG. afterok must run stepA before stepB.
  ss_wf_drain
  assert_file    "marker.txt"                              "stepA ran during drain"
  assert_file    "fromB.txt"                               "stepB ran during drain"
  assert_contains "fromB.txt" "A-was-here"                 "stepB saw stepA's output (afterok ordering honored)"
  assert_contains "fromB.txt" "B-done"                     "stepB completed"
  # execution order: stepA (lower id) ran before stepB in the drain log
  assert_contains "$wfdir/.order" "RAN"                    "drain recorded an execution order"
  order=$(awk '$2=="RAN"{print $1}' "$wfdir/.order" | tr '\n' ' ')
  first=${order%% *}
  assert_contains "$wfdir/$first.name" "stepA"             "the first job drained was stepA (topological order)"
  echo
}

# ===========================================================================
# Migration Assistant — a worker invoked as `sbatch [opts] worker.sh --flag val …`
# must keep its --flag arguments. Real sbatch passes everything after the script
# name to the script; the wf-mode mock used to parse those --flags as sbatch
# options and drop them, collapsing the values to positionals (recorded
# `bash worker.sh fileA S1 full`), which mangled DAG recovery for real pipelines
# whose workers take `--flag value` args (e.g. the ChIP-seq mapping worker).
# Regression: registration must preserve the flags verbatim, and the drained
# worker must actually receive them.
_want wf_scriptargs && {
  _begin wf_scriptargs
  ss_new_workspace "$WORK/wf_scriptargs"
  cat > worker.sh <<'WSH'
#!/bin/bash
# a realistic worker: takes --flag value args and records what it received
sn=; r1=; r2=; mode=
while [ "$1" != "" ]; do
  case $1 in
    --read1)      shift; r1=$1;;
    --read2)      shift; r2=$1;;
    --samplename) shift; sn=$1;;
    --mode)       shift; mode=$1;;
    *) echo "unknown option: $1" >&2; exit 3;;
  esac
  shift
done
printf 'sample=%s read1=%s read2=%s mode=%s\n' "$sn" "$r1" "$r2" "$mode" > "got_${sn:-NONE}.txt"
WSH
  chmod +x worker.sh
  ss_wf_activate
  # sbatch options come BEFORE the script; the worker's own --flag args come AFTER it
  jid=$(sbatch -p short -t 60 --mem 4G --parsable ./worker.sh --read1 fileA --read2 fileB --samplename S1 --mode full); rc=$?
  assert_rc "$rc" 0 "sbatch registered the worker and returned a parsable id"
  run="$MOCK_SLURM_STATE/wf/$jid.run"
  # registration fidelity: the whole --flag arg tail survives verbatim (the bug)
  assert_contains "$run" "worker.sh --read1 fileA --read2 fileB --samplename S1 --mode full" \
                         "registered payload preserves every worker --flag arg, in order"
  # the sbatch options before the script must NOT leak into the worker payload
  assert_lacks    "$run" "4G"  "sbatch options (-t 60 --mem 4G) do not leak into the payload"
  # end-to-end: on drain the worker parses its flags (not mangled to positionals)
  ss_wf_drain
  assert_file     "got_S1.txt"                                                  "worker ran during drain and used --samplename"
  assert_contains "got_S1.txt" "sample=S1 read1=fileA read2=fileB mode=full"    "drained worker received all --flag args correctly"
  echo
}

# ===========================================================================
# Migration Assistant — the full sbatch --dependency grammar. One unmodified
# submit script builds a matrix of dependency types/conjunctions; the drain must
# resolve each to the right outcome (ran / cancelled / held).
_want wf_depgrammar && {
  _begin wf_depgrammar
  ss_new_workspace "$WORK/wf_depgrammar"
  cp "$HERE"/fixtures/wf_depgrammar/* .; chmod +x ./*.sh 2>/dev/null || true
  ss_wf_activate
  bash ./submit.sh > wf.log 2>&1; rc=$?
  assert_rc "$rc" 0 "submit script exits 0"
  ss_wf_drain
  wfdir="$MOCK_SLURM_STATE/wf"
  # map each job name -> readable terminal verdict (from its .rc)
  : > states.txt
  for f in "$wfdir"/*.name; do
    jid=$(basename "$f" .name); nm=$(cat "$f"); r=$(cat "$wfdir/$jid.rc" 2>/dev/null)
    case "$r" in 0) v=ran;; cancelled) v=cancelled;; held) v=held;; *) v="exit$r";; esac
    printf '%s %s\n' "$nm" "$v" >> states.txt
  done
  assert_contains states.txt "A ran"        "base A succeeded"
  assert_contains states.txt "B exit1"      "base B ran and failed (exit 1)"
  assert_contains states.txt "C ran"        "afterany ran after B terminated (any state)"
  assert_contains states.txt "D ran"        "afternotok ran BECAUSE B failed"
  assert_contains states.txt "E cancelled"  "afternotok on a succeeded dep is cancelled (never satisfiable)"
  assert_no_file  "ran.E"                    "cancelled afternotok payload never executed"
  assert_contains states.txt "F ran"        "afterok:A:G (AND within a term) ran once both succeeded"
  assert_contains states.txt "H ran"        "afterok:A?afterok:B (OR) ran, satisfied by A despite B failing"
  assert_contains states.txt "I held"       "--kill-on-invalid-dep=no HOLDS the job (not cancelled)"
  assert_no_file  "ran.I"                    "held job never executed"
  assert_contains states.txt "J cancelled"  "default kill-on-invalid-dep cancels an unsatisfiable afterok"
  assert_contains states.txt "K ran"        "after:B+5 ran once B started/terminated (+time ignored)"
  assert_file     "ran.solo1"               "singleton: first 'solo' ran"
  assert_file     "ran.solo2"               "singleton: second 'solo' ran after the first (no deadlock)"
  echo
}

# ===========================================================================
# Migration Assistant — job arrays (--array) and array dependencies.
_want wf_arrays && {
  _begin wf_arrays
  ss_new_workspace "$WORK/wf_arrays"
  cp "$HERE"/fixtures/wf_arrays/* .; chmod +x ./*.sh 2>/dev/null || true
  ss_wf_activate
  bash ./submit.sh > wf.log 2>&1; rc=$?
  assert_rc "$rc" 0 "submit script exits 0"
  ss_wf_drain
  # success array: all three tasks ran, each with its own SLURM_ARRAY_TASK_ID
  n=$(ls arr.*.txt 2>/dev/null | wc -l); assert_rc 3 "$n" "--array=1-3 expanded and ran 3 tasks"
  assert_file "arr.2.txt"                   "array task 2 ran under its own task id"
  n=$(ls arr_*_*.out 2>/dev/null | wc -l); assert_rc 3 "$n" "per-task -o %A/%a produced one log per task"
  # afterok on the whole array waits for ALL tasks, then runs
  assert_file "postall.0.txt"               "afterok on an array runs after ALL tasks succeed"
  # aftercorr: each task ran (its corresponding array task succeeded)
  assert_file "corr.2.txt"                  "aftercorr task ran for its corresponding successful array task"
  # failing array (task 2 exits 1): tasks still ran (markers written before exit)
  assert_file "bad.2.txt"                   "failing array task ran (marker before nonzero exit)"
  # afterok on a FAILED array -> cancelled
  assert_no_file "postbad.0.txt"            "afterok on a failed array is cancelled"
  # aftercorr on the failed array: task 2 cancelled; tasks 1 & 3 ran
  assert_file "corrbad.1.txt"               "aftercorr task 1 ran (dep task 1 succeeded)"
  assert_file "corrbad.3.txt"               "aftercorr task 3 ran (dep task 3 succeeded)"
  assert_no_file "corrbad.2.txt"            "aftercorr task 2 cancelled (dep task 2 failed)"
  echo
}

# ===========================================================================
# Migration Assistant — sbatch --wait runs the job inline and returns its code.
_want wf_wait && {
  _begin wf_wait
  ss_new_workspace "$WORK/wf_wait"
  cp "$HERE"/fixtures/wf_wait/* .; chmod +x ./*.sh 2>/dev/null || true
  ss_wf_activate
  bash ./submit.sh > wf.log 2>&1; rc=$?
  assert_rc "$rc" 0 "submit script exits 0"
  # --wait executed each job INLINE (no explicit ss_wf_drain was called)
  assert_file "ran.okp"                     "sbatch --wait ran okp inline"
  assert_file "ran.badp"                    "sbatch --wait ran badp inline"
  assert_contains "waitlog.txt" "okp rc=0"  "--wait returns 0 for a successful job"
  assert_contains "waitlog.txt" "badp rc=1" "--wait returns the failing job's exit code"
  echo
}

# ===========================================================================
# Migration Assistant — scancel + squeue live state.
_want wf_scancel && {
  _begin wf_scancel
  ss_new_workspace "$WORK/wf_scancel"
  cp "$HERE"/fixtures/wf_scancel/* .; chmod +x ./*.sh 2>/dev/null || true
  ss_wf_activate
  bash ./submit.sh > wf.log 2>&1; rc=$?
  assert_rc "$rc" 0 "submit script exits 0"
  # BEFORE drain: A was scancelled, so squeue shows exactly B and C pending
  squeue -h > q_before.txt 2>&1
  assert_rc 2 "$(grep -c 'PD' q_before.txt)" "squeue shows 2 pending jobs (A was scancelled)"
  ss_wf_drain
  assert_no_file "ran.A"  "scancelled job A never ran"
  assert_file    "ran.C"  "independent job C still ran"
  assert_no_file "ran.B"  "B (afterok:A) orphaned by A's cancellation"
  # AFTER drain: nothing active -> squeue empty
  squeue -h > q_after.txt 2>&1
  assert_empty "q_after.txt" "squeue is empty after drain (all jobs settled)"
  echo
}

# ===========================================================================
# Migration Assistant — resource-sim fault injection into a user workflow.
_want wf_inject && {
  _begin wf_inject
  # All failure modes propagate identically through the DAG (the original user
  # script isn't SmartSlurm-aware, so afterok/afternotok only care success/fail),
  # but the injected mode must be faithfully LABELLED so a later differential test
  # can drive the rAP port's mode-specific recovery (more-mem vs more-time).
  for mode in oom timeout nonzero; do
    ss_new_workspace "$WORK/wf_inject_$mode"
    cp "$HERE"/fixtures/wf_inject/* .; chmod +x ./*.sh 2>/dev/null || true
    ss_wf_activate
    bash ./submit.sh > wf.log 2>&1
    export MOCK_FAIL_JOB=align MOCK_FAIL_MODE="$mode"
    ss_wf_drain
    unset MOCK_FAIL_JOB MOCK_FAIL_MODE
    wfdir="$MOCK_SLURM_STATE/wf"
    assert_file    "ran.prep"                        "[$mode] upstream step ran"
    assert_no_file "ran.align"                       "[$mode] injected step did not complete"
    assert_no_file "ran.post"                        "[$mode] afterok downstream orphaned"
    assert_file    "ran.cleanup"                     "[$mode] afternotok cleanup fired"
    assert_contains "$wfdir/.order" "FAILED($mode)"  "[$mode] drain recorded the mode faithfully"
  done
  echo
}

# ===========================================================================
# Migration Assistant — tiered execution (real payload vs stub when a tool is
# unavailable). Same fixture run three ways.
_want wf_tiered && {
  _begin wf_tiered
  # real tier: the missing tool makes toolstep fail, orphaning its consumer
  ss_new_workspace "$WORK/wf_tiered_real"
  cp "$HERE"/fixtures/wf_tiered/* .; chmod +x ./*.sh 2>/dev/null || true
  ss_wf_activate
  bash ./submit.sh > wf.log 2>&1
  ss_wf_drain
  assert_no_file "ran.toolstep"  "real tier: step with a missing tool fails"
  assert_no_file "ran.consumer"  "real tier: afterok consumer orphaned by the failure"
  # stub tier (explicit): stub toolstep so the wiring still validates
  ss_new_workspace "$WORK/wf_tiered_stub"
  cp "$HERE"/fixtures/wf_tiered/* .; chmod +x ./*.sh 2>/dev/null || true
  ss_wf_activate
  bash ./submit.sh > wf.log 2>&1
  export SS_WF_STUB="toolstep"; ss_wf_drain; unset SS_WF_STUB
  wfdir="$MOCK_SLURM_STATE/wf"
  assert_no_file "ran.toolstep"  "stub tier: real payload skipped"
  assert_file    "ran.consumer"  "stub tier: consumer runs (toolstep stubbed OK)"
  assert_contains "$wfdir/.stubs" "toolstep" "explicit stub was logged"
  # auto stub-missing: a --wrap step calling a missing tool is auto-stubbed
  ss_new_workspace "$WORK/wf_tiered_auto"
  cp "$HERE"/fixtures/wf_tiered/* .; chmod +x ./*.sh 2>/dev/null || true
  ss_wf_activate
  bash ./submit.sh > wf.log 2>&1
  export SS_WF_STUB_MISSING=1; ss_wf_drain; unset SS_WF_STUB_MISSING
  wfdir="$MOCK_SLURM_STATE/wf"
  assert_contains "$wfdir/.stubs" "wrapstep" "SS_WF_STUB_MISSING auto-stubbed the missing-tool --wrap step"
  echo
}

# ===========================================================================
# Migration Assistant — realistic end-to-end pipeline (output correctness +
# fan-in failure orphaning), beyond the single-feature wf_* unit blocks.
_want wf_pipeline && {
  _begin wf_pipeline
  # happy path: verify the COMPUTED outputs, not just that jobs ran
  ss_new_workspace "$WORK/wf_pipeline"
  cp "$HERE"/fixtures/wf_pipeline/* .; chmod +x ./*.sh 2>/dev/null || true
  ss_wf_activate
  bash ./submit.sh > wf.log 2>&1; rc=$?
  assert_rc "$rc" 0 "submit script exits 0"
  ss_wf_drain
  assert_rc 2 "$(cat wc.s1.txt 2>/dev/null)"        "s1 word count correct (Hello World -> 2)"
  assert_rc 4 "$(cat wc.s2.txt 2>/dev/null)"        "s2 word count correct (4)"
  assert_rc 1 "$(cat wc.s3.txt 2>/dev/null)"        "s3 word count correct (1)"
  assert_rc 7 "$(cat grandtotal.txt 2>/dev/null)"   "multi-id afterok fan-in summed all branches (7)"
  assert_contains report.txt "total words: 7"        "report reflects the fan-in total"
  assert_contains report.txt "samples: 3"            "report counts all samples"
  # failure path: inject OOM into one branch -> fan-in merge + report orphaned
  ss_new_workspace "$WORK/wf_pipeline_fail"
  cp "$HERE"/fixtures/wf_pipeline/* .; chmod +x ./*.sh 2>/dev/null || true
  ss_wf_activate
  bash ./submit.sh > wf.log 2>&1
  export MOCK_FAIL_JOB=wc_s2 MOCK_FAIL_MODE=oom
  ss_wf_drain
  unset MOCK_FAIL_JOB MOCK_FAIL_MODE
  assert_file    "wc.s1.txt"      "unrelated branch still ran"
  assert_no_file "wc.s2.txt"      "injected branch did not complete"
  assert_no_file "allcounts.txt"  "multi-id afterok fan-in (merge) orphaned by one failed branch"
  assert_no_file "report.txt"     "report orphaned (afterok: merge)"
  echo
}

# ===========================================================================
# Migration Assistant — arm isolation (no crosstalk). Run the sbatch (wf-mode) arm
# and the rAP arm in ONE block but SEPARATE workspaces, with SS_WF_MODE deliberately
# left set by the wf arm, and prove the rAP arm still executes (ss_run_pipeline
# clears SS_WF_MODE, so runAsPipeline -> ssbatch -> mock sbatch is not hijacked into
# register-don't-run). This guards the Scenario-2 differential-test setup.
_want wf_isolation && {
  _begin wf_isolation
  # Arm 1: an ORIGINAL user workflow in wf-mode (register-don't-run), workspace A.
  ss_new_workspace "$WORK/wf_iso_orig"
  cp "$HERE"/fixtures/wf_afterok/* .; chmod +x ./*.sh 2>/dev/null || true
  ss_wf_activate                      # sets SS_WF_MODE=1 and leaves it set
  bash ./submit.sh > wf.log 2>&1
  ss_wf_drain
  assert_file "marker.txt"            "wf arm: original afterok workflow ran"
  assert_file "fromB.txt"             "wf arm: downstream step ran"
  # Arm 2: an rAP pipeline in a SEPARATE workspace, with SS_WF_MODE STILL set.
  ss_new_workspace "$WORK/wf_iso_rap"
  cp "$HERE"/fixtures/clean_2step/* .
  ss_run_pipeline "./pipeline.sh 777" "sbatch -p short -c 1 --mem 4G -t 50:0"; rc=$?
  assert_rc "$rc" 0 "rAP arm runs despite SS_WF_MODE being set by the wf arm"
  assert_file "smartSlurmLog/2.1.mergeNumber.success" "rAP pipeline actually EXECUTED (not register-don't-run)"
  assert_nonempty "all777.txt"        "rAP produced real output (no crosstalk hijack)"
  assert_no_file "$MOCK_SLURM_STATE/wf/10001.run" "rAP workspace carries no leaked wf-mode jobs"
  echo
}

# ===========================================================================
# Migration Assistant — the mock sbatch job-id counter is atomic under concurrency.
_want wf_jobseq && {
  _begin wf_jobseq
  ss_new_workspace "$WORK/wf_jobseq"
  ss_wf_activate
  # fire many CONCURRENT submissions against ONE mock state; the atomic counter
  # must hand every one a unique id (a plain read-modify-write would collide).
  N=40
  for i in $(seq 1 "$N"); do ( sbatch -J "j$i" --wrap "true" >/dev/null 2>&1 ) & done
  wait
  wfdir="$MOCK_SLURM_STATE/wf"
  got=$(ls "$wfdir"/*.run 2>/dev/null | wc -l | tr -d ' ')
  uniq=$(ls "$wfdir"/*.run 2>/dev/null | sed 's#.*/##; s/\.run$//' | sort -u | wc -l | tr -d ' ')
  assert_rc "$N" "$got"    "all $N concurrent submissions registered (no id collision overwrote a job)"
  assert_rc "$got" "$uniq" "every registered job id is unique"
  echo
}

# ===========================================================================
# Migration Assistant — silent data failure (exit 0 but no output). Unlike a hard
# failure (which orphans the afterok downstream), the consumer here RUNS, because
# afterok is satisfied by the exit code, then fails for real on the missing input.
_want wf_exit0_empty && {
  _begin wf_exit0_empty
  ss_new_workspace "$WORK/wf_exit0_empty"
  cp "$HERE"/fixtures/wf_exit0_empty/* .; chmod +x ./*.sh 2>/dev/null || true
  ss_wf_activate
  bash ./submit.sh > wf.log 2>&1; rc=$?
  assert_rc "$rc" 0 "submit script exits 0"
  ss_wf_drain
  wfdir="$MOCK_SLURM_STATE/wf"
  : > states.txt
  for f in "$wfdir"/*.name; do jid=$(basename "$f" .name)
    printf '%s %s\n' "$(cat "$f")" "$(cat "$wfdir/$jid.rc" 2>/dev/null)" >> states.txt; done
  assert_no_file  "out.txt"                "producer exited 0 but wrote no output (silent failure)"
  assert_contains states.txt "producer 0"  "producer recorded COMPLETED (exit 0)"
  crc=$(sed -n 's/^consumer //p' states.txt)
  if [ -n "$crc" ] && [ "$crc" != cancelled ]; then
    _pass "consumer RAN — afterok satisfied by exit-0, NOT orphaned (contrast: a hard failure cancels it)"
  else _fail "consumer ran (not orphaned)" "consumer rc='$crc'"; fi
  if [ -n "$crc" ] && [ "$crc" != cancelled ] && [ "$crc" != 0 ]; then
    _pass "consumer then FAILED on the missing input (rc=$crc)"
  else _fail "consumer failed on missing input" "consumer rc='$crc'"; fi
  assert_no_file  "result.txt"             "consumer produced no result (data failure propagated)"
  echo
}

# ===========================================================================
# Migration Assistant skill — DAG extraction from a hand-rolled SLURM workflow.
# Run the submit script register-don't-run, then extract_dag.sh must recover the
# nodes and the name-resolved dependency edges (per-sample chains + fan-in + tail).
_want ma_extract_dag && {
  _begin ma_extract_dag
  ss_new_workspace "$WORK/ma_extract_dag"
  cp "$HERE"/fixtures/wf_pipeline/* .; chmod +x ./*.sh 2>/dev/null || true
  ss_wf_activate
  bash ./submit.sh > wf.log 2>&1       # register every job; nothing runs
  tool="$HERE/../skills/migration-assistant/tools/extract_dag.sh"
  bash "$tool" > dag.tsv 2>dag.err
  assert_nonempty "dag.tsv"                            "extract_dag produced a DAG"
  assert_contains "dag.tsv" "NODE.*norm_s1"            "node norm_s1 present"
  assert_contains "dag.tsv" "NODE.*merge"              "node merge present"
  assert_contains "dag.tsv" "NODE.*report"             "node report present"
  assert_contains "dag.tsv" "wc_s1.*afterok.*norm_s1"  "edge: wc_s1 depends afterok on norm_s1 (ids resolved to names)"
  assert_contains "dag.tsv" "merge.*afterok.*wc_s1"    "edge: merge fan-in includes wc_s1"
  assert_contains "dag.tsv" "merge.*afterok.*wc_s3"    "edge: merge fan-in includes wc_s3"
  assert_contains "dag.tsv" "report.*afterok.*merge"   "edge: report depends afterok on merge"
  echo
}

# ===========================================================================
# Migration Assistant skill — render a logical plan to rAP, then round-trip it:
# the emitted pipeline must run under runAsPipeline and reproduce the original's
# outputs (fan-in total, report).
_want ma_emit_rap && {
  _begin ma_emit_rap
  ss_new_workspace "$WORK/ma_emit_rap"
  cp "$HERE"/fixtures/wf_pipeline/* .; chmod +x ./*.sh 2>/dev/null || true
  cp "$HERE"/fixtures/ma_emit_rap/plan.txt .
  tool="$HERE/../skills/migration-assistant/tools/emit_rap.sh"
  bash "$tool" plan.txt > pipeline.sh
  assert_valid_bash "pipeline.sh"                           "emitted rAP script is valid bash"
  assert_contains "pipeline.sh" "#@1,0,norm,,none,sbatch"   "step 1 norm, no dep"
  assert_contains "pipeline.sh" "#@2,1,wc,,none,sbatch"     "step 2 wc depends on step 1"
  assert_contains "pipeline.sh" "#@3,2,merge,,none,sbatch"  "step 3 merge (fan-in) depends on step 2"
  assert_contains "pipeline.sh" "#@4,3,report,,none,sbatch" "step 4 report depends on step 3"
  assert_contains "pipeline.sh" "for s in"                  "fan-out steps wrapped in a per-sample loop"
  # ROUND-TRIP: the emitted rAP pipeline runs and reproduces the original's outputs
  ss_run_pipeline "./pipeline.sh" "sbatch -p short -c 1 --mem 100M -t 5:00"; rc=$?
  assert_rc "$rc" 0                                         "emitted rAP pipeline runs under runAsPipeline"
  assert_rc 7 "$(cat grandtotal.txt 2>/dev/null)"          "rAP port reproduces the fan-in total (7)"
  assert_contains report.txt "total words: 7"              "rAP port reproduces the report"
  echo
}

# ===========================================================================
# Migration Assistant skill — parse_rap_dag.sh reads an rAP driver's #@ DAG back,
# resolving dependency step-numbers to names. Round-trips emit_rap: emit a driver
# from a plan, then recover the same step graph from the driver text.
_want ma_parse_dag && {
  _begin ma_parse_dag
  ss_new_workspace "$WORK/ma_parse_dag"
  cp "$HERE"/fixtures/ma_emit_rap/plan.txt .
  emit="$HERE/../skills/migration-assistant/tools/emit_rap.sh"
  parse="$HERE/../skills/migration-assistant/tools/parse_rap_dag.sh"
  bash "$emit" plan.txt > pipeline.sh
  bash "$parse" pipeline.sh | awk '{print $2, $NF}' > dag.txt     # "<name> dep=<names>"
  assert_nonempty "dag.txt"                     "parse_rap_dag produced a step DAG"
  assert_contains "dag.txt" "norm dep=-"        "root step (norm) has no dependency"
  assert_contains "dag.txt" "wc dep=norm"       "wc's dep step-number resolved to the name norm"
  assert_contains "dag.txt" "merge dep=wc"      "merge resolved to wc"
  assert_contains "dag.txt" "report dep=merge"  "report resolved to merge"
  # a dot-joined multi-dep driver resolves EVERY id to a name
  printf '#@1,0,a,,,sbatch x\n:\n#@end\n#@2,0,b,,,sbatch x\n:\n#@end\n#@3,1.2,c,,,sbatch x\n:\n#@end\n' > multi.sh
  bash "$parse" multi.sh | awk '{print $2, $NF}' > multi_dag.txt
  assert_contains "multi_dag.txt" "c dep=a.b"   "dot-joined multi-dep resolved to both names (a.b)"
  echo
}

# ===========================================================================
# Migration Assistant skill — the identity-contract comparator (unit).
_want ma_compare && {
  _begin ma_compare
  ss_new_workspace "$WORK/ma_compare"
  mkdir -p A B
  printf 'data\n7\n' > A/result.txt; cp A/result.txt B/result.txt
  printf 'run 2024-01-01 job 123\n' > A/stamped.txt       # same data, volatile fields
  printf 'run 2025-06-06 job 999\n' > B/stamped.txt
  printf '7\n' > A/diffme.txt; printf '8\n' > B/diffme.txt  # genuine difference
  printf 'x\n' > A/onlyA.txt                                # present in one run only
  cat > contract <<'C'
FILE result.txt
FILE stamped.txt
FILE diffme.txt
FILE onlyA.txt
NORMALIZE s/[0-9]{4}-[0-9]{2}-[0-9]{2}/<DATE>/g
NORMALIZE s/job [0-9]+/job <ID>/g
C
  tool="$HERE/../skills/migration-assistant/tools/compare_outputs.sh"
  bash "$tool" A B contract > cmp.out 2>&1; crc=$?
  assert_contains cmp.out "MATCH result"      "identical file matches"
  assert_contains cmp.out "MATCH stamped"     "volatile date/job normalized away -> match"
  assert_contains cmp.out "DIFF.*diffme"      "genuinely different file flagged"
  assert_contains cmp.out "DIFF.*onlyA"       "file present in only one run flagged"
  assert_contains cmp.out "RESULT: MISMATCH"  "overall mismatch reported"
  assert_rc 1 "$crc"                           "comparator exits nonzero on mismatch"
  printf 'FILE result.txt\n' > contract2       # all-match subset
  bash "$tool" A B contract2 > cmp2.out 2>&1; crc2=$?
  assert_contains cmp2.out "RESULT: MATCH"     "clean comparison reports match"
  assert_rc 0 "$crc2"                          "comparator exits 0 when all match"
  echo
}

# ===========================================================================
# Migration Assistant skill — full Scenario 2: ORIGINAL vs rAP PORT, differential.
_want ma_differential && {
  _begin ma_differential
  # Arm 1: the ORIGINAL hand-rolled workflow (wf-mode), its own workspace.
  ss_new_workspace "$WORK/ma_diff_orig"
  cp "$HERE"/fixtures/wf_pipeline/* .; chmod +x ./*.sh 2>/dev/null || true
  ss_wf_activate
  bash ./submit.sh > wf.log 2>&1
  ss_wf_drain
  origdir="$PWD"
  # Arm 2: the rAP PORT (emit from plan) in a SEPARATE workspace. ss_run_pipeline
  # clears SS_WF_MODE, so the two arms don't crosstalk.
  ss_new_workspace "$WORK/ma_diff_rap"
  cp "$HERE"/fixtures/wf_pipeline/* .; chmod +x ./*.sh 2>/dev/null || true
  cp "$HERE"/fixtures/ma_emit_rap/plan.txt .
  bash "$HERE/../skills/migration-assistant/tools/emit_rap.sh" plan.txt > pipeline.sh
  ss_run_pipeline "./pipeline.sh" "sbatch -p short -c 1 --mem 100M -t 5:00"; rc=$?
  assert_rc "$rc" 0 "rAP port runs under runAsPipeline"
  rapdir="$PWD"
  cat > contract <<'C'
FILE grandtotal.txt
FILE report.txt
FILE wc.*.txt
C
  cmp="$HERE/../skills/migration-assistant/tools/compare_outputs.sh"
  bash "$cmp" "$origdir" "$rapdir" contract > cmp.out 2>&1; crc=$?
  sed 's/^/    /' cmp.out
  assert_contains cmp.out "RESULT: MATCH"  "original and rAP port produce identical authoritative outputs"
  assert_rc 0 "$crc"                        "differential comparison passes"
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
  # Each failure mode must be attributed distinctly end-to-end: the mock srun emits
  # a mode-specific message/exit, the mock sacct a distinct state, and cleanUp maps
  # it to the right jobStatus (OOM / OOT / Fail) -- which is what drives SmartSlurm's
  # mode-specific requeue (more mem vs more time). "spec" = mode|srun-message|jobStatus.
  for spec in "oom|Out Of Memory|OOM" "timeout|DUE TO TIME LIMIT|OOT" "nonzero|Exited with exit code 1|Fail"; do
    mode=${spec%%|*}; rest=${spec#*|}; msg=${rest%|*}; status=${rest##*|}
    ss_new_workspace "$WORK/oom_injection_$mode"
    cp "$HERE"/fixtures/clean_2step/* .; chmod +x ./*.sh 2>/dev/null || true
    MOCK_FAIL_JOB=mergeNumber MOCK_FAIL_MODE="$mode" \
      ss_run_pipeline "./pipeline.sh 777" "sbatch -p short -c 1 --mem 4G -t 50:0"; rc=$?
    assert_rc "$rc" 0 "[$mode] runAsPipeline exits 0 (submission still succeeds)"
    assert_file    "smartSlurmLog/1.0.findNumber.1.success"      "[$mode] unrelated producer still succeeds"
    assert_no_file "smartSlurmLog/2.1.mergeNumber.success"       "[$mode] failed job did NOT touch .success"
    assert_contains "smartSlurmLog/2.1.mergeNumber.out" "$msg"                 "[$mode] srun failure message present in the job log"
    assert_contains "smartSlurmLog/2.1.mergeNumber.out" "jobStatus: $status"   "[$mode] cleanUp classifies the failure correctly via sacct"
  done
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
# Fate reconciliation: a job hard-killed (OOM/timeout) bypasses cleanUp's EXIT
# trap, leaving NO .success/.failed/report -- only sacct knows. Verify the grid
# flags it as unknown (not phantom-running), diagnose is read-only by default but
# resolves the real cause on an explicit reconcile, and jobFate writes a marker for
# a real failure but NOT for CANCELLED (orphan-safe), degrading gracefully on lag.
_want fate_reconcile && {
  _begin fate_reconcile
  ss_new_workspace "$WORK/fate_reconcile"
  d="$PWD/fr"; rm -rf "$d"; mkdir -p "$d/jobs"
  export MOCK_SLURM_STATE="$d"          # mock sacct reads $MOCK_SLURM_STATE/jobs/<id>.inject
  printf '%-10s %-12s %s\n' job_id depend_on job_flag > "$d/allJobs.txt"
  printf '%s %s %s\n' 801 0   1.0.setup.S1    >> "$d/allJobs.txt"
  printf '%s %s %s\n' 802 801 2.1.oomstep.S1  >> "$d/allJobs.txt"
  printf '%s %s %s\n' 803 801 3.1.cancstep.S1 >> "$d/allJobs.txt"
  : > "$d/1.0.setup.S1.success"
  for f in 2.1.oomstep.S1 3.1.cancstep.S1; do
    printf '========== SmartSlurm step output: %s ==========\nwork\n' "$f" > "$d/$f.out"
    touch -d '-5 min' "$d/$f.out" 2>/dev/null || touch "$d/$f.out"   # stale -> "vanished", not "fresh running"
  done
  printf '# fresh squeue snapshot; dead jobs 802/803 absent\n' > "$d/squeue_cache.txt"
  echo oom > "$d/jobs/802.inject"; echo cancelled > "$d/jobs/803.inject"

  ul="$(cd "$d" && checkRunGrid --list-unknown --log . 2>/dev/null)"
  printf '%s\n' "$ul" | grep -qx '2.1.oomstep.S1' \
    && _pass "grid flags vanished OOM step as unknown (not phantom-running)" \
    || _fail "grid flags vanished OOM step as unknown" "list='$ul'"

  nb=$(ls "$d"/*.failed 2>/dev/null | wc -l)
  dd="$(cd "$d" && diagnose.sh 2.1.oomstep.S1 --log . 2>/dev/null | head -1)"
  na=$(ls "$d"/*.failed 2>/dev/null | wc -l)
  { printf '%s' "$dd" | grep -q UNCONFIRMED && [ "$nb" = "$na" ]; } \
    && _pass "default diagnose is read-only (UNCONFIRMED, no sacct, no write)" \
    || _fail "default diagnose read-only" "card='$dd' before=$nb after=$na"

  ( cd "$d" && jobFate.sh --batch 2.1.oomstep.S1 3.1.cancstep.S1 --log . >/dev/null 2>&1 )
  [ -e "$d/2.1.oomstep.S1.failed" ] \
    && _pass "reconcile marks vanished OOM as .failed" \
    || _fail "reconcile marks vanished OOM as .failed"
  [ -e "$d/3.1.cancstep.S1.failed" ] \
    && _fail "reconcile must NOT mark CANCELLED .failed (orphan-safe)" \
    || _pass "reconcile leaves CANCELLED unmarked (orphan-safe)"

  rd="$(cd "$d" && diagnose.sh 2.1.oomstep.S1 --log . 2>/dev/null | head -1)"
  printf '%s' "$rd" | grep -q 'state=OOM' \
    && _pass "diagnose resolves vanished OOM after reconcile" \
    || _fail "diagnose resolves vanished OOM after reconcile" "card='$rd'"

  printf '%s %s %s\n' 804 801 4.1.lagstep.S1 >> "$d/allJobs.txt"
  printf '========== SmartSlurm step output: 4.1.lagstep.S1 ==========\nwork\n' > "$d/4.1.lagstep.S1.out"
  touch -d '-5 min' "$d/4.1.lagstep.S1.out" 2>/dev/null || true
  echo lag > "$d/jobs/804.inject"
  lg="$(cd "$d" && jobFate.sh 4.1.lagstep.S1 --log . 2>/dev/null)"
  { [ "$lg" = UNRESOLVED ] && [ ! -e "$d/4.1.lagstep.S1.failed" ]; } \
    && _pass "sacct lag -> UNRESOLVED, no marker (graceful)" \
    || _fail "sacct lag graceful" "got '$lg'"

  unset MOCK_SLURM_STATE
  echo
}

# ===========================================================================
# Hard-kill an arbitrary pipeline step mid-run. MOCK_FAIL_MODE=sigkill runs the
# targeted job detached in its own session and SIGKILLs the whole process group
# (script + EXIT trap + monitor) -- modeling a cgroup OOM / node kill where cleanUp
# NEVER runs (the case oom/timeout injection cannot produce, since there the script
# survives). sigterm sends SIGTERM instead, so the result depends on the job's trap:
# today's EXIT-only trap is bypassed (documents the gap a TERM-trap change closes).
# This is the substrate for faithful cleanUp signal/trap development in the sandbox.
_want hardkill && {
  _begin hardkill

  ss_new_workspace "$WORK/hardkill_kill"
  cp "$HERE"/fixtures/clean_2step/* .; chmod +x ./*.sh
  MOCK_FAIL_JOB=findNumber MOCK_FAIL_MODE=sigkill \
    ss_run_pipeline "./pipeline.sh 777" "sbatch -p short -c 1 --mem 4G -t 50:0" >/dev/null 2>&1

  fjid=$(awk 'NR>1 && $3 ~ /findNumber/ {print $1; exit}' smartSlurmLog/allJobs.txt)
  fflag=$(awk 'NR>1 && $3 ~ /findNumber/ {print $3; exit}' smartSlurmLog/allJobs.txt)
  { [ -n "$fflag" ] && [ ! -e "smartSlurmLog/$fflag.success" ] && [ ! -e "smartSlurmLog/$fflag.failed" ]; } \
    && _pass "SIGKILL'd step leaves NO .success/.failed (cleanUp's EXIT trap never ran)" \
    || _fail "SIGKILL'd step leaves no markers" "flag='$fflag'"
  [ -n "$fflag" ] && [ -e "smartSlurmLog/$fflag.out" ] \
    && _pass "SIGKILL'd step has a .out (it had started)" || _fail "SIGKILL'd step has .out"
  kst=$(sacct -n -j "$fjid" -o State 2>/dev/null | head -1 | tr -d ' ')
  [ "$kst" = OUT_OF_MEMORY ] \
    && _pass "sacct reports OUT_OF_MEMORY for the SIGKILL'd (vanished) step" \
    || _fail "sacct state for SIGKILL'd step" "got '$kst'"
  # the vanished step is reconcilable from sacct
  ( cd smartSlurmLog && jobFate.sh "$fflag" --log . >/dev/null 2>&1 )
  [ -n "$fflag" ] && [ -e "smartSlurmLog/$fflag.failed" ] \
    && _pass "jobFate reconciles the vanished SIGKILL to OOM (.failed written)" \
    || _fail "jobFate reconciles the vanished SIGKILL"

  ss_new_workspace "$WORK/hardkill_term"
  cp "$HERE"/fixtures/clean_2step/* .; chmod +x ./*.sh
  MOCK_FAIL_JOB=findNumber MOCK_FAIL_MODE=sigterm \
    ss_run_pipeline "./pipeline.sh 777" "sbatch -p short -c 1 --mem 4G -t 50:0" >/dev/null 2>&1
  tflag=$(awk 'NR>1 && $3 ~ /findNumber/ {print $3; exit}' smartSlurmLog/allJobs.txt)
  { [ -n "$tflag" ] && [ ! -e "smartSlurmLog/$tflag.success" ] && [ ! -e "smartSlurmLog/$tflag.failed" ]; } \
    && _pass "SIGTERM vs today's EXIT-only trap also leaves no markers (the trap-hardening gap)" \
    || _fail "SIGTERM vs EXIT-only trap" "flag='$tflag'"
  echo
}

# ===========================================================================
# OOMOOT: a job that hit BOTH the memory AND the time limit. sacct shows TIMEOUT on
# the alloc line and OUT_OF_MEMORY on a step, so cleanUp's whole-blob check sees both
# -> jobStatus=OOMOOT, and it escalates mem AND time together. Previously unreachable
# (the mock injected a single state); covers the dual-escalation branch a cleanUp
# edit (e.g. P6) would touch.
_want oomoot_escalate && {
  _begin oomoot_escalate
  ss_new_workspace "$WORK/oomoot_escalate"
  cp "$HERE"/fixtures/clean_2step/* .; chmod +x ./*.sh
  MOCK_FAIL_JOB=findNumber MOCK_FAIL_MODE=oomoot \
    ss_run_pipeline "./pipeline.sh 777" "sbatch -p short -c 1 --mem 4G -t 50:0" >/dev/null 2>&1
  oflag=$(awk 'NR>1 && $3 ~ /findNumber/ {print $3; exit}' smartSlurmLog/allJobs.txt)

  compgen -G "smartSlurmLog/$oflag.failed.requeued.*.oomoot" >/dev/null 2>&1 \
    && _pass "cleanUp classifies OOMOOT and takes the .oomoot requeue branch" \
    || _fail "cleanUp OOMOOT requeue marker" "flag='$oflag'"

  ost=$(diagnose.sh "$oflag" --log smartSlurmLog --json 2>/dev/null | grep -o '"state":"[^"]*"' | head -1 | cut -d'"' -f4)
  [ "$ost" = OOMOOT ] && _pass "diagnose reports state=OOMOOT" || _fail "diagnose OOMOOT state" "got '$ost'"

  # dual escalation: the .adjust carries both a raised mem (>4096) and a raised time (>50)
  if [ -e "smartSlurmLog/$oflag.adjust" ]; then
    read -r _isz amem amin _rest < "smartSlurmLog/$oflag.adjust"
    { [ "${amem:-0}" -gt 4096 ] 2>/dev/null && [ "${amin:-0}" -gt 50 ] 2>/dev/null; } \
      && _pass "OOMOOT escalates BOTH mem ($amem>4096) and time ($amin>50)" \
      || _fail "OOMOOT dual escalation" "mem=$amem min=$amin"
  else
    _fail "OOMOOT dual escalation" "no .adjust file"
  fi
  echo
}

# ===========================================================================
# P6 signal/trap test-first scaffolds (OPT-IN via _optin; never in the default
# suite). The "crystal-clear, timing-deterministic" harness the signal/trap work
# needs BEFORE any cleanUp/ssbatch edit. P6 is PARKED (2026-11-11): no bin/ is
# touched. p6_signal_sequence and p6_determinism exercise a self-contained MODEL
# of the proposed rAP job body so the timing SEMANTICS are an executable spec;
# cleanup_reentrancy (below) probes the REAL cleanUp and is expected red until it
# gains the guard. Determinism comes from BARRIERS, never wall-clock: the harness
# blocks until the job provably reaches a known point, then delivers the signal,
# then releases -- so the trapped signal always fires the same way. The model runs
# the tool as `& wait` (the A0 enabler: a trapped signal interrupts `wait` at once
# instead of being deferred, which is exactly why a plain FOREGROUND tool would
# silently defeat the trap). Run: ./test/run_tests.sh p6   (serial; opt-in blocks
# are not -jN eligible). See docs/planning/p6-signal-trap-plan.md (SmartSlurm).

# Emit the model job wrapper embodying the proposed P6 design to $1.
_p6_emit_wrapper() {
  cat > "$1" <<'WRAP'
#!/bin/bash
# MODEL of the proposed rAP job body (not the real ssbatch/cleanUp). The tool runs
# as `& wait` so a trapped signal interrupts the wait immediately; a run-once-
# guarded cleanup is trapped on EXIT/TERM/USR1. Env in: LOG LEADERPID TOOL_BARRIER
# TOOL_RELEASE SUCCESS TOOL_RC [CLEAN_BARRIER CLEAN_RELEASE].
echo $$ > "$LEADERPID"
_clean_done=
cleanup() {
  if [ -n "$_clean_done" ]; then echo "cleanup-reentry-blocked reason=${1:-?}" >> "$LOG"; return; fi
  _clean_done=1
  echo "cleanup-begin reason=${1:-exit}" >> "$LOG"
  if [ -n "${CLEAN_BARRIER:-}" ]; then : > "$CLEAN_BARRIER"; while [ ! -f "${CLEAN_RELEASE:-/nonexistent}" ]; do sleep 0.02; done; fi
  echo "cleanup-end reason=${1:-exit}" >> "$LOG"
}
trap 'cleanup exit' EXIT
trap 'cleanup TERM; exit 143' TERM
trap 'cleanup USR1; exit 138' USR1
echo start >> "$LOG"
( : > "$TOOL_BARRIER"; while [ ! -f "$TOOL_RELEASE" ]; do sleep 0.02; done; exit "${TOOL_RC:-0}" ) &
_tool=$!
wait "$_tool"; _rc=$?
echo "tool-done rc=$_rc" >> "$LOG"
if [ "$_rc" -eq 0 ]; then echo success >> "$LOG"; : > "$SUCCESS"; fi
WRAP
}

# Run the model wrapper once under deterministic barrier control. Uses the global
# $WRAPPER (set by each block). Args: $1 signal at the tool barrier (TERM|USR1|NONE),
# $2 tool exit code (default 0), $3 "post" to let the tool finish THEN signal (no-op
# case). Prints the normalized event log (newlines -> '|').
_p6_run() {
  local sig="$1" trc="${2:-0}" when="${3:-at-barrier}"
  local d="$PWD/p6rundir"; rm -rf "$d"; mkdir -p "$d"
  local log="$d/log" lp="$d/leaderpid" tb="$d/tool.bar" tr="$d/tool.rel" sf="$d/success.flag"
  LOG="$log" LEADERPID="$lp" TOOL_BARRIER="$tb" TOOL_RELEASE="$tr" SUCCESS="$sf" TOOL_RC="$trc" \
    setsid bash "$WRAPPER" & local sp=$!
  # BARRIER: block until the tool provably reached its barrier (traps are installed
  # before the tool child starts, so this also guarantees the trap is armed). Bounded
  # so a broken wrapper fails the asserts instead of hanging the suite.
  local t=0; while [ ! -f "$tb" ]; do sleep 0.005; t=$((t+1)); [ "$t" -gt 1000 ] && break; done
  local leader; leader=$(cat "$lp" 2>/dev/null || echo "")
  if [ "$when" = post ]; then
    : > "$tr"; wait "$sp" 2>/dev/null                              # finish + EXIT-cleanup first
    [ "$sig" != NONE ] && [ -n "$leader" ] && kill -"$sig" "$leader" 2>/dev/null  # now a no-op
  else
    [ "$sig" != NONE ] && [ -n "$leader" ] && kill -"$sig" "$leader" 2>/dev/null  # deliver mid-tool
    : > "$tr"; wait "$sp" 2>/dev/null
  fi
  # reap any orphaned tool child (its own setsid group id == leader pid); never the
  # harness's own group.
  [ -n "$leader" ] && [ "$leader" != "$$" ] && kill -KILL -"$leader" 2>/dev/null
  [ -f "$sf" ] && echo "SUCCESSFLAG" >> "$log"
  tr '\n' '|' < "$log" 2>/dev/null
}

_begins() { printf '%s' "$1" | grep -o 'cleanup-begin' | wc -l | tr -d ' '; }

# --- p6_signal_sequence: the deterministic sequence matrix (executable spec) ---
_optin p6_signal_sequence && {
  _begin p6_signal_sequence
  ss_new_workspace "$WORK/p6_signal_sequence"
  WRAPPER="$PWD/wrapper.sh"; _p6_emit_wrapper "$WRAPPER"

  # S0 back-compat: clean success, no signal -> tool runs, .success set, cleanup once (EXIT).
  o=$(_p6_run NONE 0)
  if [[ "$o" == *"start|"* && "$o" == *"tool-done rc=0|"* && "$o" == *"success|"* \
        && "$o" == *"cleanup-begin reason=exit|"* && "$o" == *"SUCCESSFLAG"* && "$(_begins "$o")" == 1 ]]; then
    _pass "S0 back-compat: clean success -> tool runs, .success set, cleanup runs once (EXIT)"
  else _fail "S0 back-compat clean-success path" "$o"; fi

  # S0b back-compat: tool exits nonzero, no signal -> no .success, cleanup once (EXIT).
  o=$(_p6_run NONE 1)
  if [[ "$o" == *"tool-done rc=1|"* && "$o" == *"cleanup-begin reason=exit|"* \
        && "$o" != *"SUCCESSFLAG"* && "$(_begins "$o")" == 1 ]]; then
    _pass "S0b back-compat: tool nonzero -> no .success, cleanup runs once (EXIT)"
  else _fail "S0b back-compat nonzero-exit path" "$o"; fi

  # S1 TERM mid-tool -> trap fires, cleanup runs ONCE as TERM, EXIT re-entry guarded, no success.
  o=$(_p6_run TERM 0)
  if [[ "$o" == *"cleanup-begin reason=TERM|"* && "$o" == *"cleanup-reentry-blocked reason=exit"* \
        && "$o" != *"SUCCESSFLAG"* && "$(_begins "$o")" == 1 ]]; then
    _pass "S1 TERM mid-tool: trap fires, cleanup body runs exactly once, EXIT re-entry guarded"
  else _fail "S1 TERM mid-tool must run cleanup exactly once" "begins=$(_begins "$o") $o"; fi

  # S2 USR1 mid-tool (the grace signal) -> same shape, reason=USR1.
  o=$(_p6_run USR1 0)
  if [[ "$o" == *"cleanup-begin reason=USR1|"* && "$o" == *"cleanup-reentry-blocked reason=exit"* \
        && "$(_begins "$o")" == 1 ]]; then
    _pass "S2 USR1 mid-tool: grace-signal trap fires, cleanup body runs exactly once"
  else _fail "S2 USR1 mid-tool must run cleanup exactly once" "begins=$(_begins "$o") $o"; fi

  # S6 post-completion signal -> no-op; cleanup already ran once via EXIT, .success intact.
  o=$(_p6_run TERM 0 post)
  if [[ "$o" == *"cleanup-begin reason=exit|"* && "$o" == *"SUCCESSFLAG"* && "$(_begins "$o")" == 1 ]]; then
    _pass "S6 post-completion signal: no-op, cleanup ran exactly once (EXIT), .success intact"
  else _fail "S6 post-completion signal must be a no-op" "begins=$(_begins "$o") $o"; fi
  echo
}

# --- p6_determinism: the anti-flake guarantee (zero variance over many repeats) ---
_optin p6_determinism && {
  _begin p6_determinism
  ss_new_workspace "$WORK/p6_determinism"
  WRAPPER="$PWD/wrapper.sh"; _p6_emit_wrapper "$WRAPPER"
  N=20
  for sig in TERM USR1 NONE; do
    declare -A seen=()
    for i in $(seq 1 "$N"); do o=$(_p6_run "$sig" 0); seen["$o"]=$(( ${seen["$o"]:-0} + 1 )); done
    nv=${#seen[@]}
    if [ "$nv" -eq 1 ]; then
      _pass "determinism: $sig at barrier -> IDENTICAL outcome across $N runs (zero variance)"
    else
      _fail "determinism: $sig FLAKED -- $nv distinct outcomes over $N runs" \
            "$(for k in "${!seen[@]}"; do echo "   ${seen[$k]}x : $k"; done)"
    fi
    unset seen
  done
  echo
}

# ===========================================================================
# cleanUp re-entrancy SCAFFOLD (P6 test-first). cleanUp runs from the job's trap;
# today only EXIT, so it runs exactly once. A trap-hardening change (adding TERM/
# USR1) creates a SECOND path to cleanUp, so it must become idempotent or it would
# double-append jobRecords and DOUBLE-REQUEUE (resubmit the job twice). This block
# invokes cleanUp twice on one failed job and asserts single effects.
#
# EXPECTED TO FAIL until cleanUp gains a re-entrancy guard -- it is both absent from
# ALL_BLOCKS and gated by _optin, so NO default run (serial `./run_tests.sh`, bootstrap,
# or -jN) ever touches it; it runs only when named. Run it explicitly during P6 dev:
#     ./test/run_tests.sh cleanup_reentrancy      # or: ./test/run_tests.sh p6
_optin cleanup_reentrancy && {
  _begin cleanup_reentrancy
  ss_new_workspace "$WORK/cleanup_reentrancy"
  cp "$HERE"/fixtures/clean_2step/* .; chmod +x ./*.sh
  # A CLEAN run (no injection): each step completes and cleanUp runs once. We probe
  # efficiency.tsv, which cleanUp appends UNCONDITIONALLY (one row per run, no dup
  # check, never pruned -- unlike the jobRecord, whose OOM rows the too-many-fails
  # prune removes). A second cleanUp on the same job must not add a second row.
  ss_run_pipeline "./pipeline.sh 777" "sbatch -p short -c 1 --mem 4G -t 50:0" >/dev/null 2>&1
  rflag=$(awk 'NR>1 && $3 ~ /findNumber/ {print $3; exit}' smartSlurmLog/allJobs.txt)
  rjid=$(awk 'NR>1 && $3 ~ /findNumber/ {print $1; exit}' smartSlurmLog/allJobs.txt)
  eff="smartSlurmLog/efficiency.tsv"
  effcnt() { awk -F'\t' -v f="$rflag" 'NR>1 && $2==f {n++} END{print n+0}' "$eff" 2>/dev/null; }
  rows1=$(effcnt)

  # invoke cleanUp a SECOND time with the exact trap-line args the job carries
  trapline=$(grep -o 'cleanUp\.sh [^"]*' "smartSlurmLog/$rflag.sh" 2>/dev/null | head -1)
  ( cd "$PWD" && export SLURM_JOBID="$rjid" SLURM_JOB_ID="$rjid"
    eval "$trapline" >/dev/null 2>&1 || true )
  rows2=$(effcnt)

  if [ "${rows1:-0}" -lt 1 ]; then
    _fail "cleanup_reentrancy setup: no efficiency.tsv row for the step" "rows1=$rows1 flag=$rflag"
  elif [ "$rows2" = "$rows1" ]; then
    _pass "cleanUp re-run is idempotent: no extra efficiency.tsv row ($rows1 -> $rows2)"
  else
    _fail "cleanUp re-run double-writes efficiency.tsv -- needs a re-entrancy guard (P6)" "$rows1 -> $rows2"
  fi
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
  n=$(awk -F, '$12=="COMPLETED" && $2=="scaleLin"{c++} END{print c+0}' "$jr")
  assert_rc 6 "$n" "six COMPLETED records, one per size"
  # size 4000 recorded exactly per model: mem=0.5*4000+200=2200, time=0.02*4000+10=90
  # v2 columns: ...,step,ref,sample,jobID,inputSize,memUsed,memAllocated,timeUsed,...,status,...
  assert_contains "$jr" ",scaleLin,none,[^,]*,[^,]*,4000,2200,[^,]*,90,[^,]*,[^,]*,COMPLETED," "size 4000 row matches model exactly"
  ss_fit scaleLin none >/dev/null
  mstat="$smartSlurmJobRecordDir/stats/unknown.scaleLin.none.mem.stat"
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
  exact1000=$(awk -F, '$6==1000 && $2=="scaleLin"{print $7; exit}' "$jr")
  if [ "$exact1000" != "700" ]; then _pass "noise perturbs recorded mem (size 1000: $exact1000 != exact 700)"
  else _fail "noise perturbs recorded mem" "size-1000 mem is exactly 700 (no noise?)"; fi
  ss_fit scaleLin none >/dev/null
  mstat="$smartSlurmJobRecordDir/stats/unknown.scaleLin.none.mem.stat"
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
  n=$(awk -F, '$12=="COMPLETED" && $2=="scaleExp"{c++} END{print c+0}' "$jr")
  assert_rc 6 "$n" "six COMPLETED exp records"
  # recorded curve is convex (exp), not linear: the step 6000->8000 jump far
  # exceeds the step 1000->3000 jump (a straight line would make them equal-ish).
  read -r m1 m3 m6 m8 <<<"$(awk -F, '$2=="scaleExp"{a[$6]=$7} END{print a[1000],a[3000],a[6000],a[8000]}' "$jr")"
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
  clean=$(stat_val Finala "$smartSlurmJobRecordDir/stats/unknown.cjob.none.mem.stat")
  assert_near 0.5 "$clean" 0.05 "clean history fits slope ~0.5"
  # a large-input job that exits 0 having done nothing (or whose inputSize was
  # mis-measured) records a near-zero mem at a big size -> drags the slope down
  ss_seed_record cjob none 10000 60 30 COMPLETED 1
  ss_fit cjob none >/dev/null
  dirty=$(stat_val Finala "$smartSlurmJobRecordDir/stats/unknown.cjob.none.mem.stat")
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
  before=$(stat_val Finala "$smartSlurmJobRecordDir/stats/unknown.xjob.none.mem.stat")
  # add OOM rows with absurd mem at a big size: if they leaked into the fit the
  # slope would explode. They must be ignored (estimateResource greps COMPLETED).
  ss_seed_record xjob none 10000 999999 120 OOM 4
  ss_fit xjob none >/dev/null
  after=$(stat_val Finala "$smartSlurmJobRecordDir/stats/unknown.xjob.none.mem.stat")
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
  assert_rc 4 "$(awk -F, '$2=="progB"{c++}END{print c+0}' "$jr")" "progB starts with 4 unrelated records"
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
    # --pipeline unknown: the seeded history (ss_seed_record) lands in pipeline
    # 'unknown', so the triggering job must run in the SAME pipeline for its
    # fail-prune to apply to that history. (cleanUp now scopes the prune by
    # pipeline; a script-derived pipeline name would correctly NOT touch the
    # 'unknown' rows.)
    MOCK_FAIL_MODE=oom timeout 115 bash "$SMARTSLURM_ROOT/bin/runAsPipeline" \
      --script "./pipeline.sh in.dat" --pipeline unknown --sbatch-options "sbatch -p short -c 1 --mem 4G -t 60:0" \
      --tmp noTmp >rap.log 2>&1 )
  rm -f "$HERE/models/_prune.model"
  qa() { awk -F, -v p="$1" -v s="$2" '$2==p && $6==s && $12=="COMPLETED"{c++}END{print c+0}' "$jr"; }
  assert_file "$(ls "$smartSlurmJobRecordDir"/jobRecord.back.*.txt 2>/dev/null | head -1)" "fail-pruning backed up jobRecord.txt"
  # Conservative suspect-removal: keep-filter drops the offending program's failed
  # rows AND any COMPLETED row disproven by a failure (size>=failing size and
  # mem<=proven-insufficient allocation). Unrelated programs and validated
  # (consistent) COMPLETED history are kept.
  assert_rc 4 "$(awk -F, '$2=="progB"{c++}END{print c+0}' "$jr")"  "unrelated progB preserved (4)"
  assert_rc 3 "$(qa progA 1000)"   "progA legit small history preserved (size 1000 x3)"
  assert_rc 1 "$(qa progA 15000)"  "progA legit HIGH-mem large record preserved (size 15000)"
  assert_rc 0 "$(qa progA 12000)"  "progA BOGUS under-reporter removed (size 12000, ~0 mem)"
  echo
  fi
}

# ===========================================================================
# legacy_pooling -- standalone ssbatch (pipeline=unknown) must use migrated jR1
# (pipeline=legacy) history for estimation. v1 had no pipeline column, so a v1
# standalone run and a v1 pipeline step are indistinguishable -- both migrate to
# legacy; we default to inclusion so that history still sizes standalone jobs. The
# {unknown, legacy} bucket must NOT reach into a NAMED pipeline's records, and a
# fresh named pipeline must still bootstrap from legacy via its thin-group fallback.
_want legacy_pooling && {
  _begin legacy_pooling
  if ! command -v gnuplot >/dev/null 2>&1; then
    _pass "(skip) gnuplot unavailable -- size-based fit cannot build"
  else
  ss_new_workspace "$WORK/legacy_pooling" >/dev/null
  mkdir -p smartSlurmLog
  mid=3000
  # mem field (col 2) of an estimate for (prog [, pipeline]) at size $mid; clear any
  # cached stat first so each call rebuilds from the current store. No pipeline arg
  # => estimateResource defaults to pipeline=unknown (a standalone ssbatch run).
  est_mem() {
    local prog="$1" pipe="${2:-}"
    rm -f "$smartSlurmJobRecordDir/stats/"*".$prog.none.mem.stat" 2>/dev/null
    if [ -n "$pipe" ]; then
      estimateResource.sh "$prog" none ",jobSize:$mid" lp 4096 120 "" "$pipe" 2>/dev/null | cut -d' ' -f2
    else
      estimateResource.sh "$prog" none ",jobSize:$mid" lp 4096 120 2>/dev/null | cut -d' ' -f2
    fi
  }
  def=$(est_mem nohist)   # bare default: a program with NO history anywhere

  # --- 1) legacy history bootstraps a standalone (unknown) estimate -----------
  # 5 legacy COMPLETED rows (pipeline=legacy), clean linear mem = 2*size.
  for s in 1000 2000 3000 4000 5000; do ss_seed_record leg none $s $((2*s)) 30 COMPLETED 1 1 legacy; done
  unk=$(est_mem leg)          # unknown query (standalone)
  legq=$(est_mem leg legacy)  # query legacy directly (same rows)
  echo "    (info) default(no history)=$def ; standalone(unknown)=$unk ; legacy-query=$legq"
  if [ -n "$unk" ] && [ "$unk" != "$def" ]; then
    _pass "standalone (unknown) uses legacy history (est $unk != bare default $def)"
  else _fail "standalone uses legacy history" "unknown=$unk default=$def"; fi
  assert_rc "$legq" "$unk" "standalone pools exactly the legacy bucket (unknown == legacy-query)"

  # --- 2) isolation: unknown must NOT pool a NAMED pipeline's history ----------
  # a different program whose ONLY history is under a named pipeline (no legacy).
  for s in 1000 2000 3000 4000 5000; do ss_seed_record iso none $s $((3*s)) 30 COMPLETED 1 1 somepipe; done
  unkiso=$(est_mem iso)       # unknown query; must not see somepipe's rows
  assert_rc "$def" "$unkiso" "standalone (unknown) does NOT pool a named pipeline (est == bare default)"

  # --- 3) a fresh NAMED pipeline still bootstraps from legacy (unchanged) ------
  named=$(est_mem leg freshpipe)   # freshpipe has no own rows -> thin fallback pools legacy
  if [ -n "$named" ] && [ "$named" != "$def" ]; then
    _pass "a fresh named pipeline still bootstraps from legacy (est $named != default)"
  else _fail "named pipeline bootstraps from legacy" "named=$named default=$def"; fi
  assert_rc "$legq" "$named" "named-pipeline thin fallback pools the same legacy history"
  unset -f est_mem
  echo
  fi
}

# ===========================================================================
# cross_pipeline -- the jobRecords grouping key is (pipeline, step, ref). Guards
# the two data-loss regressions (records from different pipelines wrongly deduped;
# a fail-prune in one pipeline deleting another's COMPLETED history) and the
# "standalone unknown never pools into / broadens a pipeline fit" contract. The
# pre-existing suite ran everything as pipeline=unknown, so it never exercised
# cross-pipeline isolation -- which is exactly why these bugs shipped.
# ===========================================================================
_want cross_pipeline && {
  _begin cross_pipeline
  if ! command -v gnuplot >/dev/null 2>&1; then
    echo "    (skipped: gnuplot unavailable)"; echo
  else
  ss_new_workspace "$WORK/cross_pipeline" >/dev/null
  jr="$smartSlurmJobRecordDir/jobRecord.txt"

  # ---- A: dedup must not drop a distinct pipeline's row. A const model makes both
  # runs measure the SAME mem at the SAME input, so pipeB collides with pipeA's
  # (inputSize,memUsed); the old step+ref-only dedup dropped it.
  cat > "$HERE/models/_xp.model" <<EOF
MODEL_PROGRAM=prog
MODEL_REF=none
MEM_FORM=const
MEM_B=1500
TIME_FORM=const
TIME_B=20
SEED=1
EOF
  xp_run() { # $1 = pipeline name (shares the block's smartSlurmJobRecordDir)
    local pd="$WORK/cross_pipeline/$1"; mkdir -p "$pd/.mockstate/jobs"
    ( cd "$pd"; export MOCK_SLURM_STATE="$pd/.mockstate"; export SS_MODEL="$HERE/models/_xp.model"
      truncate -s 5000K in.dat
      printf '#!/bin/sh\ninput=$1\n#@1,0,prog,,input,sbatch -p short -c 1 --mem 4G -t 60:0\ncat $input >/dev/null\n\n' > pipeline.sh
      timeout 115 bash "$SMARTSLURM_ROOT/bin/runAsPipeline" \
        --script "./pipeline.sh in.dat" --pipeline "$1" \
        --sbatch-options "sbatch -p short -c 1 --mem 4G -t 60:0" --tmp noTmp >rap.log 2>&1 )
  }
  xp_run pipeA
  xp_run pipeB
  rm -f "$HERE/models/_xp.model"
  a=$(awk -F, '$1=="pipeA"&&$2=="prog"&&$12=="COMPLETED"{c++}END{print c+0}' "$jr")
  b=$(awk -F, '$1=="pipeB"&&$2=="prog"&&$12=="COMPLETED"{c++}END{print c+0}' "$jr")
  { [ "${a:-0}" -ge 1 ] && [ "${b:-0}" -ge 1 ]; } \
    && _pass "cross-pipeline dedup: pipeB row kept despite colliding (inputSize,memUsed) with pipeA (pipeA=$a pipeB=$b)" \
    || _fail "cross-pipeline dedup drops a distinct pipeline's row" "pipeA=$a pipeB=$b (both should be >=1)"

  # ---- B: a fail-prune triggered in pipeA must NOT delete pipeB's COMPLETED
  # history, and must STILL prune pipeA's own contaminated rows.
  ss_new_workspace "$WORK/cross_pipeline_prune" >/dev/null
  jr="$smartSlurmJobRecordDir/jobRecord.txt"
  ss_seed_record prog none 1000  700  15  COMPLETED 3 700  pipeA   # pipeA legit small
  ss_seed_record prog none 12000 50   30  COMPLETED 1 50   pipeA   # pipeA BOGUS under-reporter (prune target)
  ss_seed_record prog none 9000  8000 120 OOM       6 8000 pipeA   # pipeA proven-insufficient allocs
  ss_seed_record prog none 12000 50   30  COMPLETED 1 50   pipeB   # pipeB row at the colliding key -- MUST survive
  ss_seed_record prog none 1000  900  20  COMPLETED 4 900  pipeB   # pipeB other history -- MUST survive
  cat > "$HERE/models/_xpp.model" <<EOF
MODEL_PROGRAM=prog
MODEL_REF=none
MEM_FORM=const
MEM_B=80000
TIME_FORM=const
TIME_B=30
SEED=1
EOF
  pd="$WORK/cross_pipeline_prune/oomjob"; mkdir -p "$pd/.mockstate/jobs"
  ( cd "$pd"; export MOCK_SLURM_STATE="$pd/.mockstate"; export SS_MODEL="$HERE/models/_xpp.model"
    truncate -s 9000K in.dat
    printf '#!/bin/sh\ninput=$1\n#@1,0,prog,,input,sbatch -p short -c 1 --mem 4G -t 60:0\ncat $input >/dev/null\n\n' > pipeline.sh
    MOCK_FAIL_MODE=oom timeout 115 bash "$SMARTSLURM_ROOT/bin/runAsPipeline" \
      --script "./pipeline.sh in.dat" --pipeline pipeA \
      --sbatch-options "sbatch -p short -c 1 --mem 4G -t 60:0" --tmp noTmp >rap.log 2>&1 )
  rm -f "$HERE/models/_xpp.model"
  bq() { awk -F, -v P="$1" -v s="$2" '$1==P&&$2=="prog"&&$6==s&&$12=="COMPLETED"{c++}END{print c+0}' "$jr"; }
  assert_rc 1 "$(bq pipeB 12000)" "cross-pipeline prune: pipeB COMPLETED at the colliding key SURVIVES pipeA's fail-prune"
  assert_rc 4 "$(bq pipeB 1000)"  "cross-pipeline prune: pipeB's other history preserved (4)"
  assert_rc 0 "$(bq pipeA 12000)" "within-pipeline prune still works: pipeA's BOGUS under-reporter removed"

  # ---- C: standalone unknown never pools / broadens. With <3 of its OWN rows it
  # must not borrow from named pipelines (it regresses separately) -> no fit.
  ss_new_workspace "$WORK/cross_pipeline_unknown" >/dev/null
  jr="$smartSlurmJobRecordDir/jobRecord.txt"
  ss_seed_record prog none 1000 1000 10 COMPLETED 1 1000 pipeA
  ss_seed_record prog none 2000 2000 20 COMPLETED 1 2000 pipeA
  ss_seed_record prog none 3000 3000 30 COMPLETED 1 3000 pipeA
  ss_seed_record prog none 4000 4000 40 COMPLETED 1 4000 pipeA
  ss_seed_record prog none 5000 5000 50 COMPLETED 1 5000 pipeA   # 5 NAMED rows (a clean fit)
  ss_seed_record prog none 1500 1500 15 COMPLETED 2 1500 unknown # only 2 unknown rows (<3)
  rm -f "$smartSlurmJobRecordDir/stats/unknown.prog.none.mem.stat"
  estimateResource.sh prog none ",jobSize:3000" fp 4096 120 >/dev/null 2>&1   # fits pipeline=unknown
  if [ -s "$smartSlurmJobRecordDir/stats/unknown.prog.none.mem.stat" ]; then
    _fail "standalone unknown must not pool named pipelines" "unknown produced a fit from <3 own rows (pooled pipeA)"
  else
    _pass "standalone unknown regresses separately: <3 own rows => no fit (named pipelines NOT pooled)"
  fi
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
  rf=$(awk -F, -v a=ff -v b=none -v t="$(date -d '1 hour ago' +%s)" '$12!="COMPLETED"&&$15>t&&$2==a&&$3==b' "$jr" | wc -l)
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
      --tmp noTmp --pipeline "$SS_PIPELINE" >rap.log 2>&1 )
  rm -f "$HERE/models/_mrf.model"
  mkdir -p smartSlurmLog
  oomAllocs=$(awk -F, '$2=="mrf" && $12 ~ /OOM/{print $8}' "$jr" | tr '\n' ' ')
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
    # calculateMemTime is invoked here without --pipeline -> statKey = unknown.prog.ref
    printf 'Finala=%s\nFinalb=%s\nMaximum=%s\nSTDFIT=%s\nSCount=%s\n' "$3" "$4" "$5" "$6" "$7" > "$sdir/unknown.$1.$2.mem.stat"
    printf 'Finala=%s\nFinalb=%s\nMaximum=%s\nSTDFIT=%s\nSCount=%s\nRSquare=%s\n' "$3" "$4" "$5" "$6" "$7" 0.5 > "$sdir/unknown.$1.$2.time.stat"
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
  printf '%s\n' "1000 1500" "2000 1600" "4000 1550" "8000 1700" > "$sdir/unknown.P.none.mem"
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
# `smartslurm config` — the overlay config model + CLI. User edits go to a small
# ~/.smartSlurm/config/user.conf (never git-tracked), sourced LAST by the repo
# config.txt, so changes win, never clobber the repo file (Q1), and new repo
# defaults still load because the repo file is sourced too, not replaced (Q2).
_want config_overlay && {
  _begin config_overlay
  h="$WORK/config_overlay_home"; rm -rf "$h"; mkdir -p "$h"; o="$h/out.txt"
  ss_scn "$h" > "$o" 2>&1 <<'SCN'
. "$SMARTSLURM_HOME/bin/smartslurm.sh"
smartslurm config set defaultMem 8192 >/dev/null
echo "GET=$(smartslurm config get defaultMem)"
echo "OTHER=$(smartslurm config get defaultTime)"          # a repo default NOT in user.conf
echo "UCKEYS=$(sed -n 's/^export \([A-Za-z_0-9]*\)=.*/\1/p' "$HOME/.smartSlurm/config/user.conf" | tr '\n' ',')"
grep -q '^export defaultMem=4096' "$SMARTSLURM_HOME/config/config.txt" && echo "REPO_UNCHANGED=yes"
if smartslurm config set defaultMem notAnInt >/dev/null 2>&1; then echo "BADSET=accepted"; else echo "BADSET=rejected"; fi
echo "STILL=$(smartslurm config get defaultMem)"
smartslurm config set myCustomKnob hi 2>&1 | grep -q 'not a known key' && echo "CUSTOM=flagged"
smartslurm config unset defaultMem >/dev/null
echo "AFTER_UNSET=$(smartslurm config get defaultMem)"
SCN
  assert_contains "$o" "GET=8192"           "set then get returns the new value"
  assert_contains "$o" "OTHER=120"          "a repo default NOT overridden still loads (overlay, not replacement = Q2)"
  assert_contains "$o" "UCKEYS=defaultMem," "user.conf holds ONLY the changed key"
  assert_contains "$o" "REPO_UNCHANGED=yes" "repo config.txt is untouched — user edits never clobber it (Q1)"
  assert_contains "$o" "BADSET=rejected"    "a bad int value is rejected (schema validation)"
  assert_contains "$o" "STILL=8192"         "the rejected set did not change the stored value"
  assert_contains "$o" "CUSTOM=flagged"     "an unknown key is flagged custom (unvalidated)"
  assert_contains "$o" "AFTER_UNSET=4096"   "unset reverts to the repo default"
  echo
}

# ===========================================================================
# `smartslurm config` guided editors drive real behavior: refgroup -> jrRefClass,
# partition -> the partition<N>TimeLimit the scheduler reads.
_want config_guided && {
  _begin config_guided
  h="$WORK/config_guided_home"; rm -rf "$h"; mkdir -p "$h"; o="$h/out.txt"
  ss_scn "$h" > "$o" 2>&1 <<'SCN'
. "$SMARTSLURM_HOME/bin/smartslurm.sh"
smartslurm config refgroup add dm3 dm6 --group fly >/dev/null
jrc() { . "$SMARTSLURM_HOME/config/config.txt" >/dev/null 2>&1; . "$SMARTSLURM_HOME/bin/jobRecordSchema.sh" >/dev/null 2>&1; jrRefClass "$1"; }
echo "JRC_DM6=$(jrc dm6)"
echo "JRC_HG38=$(jrc hg38)"
smartslurm config partition short --time-limit 1440 >/dev/null
echo "PART1TL=$(smartslurm config get partition1TimeLimit)"
SCN
  assert_contains "$o" "JRC_DM6=fly"     "refgroup add makes jrRefClass pool dm6 under 'fly'"
  assert_contains "$o" "JRC_HG38=mammal" "an untouched ref keeps its default group (mammal)"
  assert_contains "$o" "PART1TL=1440"    "partition editor sets partition1TimeLimit (resolved from name 'short')"
  echo
}

# ===========================================================================
# `smartslurm config migrate` converts a legacy FULL ~/.smartSlurm/config/config.txt
# into the overlay: keep only the keys that differ from the repo defaults, back up
# the legacy file.
_want config_migrate && {
  _begin config_migrate
  h="$WORK/config_migrate_home"; rm -rf "$h"; mkdir -p "$h"; o="$h/out.txt"
  ss_scn "$h" > "$o" 2>&1 <<'SCN'
. "$SMARTSLURM_HOME/bin/smartslurm.sh"
mkdir -p "$HOME/.smartSlurm/config"
# a legacy full copy: identical to the repo except one changed scalar + one custom key
sed 's/^export defaultMem=4096.*/export defaultMem=16384/' "$SMARTSLURM_HOME/config/config.txt" > "$HOME/.smartSlurm/config/config.txt"
printf 'export myLegacyKnob=42\n' >> "$HOME/.smartSlurm/config/config.txt"
printf 'y\n' | smartslurm config migrate >/dev/null 2>&1
echo "MIG_DM=$(smartslurm config get defaultMem)"
echo "MIG_KNOB=$(smartslurm config get myLegacyKnob)"
echo "UCKEYS=$(sed -n 's/^export \([A-Za-z_0-9]*\)=.*/\1/p' "$HOME/.smartSlurm/config/user.conf" | sort | tr '\n' ',')"
[ -f "$HOME/.smartSlurm/config/config.txt.pre-overlay.bak" ] && echo "BACKEDUP=yes"
[ -f "$HOME/.smartSlurm/config/config.txt" ] || echo "LEGACY_GONE=yes"
SCN
  assert_contains "$o" "MIG_DM=16384"             "migrate keeps the changed scalar"
  assert_contains "$o" "MIG_KNOB=42"              "migrate keeps a custom key"
  assert_contains "$o" "UCKEYS=defaultMem,myLegacyKnob," "user.conf holds ONLY the diff, not every repo key"
  assert_contains "$o" "BACKEDUP=yes"             "the legacy full config is backed up"
  assert_contains "$o" "LEGACY_GONE=yes"          "the legacy full config is moved aside (overlay now active)"
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
      --tmp noTmp --pipeline "$SS_PIPELINE" >rap.log 2>&1 )
  rm -f "$HERE/models/_conv.model"
  ladder=$(ss_alloc_ladder conv)
  echo "    (info) auto-requeue alloc ladder (mem:status): $ladder"
  succ=$(find "$d" -name '*conv*.success' 2>/dev/null | head -1)
  if [ -n "$succ" ]; then _pass "under-provisioned job converged to success via auto-requeue"
  else _fail "auto-requeue converges to success" "no .success; ladder=$ladder"; fi
  firstoom=$(awk -F, '$2=="conv" && $12 ~ /OOM/{print $8; exit}' "$jr")
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
# fitQC -- jobRecord regression QC: assess fits, identify + remove outliers.
# Pure engine over a seeded v2 store; no real run needed. fitQC's core needs no
# gnuplot (only --plot does), so this block does not gate on it.
# ===========================================================================
_want fitqc && {
  _begin fitqc
  ss_new_workspace "$WORK/fitqc" >/dev/null
  jr="$smartSlurmJobRecordDir/jobRecord.txt"
  FITQC="$SMARTSLURM_ROOT/bin/fitQC"
  if [ ! -x "$FITQC" ]; then
    echo "    (skipped: fitQC not present in $SMARTSLURM_ROOT/bin)"; echo
  else
  HDR='pipeline,step,ref,sample,jobID,inputSize,memUsed,memAllocated,timeUsed,timeAllocated,cores,status,user,node,epoch,isoDate,flag'
  ep=$(date +%s)
  # Rebuildable store: a clean p1/bowtie line (mem=2*size+100) with ONE planted
  # outlier (jobID 99999 at size 350); a SECOND pipeline p2 running the same tool
  # (must stay a separate group); a thin group (n=2); and an OOM row (ignored).
  seed_store() {
    { echo "$HDR"
      printf 'p1,bowtie,none,,1001,100,300,512,10,60,4,COMPLETED,u,n,%s,d,1.0.bowtie.1\n'  "$ep"
      printf 'p1,bowtie,none,,1002,200,505,600,12,60,4,COMPLETED,u,n,%s,d,1.0.bowtie.2\n'  "$ep"
      printf 'p1,bowtie,none,,1003,300,695,800,14,60,4,COMPLETED,u,n,%s,d,1.0.bowtie.3\n'  "$ep"
      printf 'p1,bowtie,none,,99999,350,2000,2200,15,60,4,COMPLETED,u,n,%s,d,1.0.bowtie.X\n' "$ep"
      printf 'p1,bowtie,none,,1005,400,905,1000,16,60,4,COMPLETED,u,n,%s,d,1.0.bowtie.5\n' "$ep"
      printf 'p1,bowtie,none,,1006,500,1100,1200,18,60,4,COMPLETED,u,n,%s,d,1.0.bowtie.6\n' "$ep"
      printf 'p1,bowtie,none,,1007,600,1310,1400,20,60,4,COMPLETED,u,n,%s,d,1.0.bowtie.7\n' "$ep"
      printf 'p2,bowtie,none,,2001,100,800,900,10,60,8,COMPLETED,u,n,%s,d,1.0.bowtie.w1\n' "$ep"
      printf 'p2,bowtie,none,,2002,200,1600,1700,12,60,8,COMPLETED,u,n,%s,d,1.0.bowtie.w2\n' "$ep"
      printf 'p2,bowtie,none,,2003,300,2400,2500,14,60,8,COMPLETED,u,n,%s,d,1.0.bowtie.w3\n' "$ep"
      printf 'p1,star,none,,4001,100,4000,4200,30,120,8,COMPLETED,u,n,%s,d,3.0.star.1\n' "$ep"
      printf 'p1,star,none,,4002,200,8000,8200,40,120,8,COMPLETED,u,n,%s,d,3.0.star.2\n' "$ep"
      printf 'p1,bowtie,none,,1008,450,88888,1000,17,60,4,OOM,u,n,%s,d,1.0.bowtie.8\n' "$ep"
    } > "$jr"
  }
  seed_store

  # 1) triage flags the outlier group (one point over the cut)
  nout=$("$FITQC" --store "$jr" | awk '$1=="p1/bowtie/none"{print $5}')
  { [ -n "$nout" ] && [ "$nout" -ge 1 ]; } \
    && _pass "triage flags p1/bowtie/none as having an outlier (#out=$nout)" \
    || _fail "triage flags the outlier group" "#out=${nout:-<none>}"

  # 2) cross-pipeline isolation: p1 and p2 bowtie are SEPARATE groups
  ngrp=$("$FITQC" --store "$jr" --all | grep -cE '^(p1|p2)/bowtie/none ')
  [ "${ngrp:-0}" -eq 2 ] \
    && _pass "same tool in two pipelines -> two separate groups (no cross-contamination)" \
    || _fail "cross-pipeline isolation" "matched $ngrp group rows, expected 2"

  # 3) detail ranks the planted outlier worst (first data row = jobID 99999)
  first=$("$FITQC" --store "$jr" p1/bowtie/none | awk '/^   [0-9]/{print $1; exit}')
  [ "$first" = 99999 ] && _pass "detail ranks the planted outlier (99999) worst" \
    || _fail "detail ranks the outlier worst" "top row jobID=$first"

  # 4) thin group (n=2) reported thin, never fit
  "$FITQC" --store "$jr" --all | awk '$1=="p1/star/none"{print $7}' | grep -q thin \
    && _pass "thin group (n=2) reported thin, not fit" || _fail "thin group reported thin"

  # 5) --dry-run --remove is a no-op on the real store
  rb=$(wc -l < "$jr")
  "$FITQC" --store "$jr" p1/bowtie/none --remove 99999 --dry-run >/dev/null
  ra=$(wc -l < "$jr")
  { [ "$rb" = "$ra" ] && grep -q ',99999,' "$jr"; } \
    && _pass "--dry-run --remove writes nothing (store unchanged)" \
    || _fail "--dry-run --remove no-op" "rows $rb->$ra"

  # 6) real --remove drops the outlier, backs up, and the re-fit R2 jumps
  r2b=$("$FITQC" --store "$jr" p1/bowtie/none | sed -n 's/.*R2=\([0-9.]*\).*/\1/p' | head -1)
  "$FITQC" --store "$jr" p1/bowtie/none --remove 99999 >/dev/null
  r2a=$("$FITQC" --store "$jr" p1/bowtie/none | sed -n 's/.*R2=\([0-9.]*\).*/\1/p' | head -1)
  gone=$(awk -F, 'NR>1 && $5=="99999"{c++} END{print c+0}' "$jr")
  if [ "$gone" -eq 0 ] \
     && ls "$smartSlurmJobRecordDir"/jobRecord.back.*.txt >/dev/null 2>&1 \
     && awk -v b="$r2b" -v a="$r2a" 'BEGIN{exit !(b+0<0.9 && a+0>0.99)}'; then
    _pass "--remove drops the outlier + backs up; R2 $r2b -> $r2a"
  else _fail "--remove cleans the group" "99999 gone? backup? R2 $r2b->$r2a"; fi

  # 7) iterative --auto removes the masked outlier and converges
  seed_store
  auto=$("$FITQC" --store "$jr" p1/bowtie/none --auto)
  if echo "$auto" | grep -q 'removed jobID 99999' \
     && echo "$auto" | grep -q converged \
     && [ "$(awk -F, 'NR>1 && $5=="99999"{c++} END{print c+0}' "$jr")" -eq 0 ]; then
    _pass "--auto removes the masked outlier and converges"
  else _fail "--auto removes + converges" "$(echo "$auto" | grep -iE 'pass|converg|stop' | head -1)"; fi

  # 8) --auto refuses a HIGH-LEVERAGE worst point (never auto-remove extremes)
  { echo "$HDR"
    printf 'pL,tool,none,,5001,100,350,400,10,60,1,COMPLETED,u,n,%s,d,1.0.tool.a\n'  "$ep"
    printf 'pL,tool,none,,5002,200,650,700,10,60,1,COMPLETED,u,n,%s,d,1.0.tool.b\n'  "$ep"
    printf 'pL,tool,none,,5003,300,950,1000,10,60,1,COMPLETED,u,n,%s,d,1.0.tool.c\n' "$ep"
    printf 'pL,tool,none,,5004,400,1250,1300,10,60,1,COMPLETED,u,n,%s,d,1.0.tool.d\n' "$ep"
    printf 'pL,tool,none,,5005,500,1550,1600,10,60,1,COMPLETED,u,n,%s,d,1.0.tool.e\n' "$ep"
    printf 'pL,tool,none,,5006,2000,3000,9000,10,60,1,COMPLETED,u,n,%s,d,1.0.tool.f\n' "$ep"
  } > "$jr"
  guard=$("$FITQC" --store "$jr" pL/tool/none --auto --resid-cut 1.5)
  if echo "$guard" | grep -q 'HIGH-LEVERAGE' && grep -q ',5006,' "$jr"; then
    _pass "--auto refuses a high-leverage worst point (kept 5006)"
  else _fail "--auto high-leverage guard" "$(echo "$guard" | grep -iE 'stop|high|converg' | head -1)"; fi

  # 9) non-blocking warning when one group regresses >1 rAP step together
  #    (same step name + ref, different step NUMBERS in the flag). A single-step
  #    group must NOT trigger it.
  { echo "$HDR"
    printf 'pM,bowtie,none,,601,100,300,512,10,60,4,COMPLETED,u,n,%s,d,2.0.bowtie.a\n' "$ep"
    printf 'pM,bowtie,none,,602,200,500,600,12,60,4,COMPLETED,u,n,%s,d,2.0.bowtie.b\n' "$ep"
    printf 'pM,bowtie,none,,603,300,700,800,14,60,4,COMPLETED,u,n,%s,d,3.0.bowtie.c\n' "$ep"
    printf 'pM,bowtie,none,,604,400,900,1000,16,60,4,COMPLETED,u,n,%s,d,3.0.bowtie.d\n' "$ep"
    printf 'pM,cut,none,,611,100,150,256,5,30,1,COMPLETED,u,n,%s,d,1.0.cut.a\n' "$ep"
    printf 'pM,cut,none,,612,200,250,300,6,30,1,COMPLETED,u,n,%s,d,1.0.cut.b\n' "$ep"
    printf 'pM,cut,none,,613,300,360,400,7,30,1,COMPLETED,u,n,%s,d,1.0.cut.c\n' "$ep"
  } > "$jr"
  tri=$("$FITQC" --store "$jr" --all)
  merged_note=$(echo "$tri" | awk '$1=="pM/bowtie/none"{print $NF}')
  single_note=$(echo "$tri" | awk '$1=="pM/cut/none"{print $NF}')
  det=$("$FITQC" --store "$jr" pM/bowtie/none)
  if echo "$merged_note" | grep -q 'steps:2,3' \
     && [ "$single_note" = "-" ] \
     && echo "$det" | grep -q 'WARNING: this group regresses'; then
    _pass "merged-steps warning fires for 2-step group, silent for single-step"
  else _fail "merged-steps warning" "merged=[$merged_note] single=[$single_note]"; fi

  # 10) low-R2 surfacing: a scattered group (R2<0.80) with NO point over the cut is
  #     shown in DEFAULT triage (flagged low-R2) -- masking means the cut alone can
  #     miss a bad fit. A clean (R2~1) group stays hidden by default.
  { echo "$HDR"
    printf 'pR,aln,none,,801,100,500,600,10,60,4,COMPLETED,u,n,%s,d,2.0.aln.a\n' "$ep"
    printf 'pR,aln,none,,802,200,400,600,10,60,4,COMPLETED,u,n,%s,d,2.0.aln.b\n' "$ep"
    printf 'pR,aln,none,,803,300,900,600,10,60,4,COMPLETED,u,n,%s,d,2.0.aln.c\n' "$ep"
    printf 'pR,aln,none,,804,400,700,600,10,60,4,COMPLETED,u,n,%s,d,2.0.aln.d\n' "$ep"
    printf 'pR,aln,none,,805,500,1200,600,10,60,4,COMPLETED,u,n,%s,d,2.0.aln.e\n' "$ep"
    printf 'pR,aln,none,,806,600,800,600,10,60,4,COMPLETED,u,n,%s,d,2.0.aln.f\n' "$ep"
    printf 'pR,clean,none,,811,100,300,512,10,60,4,COMPLETED,u,n,%s,d,1.0.clean.a\n' "$ep"
    printf 'pR,clean,none,,812,200,500,600,12,60,4,COMPLETED,u,n,%s,d,1.0.clean.b\n' "$ep"
    printf 'pR,clean,none,,813,300,700,800,14,60,4,COMPLETED,u,n,%s,d,1.0.clean.c\n' "$ep"
  } > "$jr"
  tri=$("$FITQC" --store "$jr")          # default triage (no --all)
  aln=$(echo "$tri" | awk '$1=="pR/aln/none"')
  aln_out=$(echo "$aln" | awk '{print $5}')
  aln_note=$(echo "$aln" | awk '{print $NF}')
  clean_shown=$(echo "$tri" | grep -c '^pR/clean/none ')
  if echo "$aln_note" | grep -q 'low-R2' && [ "${aln_out:-x}" = 0 ] && [ "$clean_shown" -eq 0 ]; then
    _pass "low-R2 group surfaced by default (no over-cut outlier); clean group stays hidden"
  else _fail "low-R2 surfacing" "aln=[$aln] cleanShown=$clean_shown"; fi

  echo
  fi
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
# jobRecords v1 -> v2 auto-migration (the data-loss regression from ld32's real
# cluster: appending v2 rows to an un-migrated v1 store produced a MIXED file
# whose v2 rows were then deleted by the v1-positional prune and ignored by the
# v1 reader). Guards the fix in bin/jobRecordSchema.sh (jrMigrateStore/jrEmitV2/
# jrSchemaVersion), bin/cleanUp.sh (migrate-before-prune) and bin/jrRepair.sh.
# Uses the REAL store layouts seen in the wild: headerless, 21-col (trailing
# epoch) and the older 20-col variant (date, no epoch).
_want migrate_v1 && {
  _begin migrate_v1
  sch="$SMARTSLURM_ROOT/bin/jobRecordSchema.sh"
  if [ ! -r "$sch" ]; then
    _pass "(skip) jobRecordSchema.sh absent — SmartSlurm has no v2 schema on this base"
  else
  ss_new_workspace "$WORK/migrate_v1"
  jr="$smartSlurmJobRecordDir/jobRecord.txt"
  jrver()   { ( . "$sch"; jrSchemaVersion "$1" ); }
  jrmig()   { ( . "$sch"; jrMigrateStore "$1" ) 2>/dev/null; }
  jremit()  { ( . "$sch"; jrEmitV2 "$1" ); }
  jrrec()   { ( . "$sch"; jrCountRecognized "$1" ); }
  # 21-col v1 row (…,date,epoch); 20-col v1 row (…,date — NO epoch).
  v1row21() { printf '%s,%s,0,0,4096,120,%s,95,%s,tester,,%s,%s,%s,1,0,0,0,0,2026-01-01,%s\n' \
                "$1" "$2" "$3" "${6:-COMPLETED}" "$4" "${5:-none}" "1.0.$4.$1" "$(( 1735700000 + $1 ))"; }
  v1row20() { printf '%s,%s,0,0,4096,120,%s,95,%s,tester,,%s,%s,%s,1,0,0,0,0,2026-01-01\n' \
                "$1" "$2" "$3" "${6:-COMPLETED}" "$4" "${5:-none}" "1.0.$4.$1"; }

  # --- 1) headerless mixed-width v1 store migrates with zero loss + backup ------
  { v1row21 1001 5000 3800 bwa hg38; v1row20 1002 6000 7900 bowtie2 mm10; } > "$jr"
  [ "$(jrver "$jr")" = 1 ] && _pass "headerless 20+21-col store detected as v1" \
                           || _fail "headerless 20+21-col store detected as v1" "got $(jrver "$jr")"
  jrmig "$jr"
  [ "$(jrver "$jr")" = 2 ] && _pass "store migrated to v2" || _fail "store migrated to v2" "got $(jrver "$jr")"
  assert_contains "$jr" "pipeline,step,ref"            "migrated store has the v2 header"
  # both legacy jobIDs survive (by column, not substring, so size/mem can't alias)
  g1=$(awk -F, '$5==1001' "$jr" | wc -l); g2=$(awk -F, '$5==1002' "$jr" | wc -l)
  [ "$g1" -eq 1 ] && _pass "21-col legacy row 1001 preserved" || _fail "21-col legacy row 1001 preserved" "count $g1"
  [ "$g2" -eq 1 ] && _pass "20-col legacy row 1002 preserved" || _fail "20-col legacy row 1002 preserved" "count $g2"
  # 20-col row has no epoch -> epoch column blank; 21-col row keeps its epoch
  assert_contains "$jr" "legacy,bwa,hg38,"             "21-col row mapped (pipeline=legacy, step, ref)"
  e20=$(awk -F, '$5==1002{print $15}' "$jr"); [ -z "$e20" ] && _pass "20-col row epoch left blank (no epoch in source)" \
                                                          || _fail "20-col row epoch left blank" "got '$e20'"
  # 21-col row's epoch MUST survive migration (regression guard: a {n,} interval
  # regex silently blanked it on mawk). v1row21 sets epoch = 1735700000 + jobID.
  e21=$(awk -F, '$5==1001{print $15}' "$jr"); [ "$e21" = 1735701001 ] && _pass "21-col row epoch preserved (mawk-safe: 1735701001)" \
                                                          || _fail "21-col row epoch preserved" "got '$e21'"
  ls "$jr".back.*.txt >/dev/null 2>&1 && _pass "original backed up before migration" \
                                       || _fail "original backed up before migration" "no .back.*.txt"

  # --- 2) idempotent: a second migrate is a no-op, makes no new backup ----------
  nb1=$(ls "$jr".back.*.txt 2>/dev/null | wc -l); sleep 1; jrmig "$jr"
  nb2=$(ls "$jr".back.*.txt 2>/dev/null | wc -l)
  [ "$(jrver "$jr")" = 2 ] && [ "$nb1" -eq "$nb2" ] && _pass "re-migrating a v2 store is a no-op (no extra backup)" \
                                                     || _fail "re-migrating a v2 store is a no-op" "backups $nb1 -> $nb2"

  # --- 3) the exact loss vector: a MIXED file (v1 rows + an appended v2 row) ----
  # reproduces what the OLD writer produced on an un-migrated store. Migration
  # must keep ALL of them -- the v1 rows mapped and the already-v2 row verbatim.
  ss_new_workspace "$WORK/migrate_v1_mixed"; jr="$smartSlurmJobRecordDir/jobRecord.txt"
  { v1row21 2001 5000 3800 bwa hg38; v1row20 2002 6000 7900 bwa hg38
    printf 'rnaseq,align,hg38,sampC,2003,7000,9000,10000,300,600,4,COMPLETED,tester,,1735800000,2026-01-03,4.align.sampC\n'; } > "$jr"
  jrmig "$jr"
  [ "$(jrver "$jr")" = 2 ] && _pass "mixed (v1+appended-v2) store migrates to v2" \
                           || _fail "mixed store migrates to v2" "got $(jrver "$jr")"
  dat=$(( $(grep -c . "$jr") - 1 ))
  [ "$dat" -eq 3 ] && _pass "mixed store: all 3 data rows preserved (no v1-prune loss)" \
                    || _fail "mixed store: all 3 rows preserved" "got $dat"
  v2keep=$(awk -F, '$5==2003 && $1=="rnaseq"' "$jr" | wc -l)
  [ "$v2keep" -eq 1 ] && _pass "appended v2 row passed through unchanged" \
                       || _fail "appended v2 row passed through unchanged" "count $v2keep"

  # --- 4) jrRepair diverts ONLY genuine junk (not v1 or v2 rows) ----------------
  ss_new_workspace "$WORK/migrate_v1_repair"; jr="$smartSlurmJobRecordDir/jobRecord.txt"
  { v1row20 3001 5000 3800 bwa hg38
    printf 'GARBAGE,not,a,record\n'
    printf 'unknown,trim,none,,3002,900,700,800,60,90,1,COMPLETED,tester,,1735900000,2026-01-04,1.0.trim.1\n'; } > "$jr"
  bash "$SMARTSLURM_ROOT/bin/jrRepair.sh" "$smartSlurmJobRecordDir" >/dev/null 2>&1
  [ "$(jrver "$jr")" = 2 ] && _pass "jrRepair converts a mixed store to v2" \
                           || _fail "jrRepair converts a mixed store to v2" "got $(jrver "$jr")"
  sk=$(ls "$jr".jrRepair.skipped.* 2>/dev/null | head -1)
  if [ -n "$sk" ] && grep -q GARBAGE "$sk" && [ "$(wc -l < "$sk")" -eq 1 ]; then
    _pass "jrRepair diverts exactly the 1 junk row to .skipped"
  else _fail "jrRepair diverts exactly the junk row" "skipped=$(cat "$sk" 2>/dev/null)"; fi
  assert_lacks "$sk" "3001" "jrRepair did NOT divert the 20-col v1 row"
  assert_lacks "$sk" "3002" "jrRepair did NOT divert the already-v2 row"

  # --- 5) HYBRID the dangerous way: v1 rows appended UNDER a v2 header ----------
  # (an old node writing v1 into a store a new node already migrated, mid-rollout)
  # A line-1-only schema check called this "pure v2" and skipped migration, then
  # the v2 prune read the v1 rows through the wrong columns and DELETED COMPLETED
  # history. jrSchemaVersion must treat any v1 row as "needs migration".
  # v2row JOBID STEP REF INPUTSIZE MEMUSED [SAMPLE] -> one canonical 17-col v2 row
  v2row() { local id=$1 st=$2 rf=$3 sz=$4 mu=$5 sm=${6:-}
    printf 'rnaseq,%s,%s,%s,%s,%s,%s,%s,95,120,1,COMPLETED,tester,n,%s,2026-01-05,1.0.%s.%s\n' \
      "$st" "$rf" "$sm" "$id" "$sz" "$mu" "$mu" "$(( 1735900000 + id ))" "$st" "$id"; }
  hdr2='pipeline,step,ref,sample,jobID,inputSize,memUsed,memAllocated,timeUsed,timeAllocated,cores,status,user,node,epoch,isoDate,flag'
  # v1-prune-sim: replicate cleanUp's EXACT branch selection, then assert no loss.
  prunesim() { local f=$1 t; t=$(date -d '7 days ago' +%s)
    if [ "$(jrver "$f")" = 2 ]; then
      awk -F, -v t="$t" 'NR==1{for(i=1;i<=NF;i++)c[$i]=i;print;next} ($(c["epoch"])+0>t)||($(c["status"])=="COMPLETED")' "$f" > "$f.n" && mv "$f.n" "$f"
    else
      awk -F, -v t="$t" '{if($21>t||$9=="COMPLETED")print}' "$f" > "$f.n" && mv "$f.n" "$f"
    fi; }
  for order in v2first v1first; do
    ss_new_workspace "$WORK/migrate_v1_hybrid_$order"; jr="$smartSlurmJobRecordDir/jobRecord.txt"
    if [ "$order" = v2first ]; then { echo "$hdr2"; v2row 801 bwa hg38 5000 3800 sampA; v1row21 802 6000 7900 bwa hg38; } > "$jr"
    else { echo "$hdr2"; v1row20 901 6000 7900 bwa hg38; v2row 902 bwa hg38 5000 3800 sampB; } > "$jr"; fi
    [ "$(jrver "$jr")" = 1 ] && _pass "hybrid($order): v2-headed store with a v1 row flagged for migration" \
                             || _fail "hybrid($order): flagged for migration" "jrver=$(jrver "$jr")"
    jrmig "$jr"; prunesim "$jr"
    n=$(awk -F, 'NR>1 && NF==17 && $1!~/^[0-9]+$/' "$jr" | wc -l)
    [ "$n" -eq 2 ] && _pass "hybrid($order): both COMPLETED rows survive migrate+prune (no self-destruct)" \
                   || _fail "hybrid($order): both rows survive" "survivors=$n"
    norawv1=$(awk -F, 'NR>1 && $1~/^[0-9]+$/ && (NF==20||NF==21)' "$jr" | wc -l)
    [ "$norawv1" -eq 0 ] && _pass "hybrid($order): no raw v1 row left for the prune to mis-read" \
                         || _fail "hybrid($order): raw v1 row remains" "count=$norawv1"
  done

  # --- 6) headerless v2 store gets a leading header (prune eats line 1 as hdr) --
  ss_new_workspace "$WORK/migrate_v1_hlv2"; jr="$smartSlurmJobRecordDir/jobRecord.txt"
  { v2row 1001 bwa hg38 5000 3800 sampA; v2row 1002 bwa hg38 6000 4200 sampB; } > "$jr"
  jrmig "$jr"
  assert_contains "$jr" "pipeline,step,ref" "headerless v2 store gets a leading v2 header on migrate"
  hv=$(awk -F, '$5==1001 && $1=="rnaseq"' "$jr" | wc -l)
  [ "$hv" -eq 1 ] && _pass "headerless v2: first data row preserved (not eaten as header)" \
                  || _fail "headerless v2: first row preserved" "count=$hv"

  # --- 7) NON-SILENT diagnostics: malformed rows must be reported, not dropped --
  # quietly. A store with a genuinely malformed row (unknown column count) must:
  # warn on stderr naming the count, keep the row in the backup, and be surfaced
  # by `jrRepair.sh --check` (schema + unrecognized rows, non-zero exit).
  ss_new_workspace "$WORK/migrate_v1_nonsilent"; jr="$smartSlurmJobRecordDir/jobRecord.txt"
  { v1row21 5001 5000 3800 bwa hg38
    printf '5002,5000,0,0,4096,120,3800,95,COMPLETED,u,,bwa,hg38,1.0.bwa.s2,1,0,0\n'   # 19-col: malformed
    v2row 5003 bwa hg38 6000 4000 sampC; } > "$jr"
  warn=$( ( . "$sch"; jrMigrateStore "$jr" ) 2>&1 1>/dev/null )
  echo "$warn" | grep -qiE 'WARNING.*1 row.*no known schema' && _pass "migration WARNS about the dropped malformed row (not silent)" \
                                                             || _fail "migration warns about dropped row" "stderr: $warn"
  bk=$(ls "$jr".back.*.txt 2>/dev/null | head -1)
  if [ -n "$bk" ] && awk -F, '$1==5002' "$bk" | grep -q .; then _pass "malformed row preserved in backup (recoverable)"
  else _fail "malformed row preserved in backup" "backup=$bk"; fi
  uf=$(ls "$jr".unrecognized.*.txt 2>/dev/null | head -1)
  if [ -n "$uf" ] && grep -q 5002 "$uf" && [ "$(grep -c . "$uf")" -eq 1 ]; then
    _pass "dropped row also saved to a .unrecognized file (just the bad line)"
  else _fail "dropped row saved to .unrecognized file" "file=$uf"; fi
  noraw=$(awk -F, 'NR>1 && NF!=17' "$jr" | grep -c .)
  [ "$noraw" -eq 0 ] && _pass "live store left clean v2 (malformed row not carried in)" \
                     || _fail "live store clean after migrate" "stray rows=$noraw"
  # --- integrity gate: a broken conversion aborts and leaves the store untouched -
  ss_new_workspace "$WORK/migrate_v1_gate"; jr="$smartSlurmJobRecordDir/jobRecord.txt"
  v1row21 6001 5000 3800 bwa hg38 > "$jr"; v1row21 6002 6000 4000 bwa hg38 >> "$jr"
  g_before=$(md5sum "$jr" | cut -d' ' -f1)
  gerr=$( ( . "$sch"
            jrEmitV2() { printf '%s\n' "$JR_HEADER_V2"; echo "rnaseq,bwa,hg38,,6001,5000,3800,4096,95,120,1,COMPLETED,u,n,1,2026-01-05,1.0.bwa.s"; }  # drops 6002
            jrMigrateStore "$jr" ) 2>&1 1>/dev/null )
  if [ "$(md5sum "$jr" | cut -d' ' -f1)" = "$g_before" ] && echo "$gerr" | grep -qiE 'ERROR.*integrity.*UNTOUCHED'; then
    _pass "integrity gate: a conversion that loses a recognizable row aborts, store untouched + errors"
  else _fail "integrity gate aborts a lossy conversion" "changed=$([ "$(md5sum "$jr"|cut -d' ' -f1)" = "$g_before" ] && echo no || echo YES) err=$gerr"; fi
  [ -z "$(ls "$jr".back.*.txt 2>/dev/null)" ] && _pass "integrity gate: redundant backup removed on abort" \
                                              || _fail "integrity gate: backup removed on abort" "backup left"
  # jrRepair --check surfaces an unrecognized store (seed a fresh junk store)
  printf 'total,garbage\n99,bad,row,here\n' > "$smartSlurmJobRecordDir/jobRecord.oom.oot.over.txt"
  chk_out=$(bash "$SMARTSLURM_ROOT/bin/jrRepair.sh" --check "$smartSlurmJobRecordDir" 2>&1); chk_rc=$?
  echo "$chk_out" | grep -q 'UNRECOGNIZED' && [ "$chk_rc" -ne 0 ] \
    && _pass "jrRepair --check reports UNRECOGNIZED rows and exits non-zero" \
    || _fail "jrRepair --check surfaces unrecognized store" "rc=$chk_rc out=$chk_out"
  echo "$chk_out" | grep -qE 'schema version' && _pass "jrRepair --check prints a schema/row census" \
                                              || _fail "jrRepair --check prints a census" "out=$chk_out"

  # --- 8) property/fuzz harness: strong invariants over randomized shapes -------
  if [ -x "$HERE/jr_conversion_stress.sh" ]; then
    for seed in 12345 2026; do
      if SMARTSLURM_ROOT="$SMARTSLURM_ROOT" "$HERE/jr_conversion_stress.sh" "$seed" >/dev/null 2>&1; then
        _pass "conversion stress/fuzz harness all-green (seed $seed)"
      else _fail "conversion stress/fuzz harness (seed $seed)" "see test/jr_conversion_stress.sh $seed"; fi
    done
  else _pass "(skip) stress harness not present"; fi

  # --- 8b) migration coverage gaps: migrate-time diversion, lock reclaim/skip, -----
  # crash-temp tolerance, and estimate-equivalence across migration (paths the fuzz
  # had code for but no assertion -- each check is teeth-proven by mutation).
  if [ -x "$HERE/jr_migration_gaps_guard.sh" ]; then
    if SMARTSLURM_ROOT="$SMARTSLURM_ROOT" "$HERE/jr_migration_gaps_guard.sh" >/dev/null 2>&1; then
      _pass "migration-gaps guard all-green"
    else _fail "migration-gaps guard" "see test/jr_migration_gaps_guard.sh"; fi
  else _pass "(skip) migration-gaps guard not present"; fi

  # --- 8c) mawk enforcement: re-run the migration suites under mawk (Debian/Ubuntu ---
  # default awk). mawk ignores {n,} interval regexes, which would silently blank every
  # migrated epoch; gawk here would never catch that. Gated on mawk being installed.
  if command -v mawk >/dev/null 2>&1; then
    _mawkd=$(mktemp -d); ln -sf "$(command -v mawk)" "$_mawkd/awk"
    _mawk_ok=1
    PATH="$_mawkd:$PATH" SMARTSLURM_ROOT="$SMARTSLURM_ROOT" "$HERE/jr_conversion_stress.sh" 12345 >/dev/null 2>&1 || _mawk_ok=0
    if [ -x "$HERE/jr_migration_gaps_guard.sh" ]; then
      PATH="$_mawkd:$PATH" SMARTSLURM_ROOT="$SMARTSLURM_ROOT" "$HERE/jr_migration_gaps_guard.sh" >/dev/null 2>&1 || _mawk_ok=0
    fi
    rm -rf "$_mawkd"
    [ "$_mawk_ok" = 1 ] && _pass "migration suites all-green under mawk" \
                        || _fail "migration suites under mawk" "re-run with awk->mawk; see jr_conversion_stress.sh / jr_migration_gaps_guard.sh"
  else _pass "(skip) mawk not installed -- cannot enforce mawk parity"; fi

  # --- 9) pipeline literally named "pipeline" (runAsPipeline's default id) -------
  # The v2 header was detected by $1=="pipeline" alone, but that is also the
  # generic pipeline id runAsPipeline assigns when --pipeline is omitted. So every
  # data row of a default-named pipeline looked like a header line and was dropped
  # by jrEmitV2 / ignored by jrSchemaVersion / jrCountRecognized -- silently
  # erasing a whole pipeline's history (and blinding its resource estimation).
  # The real header's 2nd field is the literal "step"; a data row's 2nd field is
  # the step name. These must round-trip with ZERO loss.  (Pre-fix: jremit would
  # emit only the header and jrrec would count 0; both FAIL here.)
  # ppl ROW: a canonical 17-col v2 row whose PIPELINE is "pipeline".
  ppl() { local id=$1 st=$2 sz=$3 mu=$4 status=${5:-COMPLETED}
    printf 'pipeline,%s,none,,%s,%s,%s,%s,30,60,1,%s,tester,n,%s,2026-01-06,1.0.%s\n' \
      "$st" "$id" "$sz" "$mu" "$mu" "$status" "$(( 1735900000 + id ))" "$st"; }
  ss_new_workspace "$WORK/migrate_v1_pplname"; jr="$smartSlurmJobRecordDir/jobRecord.txt"
  { echo "$hdr2"; ppl 7101 mseq 5000 4096 OOM; ppl 7102 mseq 5000 9840 OOM
    ppl 7103 mseq 5000 150000 COMPLETED; } > "$jr"
  # header-led store of pipeline="pipeline" rows is canonical v2 (a migration no-op)
  [ "$(jrver "$jr")" = 2 ] && _pass "pplname: header-led pipeline=\"pipeline\" store detected as canonical v2" \
                           || _fail "pplname: detected as v2" "got $(jrver "$jr")"
  # the crux: all 3 data rows are RECOGNIZED, not misread as header lines
  rec=$(jrrec "$jr")
  [ "$rec" -eq 3 ] && _pass "pplname: all 3 pipeline=\"pipeline\" rows counted (not skipped as header)" \
                   || _fail "pplname: rows counted as recognized" "jrCountRecognized=$rec (pre-fix bug counts 0)"
  # jrEmitV2 preserves every row (header + 3), not just the header line
  em=$(jremit "$jr" | grep -c .)
  [ "$em" -eq 4 ] && _pass "pplname: jrEmitV2 preserves header + all 3 rows (no header-collision drop)" \
                  || _fail "pplname: jrEmitV2 preserves all rows" "emitted $em lines (pre-fix bug emits 1)"
  # by-column survival: each jobID present exactly once after an emit round-trip
  for id in 7101 7102 7103; do
    n=$(jremit "$jr" | awk -F, -v j="$id" '$5==j && $1=="pipeline"' | wc -l)
    [ "$n" -eq 1 ] && _pass "pplname: row $id survives emit (pipeline col intact)" \
                   || _fail "pplname: row $id survives emit" "count=$n"
  done
  # headerless pipeline="pipeline" store: the dangerous path. A leading v2 DATA row
  # with pipeline="pipeline" must NOT be mistaken for a canonical v2 header (which
  # would skip migration and let the header-eating prune drop line 1). It must flag
  # for migration, then migrate to a header-led store with all rows intact.
  ss_new_workspace "$WORK/migrate_v1_pplname_hl"; jr="$smartSlurmJobRecordDir/jobRecord.txt"
  { ppl 7201 mseq 5000 4096 OOM; ppl 7202 mseq 5000 150000 COMPLETED; } > "$jr"
  [ "$(jrver "$jr")" = 1 ] && _pass "pplname/headerless: leading pipeline=\"pipeline\" row flagged for migration" \
                           || _fail "pplname/headerless: flagged for migration" "jrver=$(jrver "$jr") (pre-fix bug says 2 -> no migrate)"
  jrmig "$jr"
  assert_contains "$jr" "pipeline,step,ref" "pplname/headerless: gets a leading v2 header on migrate"
  dat=$(awk -F, 'NR>1 && $1=="pipeline" && NF==17' "$jr" | wc -l)
  [ "$dat" -eq 2 ] && _pass "pplname/headerless: both pipeline=\"pipeline\" data rows preserved on migrate" \
                   || _fail "pplname/headerless: both rows preserved" "survivors=$dat"

  unset -f jrver jrmig jremit jrrec v1row21 v1row20 v2row prunesim ppl
  fi
  echo
}

# ===========================================================================
# Schema-drift guard: proves no tool hand-rolls the jobRecord format (the root
# cause of the original 18-vs-21-col drift). Runs the committed guard against the
# SmartSlurm checkout, and self-tests that the guard actually catches a violation
# (a guard that can't fail is worthless).
_want jr_drift_guard && {
  _begin jr_drift_guard
  if [ ! -r "$SMARTSLURM_ROOT/bin/jobRecordSchema.sh" ]; then
    _pass "(skip) jobRecordSchema.sh absent — the v2 I/O API does not exist on this base"
  elif [ -x "$HERE/jr_drift_guard.sh" ]; then
    if SMARTSLURM_ROOT="$SMARTSLURM_ROOT" "$HERE/jr_drift_guard.sh" >/dev/null 2>&1; then
      _pass "no tool hand-rolls the jobRecord format (all store I/O via the schema API)"
    else
      _fail "schema-drift guard" "$(SMARTSLURM_ROOT="$SMARTSLURM_ROOT" "$HERE/jr_drift_guard.sh" 2>&1 | head -12)"
    fi
    # self-test: a bin/ with a raw append + header literal MUST trip the guard
    gd=$(mktemp -d); mkdir "$gd/bin"
    cp "$SMARTSLURM_ROOT/bin/jobRecordSchema.sh" "$gd/bin/" 2>/dev/null
    { echo '#!/bin/bash'
      echo 'echo "$r" >> $smartSlurmJobRecordDir/jobRecord.txt'
      echo 'echo 1jobID,2inputSize,3memDefault > $smartSlurmJobRecordDir/jobRecord.txt'; } > "$gd/bin/evil.sh"
    if SMARTSLURM_ROOT="$gd" "$HERE/jr_drift_guard.sh" >/dev/null 2>&1; then
      _fail "guard self-test: catches a planted violation" "guard passed a file with a raw append + header literal"
    else
      _pass "guard self-test: a planted raw-append + header-literal is caught"
    fi
    rm -rf "$gd"
  else _pass "(skip) drift guard not present"; fi
  echo
}

# ===========================================================================
ss_report
[ "$FAILS" -eq 0 ] && exit 0 || exit 1
