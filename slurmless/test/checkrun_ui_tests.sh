#!/bin/bash
# ===========================================================================
# checkrun_ui_tests.sh -- UX regression tests for the checkRun rework (M3).
#
#   ./test/checkrun_ui_tests.sh              run all UI blocks
#   ./test/checkrun_ui_tests.sh grid         run only blocks matching "grid"
#
# Builds REAL fixtures with runAsPipeline under the mock scheduler (same M0
# harness as run_tests.sh), then exercises:
#   * the checkRunGrid ENGINE  -- framed grid render + address/number resolve
#   * the checkRun CONTROLLER   -- grid-as-home, <verb> <addr> dispatch, list
#     view, raw browser, and in-place (non-flooding) redraw -- driven over a
#     pseudo-terminal because checkRun reads from /dev/tty.
#
# Exit status is non-zero if any assertion fails (CI-friendly).
# ===========================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/harness.sh
source "$HERE/lib/harness.sh"
ss_activate

FILTER="${1:-}"
WORK="${SS_SCRATCH:-/tmp/ss-scratch}/work"
_want() { [ -z "$FILTER" ] || [[ "$1" == *"$FILTER"* ]]; }

# --- stubs + pty driver so the interactive controller never blocks ----------
STUBS="${SS_SCRATCH:-/tmp/ss-scratch}/ui-stubs"
mkdir -p "$STUBS"
printf '#!/bin/sh\nexec cat "$@"\n' > "$STUBS/less"      # pager -> dump, no block
printf '#!/bin/sh\nexit 0\n'        > "$STUBS/display"    # X viewer -> no-op
chmod +x "$STUBS"/*

PTY="$STUBS/pty_drive.py"
cat > "$PTY" <<'PYEOF'
import os, pty, select, time, sys, re
workdir, extra_path = sys.argv[1], sys.argv[2]
cmd = sys.argv[3:]
inputs = [l.rstrip("\n") for l in sys.stdin]     # each line: "TEXT" or "DELAY\tTEXT"
pid, fd = pty.fork()
if pid == 0:
    env = dict(os.environ)
    env["PATH"] = extra_path + ":" + env.get("PATH", "")
    env["TERM"] = "xterm-256color"
    os.chdir(workdir); os.execvpe(cmd[0], cmd, env); os._exit(127)
buf = []
def drain(t):
    end = time.time() + t
    while time.time() < end:
        r, _, _ = select.select([fd], [], [], 0.15)
        if r:
            try: d = os.read(fd, 65536)
            except OSError: return
            if not d: return
            buf.append(d.decode(errors="replace")); end = time.time() + t
drain(0.8)
for line in inputs:
    delay, text = ("0.4", line)
    if "\t" in line: delay, text = line.split("\t", 1)
    if text.startswith("@WAIT"):                # observe output without sending keys
        secs = float(text.split(":", 1)[1]) if ":" in text else 4.0
        end = time.time() + secs
        while time.time() < end: drain(0.3)
        continue
    os.write(fd, (text + "\n").encode()); time.sleep(float(delay)); drain(0.5)
drain(0.8)
try: os.close(fd)
except OSError: pass
try: os.waitpid(pid, 0)
except Exception: pass
raw = "".join(buf)
# metrics to stderr; ANSI-stripped transcript to stdout
sys.stderr.write("##esc_home=%d\n"  % raw.count("\x1b[H"))
sys.stderr.write("##esc_clear=%d\n" % (raw.count("\x1b[2J") + raw.count("\x1b[3J")))
sys.stdout.write(re.sub(r"\x1b\[[0-9;?]*[A-Za-z]", "", raw))
PYEOF

# drive_checkrun OUTFILE  <- inputs on stdin ; transcript->OUTFILE, metrics->OUTFILE.meta
drive_checkrun() {
  python3 "$PTY" "$PWD" "$STUBS" bash "$SMARTSLURM_ROOT/bin/checkRun" > "$1" 2> "$1.meta"
}

# ===========================================================================
_want grid && {
  echo "== checkRunGrid engine: framed grid render + address resolve =="
  ss_new_workspace "$WORK/ui_grid"
  cp "$TEST_DIR"/fixtures/clean_2step/* .; chmod +x ./*.sh
  ss_run_pipeline "./pipeline.sh 777" "sbatch -p short -c 1 --mem 4G -t 50:0" >/dev/null 2>&1

  checkRunGrid --render --no-color --no-hint --log smartSlurmLog > grid.txt 2>&1
  if grep -qE 'sample .*row' grid.txt; then _pass "grid: sample/row header present"; else _fail "grid: sample/row header present"; fi
  assert_contains grid.txt "findNum"     "grid: step NAME row present (truncated)"
  assert_contains grid.txt "┌"           "grid: top frame line drawn"
  assert_contains grid.txt "│"           "grid: cell frame drawn"
  assert_contains grid.txt "└"           "grid: bottom frame line drawn"
  # numbers row sits on its own line (just the step number, right above the frame)
  if grep -qE '^[[:space:]]+1[[:space:]]*$' grid.txt; then _pass "grid: step-number row is separate from names"; else _fail "grid: step-number row is separate from names"; fi

  chk() { local got; got="$(checkRunGrid --resolve "$1" --log smartSlurmLog 2>/dev/null)"
          [ "$got" = "$2" ] && _pass "resolve $1 -> $2" || _fail "resolve $1 -> $2" "got '$got'"; }
  chk a1 1.0.findNumber.1
  chk b1 1.0.findNumber.2
  chk e2 2.1.mergeNumber
  # reversed order must resolve to the same cell (verbs are separate tokens now)
  chk 1a 1.0.findNumber.1
  chk 1b 1.0.findNumber.2
  chk 2e 2.1.mergeNumber
  checkRunGrid --resolve z9 --log smartSlurmLog >/dev/null 2>&1; assert_rc $? 1 "resolve unknown address exits non-zero"
  # the boxed legend must keep its colored glyphs under --color (regression guard)
  if checkRunGrid --render --color --no-hint --log smartSlurmLog 2>&1 | grep -a 'multi-job' | grep -q "$(printf '\033')\["; then
    _pass "legend glyphs are colored with --color"
  else
    _fail "legend glyphs are colored with --color"
  fi
  echo
}

# ===========================================================================
_want linear && {
  echo "== checkRunGrid engine: linear view + integer-address agreement =="
  ss_new_workspace "$WORK/ui_linear"
  cp "$TEST_DIR"/fixtures/clean_2step/* .; chmod +x ./*.sh
  ss_run_pipeline "./pipeline.sh 777" "sbatch -p short -c 1 --mem 4G -t 50:0" >/dev/null 2>&1

  checkRunGrid --render --linear --no-color --no-hint --log smartSlurmLog > lin.txt 2>&1
  assert_contains lin.txt "linear view" "linear: linear view is produced"
  # regression: a successful render must exit 0 even with --no-hint (a trailing
  # '[ -z "$NO_HINT" ] && echo' once leaked exit 1, tripping checkRun's fallback)
  checkRunGrid --render --linear --no-hint --log smartSlurmLog >/dev/null 2>&1; assert_rc $? 0 "linear render exits 0 with --no-hint"
  checkRunGrid --render         --no-hint --log smartSlurmLog >/dev/null 2>&1; assert_rc $? 0 "grid render exits 0 with --no-hint"
  # the number shown for each job MUST resolve to that same job
  ok=1
  while read -r num _glyph flag; do
    [ -z "${num:-}" ] && continue
    n="$(checkRunGrid --resolve "$num" --log smartSlurmLog 2>/dev/null)"
    [ "$n" = "$flag" ] || { ok=0; echo "       mismatch: N=$num render=$flag resolve=$n"; }
  done < <(grep -E '^[[:space:]]+[0-9]+[[:space:]]' lin.txt | awk '{print $1, $2, $3}')
  [ "$ok" = 1 ] && _pass "linear: render number N == resolve N for every job" \
                || _fail "linear: render/resolve numbering disagree"
  echo
}

# ===========================================================================
_want linear_order && {
  echo "== checkRunGrid linear view: SUBMISSION order preserved (step 1 above step 3), not status-grouped =="
  ss_new_workspace "$WORK/ui_linorder"
  L=smartSlurmLog; mkdir -p "$L"
  # a ragged run (disjoint samples per step) -> linear view. step 1 is DONE,
  # step 3 is RUNNING: the old status-grouped order floated the running step-3
  # jobs to the top (numbered before the done step-1 jobs). They must now read in
  # submission order: 1.0.* , then 2.1.* , then 3.2.*.
  { printf 'job_id\tdepend_on\tjob_flag\n'; i=1000
    for n in 1 2 3; do printf '%s\tnull\t1.0.findNumber.%s\n' "$((i++))" "$n"; : > "$L/1.0.findNumber.$n.success"; done
    printf '%s\tnull\t2.1.mergeNumber\n' "$((i++))"; : > "$L/2.1.mergeNumber.success"
    for n in 4 5 6; do printf '%s\tnull\t3.2.findNumber1.%s\n' "$((i++))" "$n"; printf 'started\n' > "$L/3.2.findNumber1.$n.out"; done
  } > "$L/allJobs.txt"

  checkRunGrid --render --no-color --no-hint --log "$L" > lo.txt 2>&1
  assert_contains lo.txt "showing linear view" "ragged run renders the linear view"
  s1=$(grep -n '1\.0\.findNumber\.1'  lo.txt | head -1 | cut -d: -f1)
  s2=$(grep -n '2\.1\.mergeNumber'    lo.txt | head -1 | cut -d: -f1)
  s3=$(grep -n '3\.2\.findNumber1\.4' lo.txt | head -1 | cut -d: -f1)
  if [ -n "$s1" ] && [ -n "$s2" ] && [ -n "$s3" ] && [ "$s1" -lt "$s2" ] && [ "$s2" -lt "$s3" ]; then
    _pass "linear list is in submission order (step 1 < step 2 < step 3)"
  else _fail "linear list is in submission order" "lines s1=$s1 s2=$s2 s3=$s3"; fi

  # numbering follows submission order, and verb-forms + reversed address resolve
  chk() { local got; got="$(checkRunGrid --resolve "$1" --log "$L" 2>/dev/null)"
          [ "$got" = "$2" ] && _pass "resolve $1 -> $2" || _fail "resolve $1 -> $2" "got '$got'"; }
  chk 1 1.0.findNumber.1          # first submitted job is #1 (was a running step-3 job before)
  chk 4 2.1.mergeNumber
  chk 5 3.2.findNumber1.4
  chk 7 3.2.findNumber1.6
  echo
}

# ===========================================================================
_want linear_hint && {
  echo "== checkRun controller (pty): command hint matches the view (numbers in linear, b2 in grid) =="
  if command -v python3 >/dev/null 2>&1; then
    # (a) ragged -> linear view -> hint must instruct a job NUMBER, never 'b2'
    ss_new_workspace "$WORK/ui_linhint"
    : > .smartSlurm.log            # checkRun refuses to launch without the submission log
    L=smartSlurmLog; mkdir -p "$L"
    { printf 'job_id\tdepend_on\tjob_flag\n'; i=1000
      for n in 1 2 3; do printf '%s\tnull\t1.0.findNumber.%s\n' "$((i++))" "$n"; : > "$L/1.0.findNumber.$n.success"; done
      for n in 4 5 6; do printf '%s\tnull\t3.2.findNumber1.%s\n' "$((i++))" "$n"; printf 'started\n' > "$L/3.2.findNumber1.$n.out"; done
    } > "$L/allJobs.txt"
    export COLUMNS=100 LINES=44
    printf '%s\n' "@WAIT:0.6" "qq" | drive_checkrun lh.txt
    unset COLUMNS LINES
    assert_contains lh.txt "showing linear view" "ragged run renders the linear view in checkRun"
    assert_contains lh.txt "type a job number"   "linear hint tells the user to type a number"
    if grep -q 'b2' lh.txt; then _fail "linear hint must NOT mention the grid address b2"; else _pass "linear hint does not mention grid address b2"; fi

    # (b) clean grid -> grid view -> hint still uses the b2 grid address
    ss_new_workspace "$WORK/ui_gridhint"
    cp "$TEST_DIR"/fixtures/clean_2step/* .; chmod +x ./*.sh
    ss_run_pipeline "./pipeline.sh 777" "sbatch -p short -c 1 --mem 4G -t 50:0" >/dev/null 2>&1
    export COLUMNS=100 LINES=44
    printf '%s\n' "@WAIT:0.6" "qq" | drive_checkrun gh.txt
    unset COLUMNS LINES
    assert_contains gh.txt "b2" "grid hint uses the b2 grid address"
    if grep -q 'type a job number' gh.txt; then _fail "grid hint must NOT use the linear number wording"; else _pass "grid hint does not use the linear number wording"; fi
  else
    _fail "linear_hint: python3 required for the pty test (not found)"
  fi
  echo
}

# ===========================================================================
_want interactive && {
  echo "== checkRun controller (pty): full traversal — verbs, raw-log descent, list, redraw =="
  if command -v python3 >/dev/null 2>&1; then
    ss_new_workspace "$WORK/ui_interactive"
    cp "$TEST_DIR"/fixtures/clean_2step/* .; chmod +x ./*.sh
    MOCK_FAIL_JOB=mergeNumber MOCK_FAIL_MODE=oom \
      ss_run_pipeline "./pipeline.sh 777" "sbatch -p short -c 1 --mem 4G -t 50:0" >/dev/null 2>&1

    # Traversal script:
    #   d e2   diagnose the failed experiment step   (+ "" to ack the pause)
    #   a1     view a good job's log
    #   r a1   descend into the raw-file browser
    #   1      OPEN a raw file (pager) and return to the browser
    #   q      back up to the grid home
    #   l      list view  ->  q back to home
    #   qq     quit
    printf '%s\n' "d e2" "" "a1" "r a1" "1" "q" "l" "q" "qq" | drive_checkrun t.txt

    assert_contains t.txt "Showing run: /"          "home LOCATION BAR shows the current run dir as a full path"
    assert_contains t.txt "┌"                          "home shows the framed grid"
    assert_contains t.txt "diagnose"                  "home shows the command area (usage)"
    assert_contains t.txt "state="                     "'d e2' produced a diagnosis card"
    assert_contains t.txt "Raw files for 1.0.findNumber.1" "'r a1' descended into the raw-file browser"
    # opening a file re-renders the browser afterward: >=2 browser headers means
    # we selected a file, viewed it (pager), and returned to the browser
    nlist="$(grep -c 'Raw files for' t.txt)"
    [ "${nlist:-0}" -ge 2 ] && _pass "opened a raw file and returned to the browser" \
                            || _fail "opened a raw file and returned to the browser" "browser rendered ${nlist:-0}x"
    assert_contains t.txt "Available logs"             "'l' opened the list view"
    # after descending into raw logs and back, we must be at the grid home again
    if [ "$(grep -c 'Showing run: /' t.txt)" -ge 2 ]; then _pass "returned to grid home after descending into raw logs"; else _fail "returned to grid home after descending into raw logs"; fi

    home="$(grep -o '##esc_home=[0-9]*'  t.txt.meta | head -1 | cut -d= -f2)"
    [ "${home:-0}" -ge 3 ] && _pass "redraw is in-place (${home} cursor-home paints, no flooding)" \
                           || _fail "redraw is in-place" "esc_home=${home:-0}"
  else
    _fail "interactive: python3 required for the pty test (not found)"
  fi
  echo
}

# ===========================================================================
_want traversal && {
  echo "== checkRun lateral movement: switch between run directories via 'o' =="
  if command -v python3 >/dev/null 2>&1; then
    ss_new_workspace "$WORK/ui_traversal"
    cp "$TEST_DIR"/fixtures/clean_2step/* .; chmod +x ./*.sh
    ss_run_pipeline "./pipeline.sh 777" "sbatch -p short -c 1 --mem 4G -t 50:0" >/dev/null 2>&1
    # fabricate a SECOND, older run directory so there is somewhere to switch to
    cp -r smartSlurmLog smartSlurmLog.old
    touch -d '2 hours ago' smartSlurmLog.old smartSlurmLog.old/* 2>/dev/null

    # launch -> home should say 1 other run is available; open switcher; pick the
    # older run (#2 by recency); confirm the location bar now names it; quit.
    printf '%s\n' "@WAIT:0.5" "o" "2" "@WAIT:0.5" "qq" | drive_checkrun v.txt

    if grep -qE 'Showing run: .*/smartSlurmLog($|[^.])' v.txt; then _pass "launch shows the most-recent run as current"; else _fail "launch shows the most-recent run as current"; fi
    assert_contains v.txt "other run(s) available"            "home advertises that other runs exist"
    assert_contains v.txt "Runs and logs in"                  "'o' opened the run switcher"
    assert_contains v.txt "smartSlurmLog.old"                 "switcher lists the older run directory"
    if grep -q '/smartSlurmLog.old' v.txt; then _pass "after switching, location bar names the OTHER run dir"; else _fail "after switching, location bar names the OTHER run dir"; fi
  else
    _fail "traversal: python3 required for the pty test (not found)"
  fi
  echo
}

# ===========================================================================
_want waiting && {
  echo "== checkRun launch race: waiting screen auto-refreshes until run is ready =="
  if command -v python3 >/dev/null 2>&1; then
    ss_new_workspace "$WORK/ui_waiting"
    # Simulate 'checkRun run immediately after runAsPipeline': the submission log
    # exists, but the run's folder/allJobs.txt does not yet.
    : > .smartSlurm.log
    # a background task makes the run "ready" ~1s in, while checkRun waits
    ( sleep 1
      mkdir -p smartSlurmLog
      printf 'job_id\tdepend_on\tjob_flag\n101\tnull\t1.0.step.a\n' > smartSlurmLog/allJobs.txt
      : > smartSlurmLog/1.0.step.a.success
    ) &
    # short countdown; @WAIT lets the auto-refresh fire (no keypress), then quit
    printf '%s\n' "@WAIT:5" "qq" | CHECKRUN_WAIT_SECS=2 drive_checkrun w.txt

    assert_contains w.txt "Waiting for the run"     "shows the waiting screen when nothing is ready"
    assert_contains w.txt "refresh in"              "waiting screen shows the countdown"
    assert_contains w.txt "Showing run:"             "auto-entered the grid once the run appeared (no keypress)"
  else
    _fail "waiting: python3 required for the pty test (not found)"
  fi
  echo
}

# ===========================================================================
_want readycheck && {
  echo "== checkRun launch gate: a COMPLETED run must open the grid, never wait =="
  if command -v python3 >/dev/null 2>&1; then
    ss_new_workspace "$WORK/ui_ready"
    cp "$TEST_DIR"/fixtures/clean_2step/* .; chmod +x ./*.sh
    ss_run_pipeline "./pipeline.sh 777" "sbatch -p short -c 1 --mem 4G -t 50:0" >/dev/null 2>&1
    # regression for the cluster bug: readiness must NOT depend on the allJobs.txt
    # header text (a non-standard first column stranded a real completed run on
    # the waiting screen).
    sed -i '1s/^job_id/RUNID/' smartSlurmLog/allJobs.txt
    printf '%s\n' "@WAIT:0.6" "qq" | drive_checkrun ready.txt
    assert_contains    ready.txt "Showing run:"      "completed run opens the grid at launch"
    if grep -q "Waiting for the run" ready.txt; then _fail "completed run does NOT show the waiting screen"; else _pass "completed run does NOT show the waiting screen"; fi
  else
    _fail "readycheck: python3 required for the pty test (not found)"
  fi
  echo
}

# ===========================================================================
_want dryrun && {
  echo "== checkRun launch gate: a DRY-RUN folder browses scripts, never waits =="
  if command -v python3 >/dev/null 2>&1; then
    ss_new_workspace "$WORK/ui_dryrun"
    : > .smartSlurm.log
    mkdir -p smartSlurmLog
    printf '#!/bin/bash\necho dry\n' > smartSlurmLog/slurmPipeLine.deadbeef.dryrun.sh
    printf '%s\n' "@WAIT:0.6" "q" | drive_checkrun dry.txt
    assert_contains dry.txt "Dry-run log folder"                "dry-run folder shows the dry-run panel"
    assert_contains dry.txt "slurmPipeLine.deadbeef.dryrun.sh"  "dry-run panel lists the generated script"
    if grep -q "Waiting for the run" dry.txt; then _fail "dry-run folder does NOT show the waiting screen"; else _pass "dry-run folder does NOT show the waiting screen"; fi
  else
    _fail "dryrun: python3 required for the pty test (not found)"
  fi
  echo
}

# ===========================================================================
_want waitcap && {
  echo "== checkRun launch gate: bounded wait then give-up guidance (no run appears) =="
  if command -v python3 >/dev/null 2>&1; then
    ss_new_workspace "$WORK/ui_waitcap"
    : > .smartSlurm.log
    export CHECKRUN_WAIT_SECS=2 CHECKRUN_WAIT_MAX=3
    printf '%s\n' "@WAIT:6" | drive_checkrun cap.txt
    unset CHECKRUN_WAIT_SECS CHECKRUN_WAIT_MAX
    assert_contains cap.txt "Waiting for the run"        "shows the waiting screen while checking"
    assert_contains cap.txt "None of the expected files" "gives up with guidance after the wait cap"
  else
    _fail "waitcap: python3 required for the pty test (not found)"
  fi
  echo
}

# ===========================================================================
_want height && {
  echo "== checkRun height-adaptive: a tall frame collapses usage to one line =="
  if command -v python3 >/dev/null 2>&1; then
    ss_new_workspace "$WORK/ui_height"
    cp "$TEST_DIR"/fixtures/clean_2step/* .; chmod +x ./*.sh
    ss_run_pipeline "./pipeline.sh 777" "sbatch -p short -c 1 --mem 4G -t 50:0" >/dev/null 2>&1

    export COLUMNS=90 LINES=12
    printf '%s\n' "@WAIT:0.6" "qq" | drive_checkrun tall.txt
    unset LINES COLUMNS
    if grep -q 'b2 view  d b2 diagnose' tall.txt; then _pass "tall frame uses the compact one-line usage"; else _fail "tall frame uses the compact one-line usage"; fi
    if grep -q 'view log        (just type' tall.txt; then _fail "tall frame omits the verbose usage block"; else _pass "tall frame omits the verbose usage block"; fi

    export COLUMNS=90 LINES=44
    printf '%s\n' "@WAIT:0.6" "qq" | drive_checkrun norm.txt
    unset LINES COLUMNS
    assert_contains norm.txt "view log        (just type" "normal height shows the full usage block"
    assert_contains norm.txt "┌─ legend"                  "grid shows the boxed legend"
  else
    _fail "height: python3 required for the pty test (not found)"
  fi
  echo
}

# ===========================================================================
_want addressverb && {
  echo "== checkRun input grammar: a verb never collides with an address (the '3a' bug) =="
  if command -v python3 >/dev/null 2>&1; then
    ss_new_workspace "$WORK/ui_addrverb"
    cp "$TEST_DIR"/fixtures/clean_2step/* .; chmod +x ./*.sh
    ss_run_pipeline "./pipeline.sh 777" "sbatch -p short -c 1 --mem 4G -t 50:0" >/dev/null 2>&1
    export COLUMNS=90 LINES=44
    printf '%s\n' "@WAIT:0.5" "1a" "@WAIT:0.5" "qq" | drive_checkrun av.txt
    unset COLUMNS LINES
    if grep -q "Asking AI" av.txt; then _fail "'1a' must NOT trigger AI (no glued verb suffix)"; else _pass "'1a' does not silently trigger AI"; fi
    if grep -q "Did not understand '1a'" av.txt; then _fail "'1a' is accepted as an address (reversed order)"; else _pass "'1a' is accepted as an address (reversed order)"; fi
  else
    _fail "addressverb: python3 required for the pty test (not found)"
  fi
  echo
}

# ===========================================================================
_want experiment && {
  echo "== checkRunGrid: e-disambiguation (experiment 'E' vs row 'e') + capital-E display =="
  ss_new_workspace "$WORK/ui_exp"
  # hand-built 5-sample run: step 1 per-sample (rows a..e), step 2 experiment-level
  L=smartSlurmLog; mkdir -p "$L"
  { printf 'job_id\tdepend_on\tjob_flag\n'; i=1000
    for s in s1 s2 s3 s4 s5; do i=$((i+1)); printf '%s\tnull\t1.0.align.%s\n' "$i" "$s"; : > "$L/1.0.align.$s.success"; done
    printf '1099\tnull\t2.1.summary\n'; } > "$L/allJobs.txt"; : > "$L/2.1.summary.success"

  checkRunGrid --render --no-color --no-hint --log "$L" > exp.txt 2>&1
  assert_contains exp.txt "E2" "experiment step displays with a capital E"

  chkx() { local got; got="$(checkRunGrid --resolve "$1" --log "$L" 2>/dev/null)"
           [ "$got" = "$2" ] && _pass "resolve $1 -> $2" || _fail "resolve $1 -> $2" "got '$got'"; }
  chkx e1 1.0.align.s5     # row e (5th sample), step 1 per-sample -> reachable now
  chkx 1e 1.0.align.s5     # reversed
  chkx e2 2.1.summary      # experiment step 2 (lowercase accepted)
  chkx E2 2.1.summary      # capital accepted
  chkx 2e 2.1.summary      # reversed
  echo
}

# ===========================================================================
_want gapsteps && {
  echo "== checkRunGrid: grid uses ACTUAL #@ step numbers (gaps preserved, not renumbered) =="
  ss_new_workspace "$WORK/ui_gap"
  # a pipeline whose step numbers skip 2 and 3 (rAP allows increasing-with-gaps)
  L=smartSlurmLog; mkdir -p "$L"
  { printf 'job_id\tdepend_on\tjob_flag\n'; i=1000
    for s in alpha bravo; do
      i=$((i+1)); printf '%s\tnull\t1.0.align.%s\n'  "$i" "$s"; : > "$L/1.0.align.$s.success"
      i=$((i+1)); printf '%s\t0\t4.1.call.%s\n'      "$i" "$s"; : > "$L/4.1.call.$s.success"
      i=$((i+1)); printf '%s\t0\t5.4.report.%s\n'    "$i" "$s"; : > "$L/5.4.report.$s.success"
    done; } > "$L/allJobs.txt"

  checkRunGrid --render --no-color --no-hint --log "$L" > gap.txt 2>&1
  # the step-number header row must read "1 4 5", never "1 2 3"
  if grep -qE '^[[:space:]]+1[[:space:]]+4[[:space:]]+5[[:space:]]*$' gap.txt; then _pass "column header shows real step numbers 1 4 5"; else _fail "column header shows real step numbers 1 4 5"; fi
  if grep -qE '^[[:space:]]+1[[:space:]]+2[[:space:]]+3[[:space:]]*$' gap.txt; then _fail "columns must NOT be renumbered to 1 2 3"; else _pass "columns are not renumbered to 1 2 3"; fi

  chkg() { local got; got="$(checkRunGrid --resolve "$1" --log "$L" 2>/dev/null)"
           [ "$got" = "$2" ] && _pass "resolve $1 -> $2" || _fail "resolve $1 -> $2" "got '$got'"; }
  chkg a4 4.1.call.alpha       # addressed by the ACTUAL step number
  chkg b5 5.4.report.bravo
  checkRunGrid --resolve a2 --log "$L" >/dev/null 2>&1; assert_rc $? 1 "a2 (skipped step) resolves to nothing"
  echo
}

# ===========================================================================
_want listview && {
  echo "== checkRun list view: consecutive numbering + untruncated summaries =="
  if command -v python3 >/dev/null 2>&1; then
    ss_new_workspace "$WORK/ui_list"
    : > .smartSlurm.log
    L=smartSlurmLog; mkdir -p "$L"
    { printf 'job_id\tdepend_on\tjob_flag\n1\tnull\t1.0.align.s1\n2\tnull\t1.0.align.s2\n3\tnull\t2.1.count.s1\n'; } > "$L/allJobs.txt"
    : > "$L/1.0.align.s1.success"                 # job 3 has no flags -> UNCONFIRMED (summary contains '/')
    export COLUMNS=120 LINES=44
    printf '%s\n' "@WAIT:0.6" "l" "@WAIT:0.4" "q" "qq" | drive_checkrun lv.txt
    unset COLUMNS LINES
    assert_contains lv.txt "check the cluster / requeue" "list view shows full summaries (no '/' truncation)"
    if grep -qE '^1 ' lv.txt && grep -qE '^2 ' lv.txt && grep -qE '^3 ' lv.txt; then _pass "list view numbers rows consecutively (1,2,3)"; else _fail "list view numbers rows consecutively"; fi
  else
    _fail "listview: python3 required for the pty test (not found)"
  fi
  echo
}

# ===========================================================================
_want multidep && {
  echo "== checkRunGrid: multi-dependency flags parse to the right sample (grid, not linear) =="
  ss_new_workspace "$WORK/ui_multidep"
  L=smartSlurmLog; mkdir -p "$L"
  # flags with a VARIABLE number of numeric dependency fields between step & name
  { printf 'job_id\tdepend_on\tjob_flag\n'; i=1000
    for s in sampA sampB; do
      i=$((i+1)); printf '%s\t0\t2.0.trim.%s.fastq.gz\n'       "$i" "$s"   # 1 dep
      i=$((i+1)); printf '%s\t0\t3.1.2.map.%s.fastq.gz\n'      "$i" "$s"   # 2 deps
      i=$((i+1)); printf '%s\t0\t7.4.5.6.dedup.%s.fastq.gz\n'  "$i" "$s"   # 3 deps
    done
    printf '1099\t0\t11.10.multistats\n'; } > "$L/allJobs.txt"
  checkRunGrid --render --no-color --no-hint --log "$L" > md.txt 2>&1
  if grep -q "┌" md.txt && ! grep -q "does not flatten" md.txt; then _pass "multi-dep flags still flatten to a grid (not linear)"; else _fail "multi-dep flags still flatten to a grid (not linear)"; fi
  chkm() { local got; got="$(checkRunGrid --resolve "$1" --log "$L" 2>/dev/null)"
           [ "$got" = "$2" ] && _pass "resolve $1 -> $2" || _fail "resolve $1 -> $2" "got '$got'"; }
  chkm a3 3.1.2.map.sampA.fastq.gz       # step 3 (2 deps) -> sampA
  chkm b7 7.4.5.6.dedup.sampB.fastq.gz   # step 7 (3 deps) -> sampB
  chkm E11 11.10.multistats              # experiment step (no sample)
  echo
}

# ===========================================================================
_want sampletrim && {
  echo "== checkRunGrid: sample labels are trimmed to DISTINCT strings (prefix/suffix+window) =="
  ss_new_workspace "$WORK/ui_trim"
  L=smartSlurmLog; mkdir -p "$L"
  # three names sharing a long prefix AND suffix, diverging deep in the middle
  { printf 'job_id\tdepend_on\tjob_flag\n'; i=1000
    for s in PFX_LIB066167_SRN00305666_SFX PFX_LIB066167_SRN00305667_SFX PFX_LIB069225_SRN00322925_SFX; do
      i=$((i+1)); printf '%s\t0\t2.0.trim.%s\n' "$i" "$s"
    done; } > "$L/allJobs.txt"
  export CHECKRUN_SAMPLE_W=8
  local nrows ndist
  nrows=$(checkRunGrid --samples --log "$L" | grep -c .)
  ndist=$(checkRunGrid --samples --log "$L" | cut -f2 | sort -u | grep -c .)
  unset CHECKRUN_SAMPLE_W
  [ "$nrows" -eq 3 ] && _pass "--samples emits one row per sample (a,b,c)" || _fail "--samples row count" "got $nrows"
  [ "$ndist" -eq 3 ] && _pass "hierarchically-similar names get DISTINCT trimmed labels at W=8" || _fail "distinct labels" "got $ndist distinct at W=8"
  # the full name column carries the untrimmed sample back
  if checkRunGrid --samples --log "$L" | grep -q 'PFX_LIB069225_SRN00322925_SFX'; then _pass "--samples carries the full sample name"; else _fail "--samples carries the full sample name"; fi
  echo
}

echo "==========================================================="
echo "checks: $CHECKS   failures: $FAILS"
[ "$FAILS" -eq 0 ] && { echo "RESULT: PASS"; exit 0; } || { echo "RESULT: FAIL"; exit 1; }