---
name: slurmless-harness
description: Run or develop with the slurmless Slurm-free SmartSlurm test harness. Use when testing or developing SmartSlurm (runAsPipeline, ssbatch, cleanUp.sh, checkRun, checkRunGrid, diagnose.sh, resource estimation) without a Slurm scheduler, running the mock-Slurm test suite, adding a test fixture, reproducing an OOM/timeout/failure, or converting a pipeline to runAsPipeline and testing it off-cluster.
---

# slurmless harness

slurmless runs SmartSlurm's **unmodified** tools with **no Slurm**, by interposing
mock Slurm commands on `$PATH`. The mocks execute each job's real payload, so
`.success`/`.failed` flags, logs, and failure attribution behave as on a cluster.
Full operator detail is in `AGENTS.md`, `CONTEXT.md`, and `test/README.md`.

## First move, always

```bash
./bootstrap.sh            # resolve the SmartSlurm checkout + run the clean_2step baseline
```
Green baseline → proceed. Red → the environment isn't set up; diagnose before anything else.
`SMARTSLURM_ROOT` auto-resolves in-tree (`SmartSlurm/slurmless/`); for a separate
checkout, `export SMARTSLURM_ROOT=/path/to/SmartSlurm`.

## Running the suite

```bash
./test/run_tests.sh               # WHOLE suite — PARALLEL by default (one worker per core, capped at 8)
./test/run_tests.sh <block> …     # named block(s), substring match — serial inline, live output
./test/run_tests.sh -j1           # whole suite, forced serial (debugging, or a low-core host)
./test/run_tests.sh -jN …         # N workers explicitly
./test/run_tests.sh --list        # catalog: every block + the one line it verifies
```
- **One** selected block (and `-j1`) runs on the **serial inline path** and streams
  output live; **two or more** fan out to one worker each. A worker runs **exactly its
  one block** — names stay isolated even when one is a substring of another (e.g.
  `rerun` vs `rerun_lifecycle`); the `parallel_worker_isolation` block guards this.
- Scratch defaults to **RAM** (`/dev/shm`, with a disk fallback when tmpfs is absent or
  too small). The activation line prints the path and kind it chose; override with
  `SS_SCRATCH=<path>`.
- Opt-in scaffolds (`config_*`, `p6_*`, `cleanup_reentrancy`) are **not** in the
  default suite — they run only when named explicitly.

## Common tasks

**Reproduce a failure mode**
```bash
MOCK_FAIL_MODE=oom|timeout|nonzero  MOCK_FAIL_JOB=<jobName-substring>  ./test/run_tests.sh oom_injection
```

**Add a fixture + block**
1. `test/fixtures/<name>/pipeline.sh` (annotated `#@` script) + any bare-name step tools.
2. A block in `test/run_tests.sh`: `_begin <name>`, `ss_new_workspace`, copy the fixture,
   `ss_run_pipeline "./pipeline.sh …" "sbatch …"`, assert with the helpers.
3. Add `<name>` to `ALL_BLOCKS` and an `ss_purpose` line. Keep payloads tiny and
   deterministic. **Changing the runner itself** (dispatch, `ss_parallel_run`, scratch,
   `ss_activate`)? It needs its own self-test and a bounded process-count check — see
   the `AGENTS.md` hard rules; "the suite passed" can't certify the thing that runs it.

**Exercise the resource-estimation learning loop** (needs `gnuplot`): use
`ss_use_model` / `ss_run_sweep` / `ss_fit` / `ss_seed_record` / `assert_near`
(see `test/README.md` → "M1").

**Check sacct-format drift** (after a Slurm upgrade; run on the cluster):
```bash
./test/mockSlurm/check_sacct_format.sh <completedJobId>
```

## Guardrails

- Do **not** edit SmartSlurm `bin/` to force a green test. A failing test that reveals a
  real tool bug is a finding — report it. Change a tool only when the task explicitly is
  to change that tool, then update the affected fixtures/asserts.
- Scratch only. Never touch a real `~/.smartSlurm` or a real cluster; the harness
  confines everything to `$SS_SCRATCH`.
- On a shared host, opt out of env mutations with `SS_NO_INSTALL=1` and
  `SS_NO_SH_REPOINT=1` (these also skip per-worker setup; `ss_parallel_run` sets them
  for its workers).
- Source of truth is **`SmartSlurm/slurmless/` (in-tree)** — develop there directly,
  versioned with the engine. The former standalone repo is archived; there is no
  re-snapshot step.
