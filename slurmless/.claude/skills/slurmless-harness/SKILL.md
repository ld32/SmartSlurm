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
(`SMARTSLURM_ROOT` auto-resolves inside SmartSlurm at `SmartSlurm/slurmless/`;
to test a different checkout, `export SMARTSLURM_ROOT=/path/to/SmartSlurm`.)

## Common tasks

**Run the suite / a block**
```bash
./test/run_tests.sh --list        # catalog of blocks + what each verifies
./test/run_tests.sh <block>       # one block (substring match); run one at a time — the suite is slow
./test/checkrun_ui_tests.sh       # interactive-monitor UI suite (needs python3)
```

**Reproduce a failure mode**
```bash
MOCK_FAIL_MODE=oom|timeout|nonzero  MOCK_FAIL_JOB=<jobName-substring>  ./test/run_tests.sh oom_injection
```

**Add a fixture + block**
1. `test/fixtures/<name>/pipeline.sh` (annotated `#@` script) + any bare-name step tools.
2. A block in `test/run_tests.sh`: `_begin <name>`, `ss_new_workspace`, copy fixture,
   `ss_run_pipeline "./pipeline.sh …" "sbatch …"`, assert with the helpers.
3. Add `<name>` to `ALL_BLOCKS` and an `ss_purpose` line. Keep payloads tiny and deterministic.

**Exercise the resource-estimation learning loop** (needs `gnuplot`): use
`ss_use_model` / `ss_run_sweep` / `ss_fit` / `ss_seed_record` / `assert_near`
(see `test/README.md` → "M1").

**Check sacct-format drift** (after a Slurm upgrade; run on the cluster):
```bash
./test/mockSlurm/check_sacct_format.sh <completedJobId>
```

## Guardrails

- Do **not** edit SmartSlurm `bin/` to force a green test. A failing test that
  reveals a real tool bug is a finding — report it. Change a tool only when the
  task is explicitly to change that tool, then update the affected fixtures/asserts.
- Scratch only (`$SS_SCRATCH`); never touch a real `~/.smartSlurm` or cluster.
- On a shared host, opt out of env mutations: `SS_NO_INSTALL=1`, `SS_NO_SH_REPOINT=1`.
- Develop slurmless in-tree: edit `slurmless/…` and open a pull request to
  SmartSlurm, like any other part of the repo (no separate slurmless repository).
