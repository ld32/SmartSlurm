# AGENTS.md — working in slurmless

Model-neutral session-startup guide for any coding agent (Anthropic, OpenAI,
Cursor, …) operating in this directory. Humans: this doubles as a quick operator
guide; see `CONTEXT.md` and `test/README.md` for depth.

## What this is

**slurmless** is a **Slurm-free test harness** for
[SmartSlurm](https://github.com/ld32/smartSlurm). It runs SmartSlurm's
**unmodified** tools (`runAsPipeline`, `ssbatch`, `cleanUp.sh`, `checkRun`,
`checkRunGrid`, `diagnose.sh`) with **no Slurm scheduler**, by interposing mock
Slurm commands (`sbatch`, `sacct`, `srun`, …) on `$PATH`. The mocks actually run
each job's payload, so real `.success`/`.failed` flags and logs are produced —
failure attribution behaves as on a cluster.

It is a **component of SmartSlurm**, living at `SmartSlurm/slurmless/` as ordinary
files and developed in-tree — edit `slurmless/…` and open a pull request to
SmartSlurm, like any other part of the repo.

## Session startup (do this first)

1. **Confirm the SmartSlurm checkout.** In-tree (`SmartSlurm/slurmless/`):
   `SMARTSLURM_ROOT` auto-resolves. Testing a different SmartSlurm checkout:
   `export SMARTSLURM_ROOT=/path/to/SmartSlurm`.
2. **Do not modify SmartSlurm's `bin/` to make a test pass.** The harness tests
   the *unmodified* tools. If a test reveals a tool bug, report it; change a tool
   only when the task explicitly is to change that tool (then update the affected
   fixtures/asserts).
3. **Prove the baseline before anything else:** `./bootstrap.sh` (runs the
   `clean_2step` block). If that is red, stop and diagnose — the environment is
   not set up.
4. **Run one block at a time.** The mock scheduler sleeps through SmartSlurm's
   cleanup counters, so the full suite is slow; a single block is usually quick.
   Curve-fit blocks (`scaling_*`, `noise_fit`, `contaminate_*`) need `gnuplot`.

## Golden-path commands

```bash
./bootstrap.sh                         # one-command: resolve SmartSlurm + run baseline
./test/run_tests.sh --list             # every block and what it verifies
./test/run_tests.sh <block>            # run one block (substring match)
./test/checkrun_ui_tests.sh            # the interactive-monitor UI suite (needs python3)
./test/mockSlurm/check_sacct_format.sh # sacct-format drift check (CI-friendly)
```

Interactive:

```bash
source test/lib/harness.sh && ss_activate
ss_new_workspace /tmp/play
cp test/fixtures/clean_2step/* .
ss_run_pipeline "./pipeline.sh 777" "sbatch -p short -c 1 --mem 4G -t 50:0"
```

## What `ss_activate` changes (and how to opt out)

It is built for a disposable sandbox / CI runner. It will: install `gnuplot` if
missing, repoint `/bin/sh → bash` (SmartSlurm's `#!/bin/sh` scripts need bash;
production RHEL `/bin/sh` *is* bash), stage the mocks on `$PATH`, and redirect
`HOME` into scratch (so nothing touches a real `~/.smartSlurm`). On a shared
host, opt out: `SS_NO_INSTALL=1`, `SS_NO_SH_REPOINT=1`. It prints a one-line
summary of what it touched.

## Adding a test

1. `mkdir test/fixtures/<name>/` with a `pipeline.sh` (annotated `#@` script) and
   any step tools it calls by bare name.
2. Add a block to `test/run_tests.sh`: `_begin <name>`, `ss_new_workspace`, copy
   the fixture, `ss_run_pipeline`, assert with the helpers (`assert_file`,
   `assert_nonempty`, `assert_contains`, `assert_valid_bash`, `assert_rc`, …).
3. Add `<name>` to `ALL_BLOCKS` and give it an `ss_purpose` line. Keep payloads
   tiny and deterministic.

## Simulating failures

```bash
MOCK_FAIL_MODE=oom|timeout|nonzero  MOCK_FAIL_JOB=<jobName-substring> ss_run_pipeline ...
```
Injected jobs are killed before `.success` is touched, so `cleanUp.sh` falls
through to the mock `sacct` state exactly as on a real cluster.

## Drift — the thing to protect

- **Tools vs harness:** in-tree, they share one commit (can't drift); run the
  harness before pushing a tool change (`./bootstrap.sh`, then the relevant
  blocks). `test/version.txt` records the validated baseline for runs against an
  external checkout.
- **sacct format** is pinned to production Slurm (25.11.7); it is the most
  version-sensitive file. If a Slurm upgrade changes `sacct` output, run
  `test/mockSlurm/check_sacct_format.sh` on the cluster and update
  `test/mockSlurm/sacct`.
- **Trap-line coupling:** the resource model reads the generated `cleanUp.sh`
  trap line positionally; the mock `sbatch` now warns loudly if that shape
  changes. Heed the warning rather than trusting a modeled result.

## Hard rules

- Never edit SmartSlurm `bin/` to force a green test (see startup #2).
- Never point the harness at a real `~/.smartSlurm` or a real cluster; scratch
  only (`$SS_SCRATCH`, default `/tmp/ss-scratch`).
- Develop slurmless in-tree: edit `slurmless/…` and open a pull request to
  SmartSlurm, like any other part of the repo.
