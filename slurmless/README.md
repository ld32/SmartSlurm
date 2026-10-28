# slurmless

**A Slurm-free test harness for [SmartSlurm](https://github.com/ld32/smartSlurm).**

slurmless lets you run SmartSlurm's **unmodified** tools — `runAsPipeline`,
`ssbatch`, `cleanUp.sh`, `checkRun`, `checkRunGrid`, `diagnose.sh` — on a machine
with **no Slurm scheduler** (a laptop, a CI runner, or an AI-agent sandbox), while
still executing real job payloads. That makes it possible to develop and
regression-test SmartSlurm — pipeline conversion, dependency wiring, failure
attribution, the status grid, the rerun menu, and the resource-estimation
learning loop — entirely off-cluster.

It works by interposing **mock Slurm commands** (`sbatch`, `sacct`, `srun`, …) on
`$PATH` ahead of the real binaries. SmartSlurm is **not modified**; the mocks
actually execute each job's payload, so the real `cleanUp.sh` trap and the real
`.success`/`.failed` flags are produced, exactly as on a cluster.

## Who it's for

- **SmartSlurm developers** — a fast, Slurm-free inner loop and regression net for
  work on `runAsPipeline`, resource estimation, logging, and the monitor UI.
- **Pipeline builders (bioinformaticians)** — a safe place to convert a pipeline
  to `runAsPipeline` and watch how it behaves before running it on real HPC.
- **AI agents** — a self-contained environment to develop and test SmartSlurm
  changes without a scheduler. (Agent onboarding: see `AGENTS.md` and the
  `/slurmless-harness` skill.)

## Quick start

slurmless is a component of SmartSlurm, living at `SmartSlurm/slurmless/`, so from
a SmartSlurm checkout:

```bash
cd SmartSlurm
./slurmless/test/run_tests.sh --list          # see every test block and what it checks
./slurmless/test/run_tests.sh clean_2step     # run the baseline block (must pass)
```

To run it against a different SmartSlurm checkout (e.g. a separate clone):

```bash
export SMARTSLURM_ROOT=/path/to/SmartSlurm
/path/to/slurmless/test/run_tests.sh clean_2step
```

A single block is usually quick; the full suite is slow (the mock scheduler
sleeps through SmartSlurm's cleanup counters), so pass a block name in a
time-boxed environment.

## Layout

| Path | What it is |
|------|------------|
| `test/run_tests.sh` | The pipeline + resource-model test suite (CI exit code). |
| `test/checkrun_ui_tests.sh` | The interactive-monitor (`checkRun`/`checkRunGrid`) UI suite, driven over a pseudo-terminal. |
| `test/lib/harness.sh` | Environment activation, the resource-model API, and assertion helpers. |
| `test/mockSlurm/` | The mock Slurm commands (`sbatch`, `sacct`, `srun`, `squeue`, …). |
| `test/fixtures/` | Tiny annotated pipelines exercising specific behaviors. |
| `test/models/` | Resource-curve models for the estimator learning-loop tests. |
| `test/version.txt` | The SmartSlurm baseline this harness was last validated against. |

## Learn more

- **[CONTEXT.md](CONTEXT.md)** — setup, how it stays in sync with SmartSlurm, the
  harness API, failure injection, and gotchas.
- **[test/README.md](test/README.md)** — the full block catalog, the resource-model
  layer, fidelity tiers, and how to add a fixture.
- **[test/mockSlurm/README.md](test/mockSlurm/README.md)** — the per-command mock
  contracts.
