# slurmless

**A Slurm-free test harness for [SmartSlurm](https://github.com/ld32/smartSlurm).**

slurmless runs SmartSlurm's **unmodified** tools — `runAsPipeline`, `ssbatch`,
`cleanUp.sh`, `checkRun`, `checkRunGrid`, `diagnose.sh` — on a machine with **no
Slurm scheduler** (a laptop, a CI runner, or an AI-agent sandbox), while still
executing real job payloads. It interposes **mock Slurm commands** (`sbatch`,
`sacct`, `srun`, …) on `$PATH`; SmartSlurm itself is **not modified**, and the mocks
actually execute each job's payload, so the real `cleanUp.sh` trap and genuine
`.success`/`.failed` flags are produced, exactly as on a cluster.

That makes it possible to develop and regression-test SmartSlurm **off-cluster**,
end to end:

- pipeline conversion, dependency wiring, and failure attribution (`runAsPipeline`)
- the resource-estimation learning loop — curve fitting, extrapolation, and
  auto / manual resubmit
- the `checkRun` / `checkRunGrid` status grid, `diagnose.sh`, and the rerun menu
- resource-efficiency reporting and consolidated email
- the `smartslurm` activation command (activate / deactivate / status)

The current suite is ~40 test blocks spanning all of the above.

## Who it's for

- **SmartSlurm developers** — a fast, Slurm-free inner loop and regression net.
- **Pipeline builders (bioinformaticians)** — a safe place to convert a pipeline to
  `runAsPipeline` and watch how it behaves before running it on real HPC.
- **AI agents** — a self-contained environment to develop and test SmartSlurm
  changes without a scheduler. (Agent onboarding: see `AGENTS.md` and the
  `/slurmless-harness` skill.)

## Quick start

From a SmartSlurm checkout, one command sets things up and proves the baseline:

```bash
./slurmless/bootstrap.sh                       # resolve the checkout, activate, run clean_2step
```

Then run whatever you need:

```bash
./slurmless/test/run_tests.sh --list           # every block and what it verifies
./slurmless/test/run_tests.sh clean_2step       # one block (the baseline; must pass)
./slurmless/test/checkrun_ui_tests.sh           # the checkRun / checkRunGrid UI suite
```

A single block is usually quick; the full suite is slow (the mock scheduler sleeps
through SmartSlurm's cleanup counters), so pass a block name in a time-boxed
environment. To run against a different SmartSlurm checkout, point it there:

```bash
export SMARTSLURM_ROOT=/path/to/SmartSlurm
/path/to/slurmless/test/run_tests.sh clean_2step
```

> **Heads-up:** activation mutates the environment for faithful emulation — it
> repoints `/bin/sh`→bash, installs `gnuplot`, and redirects `HOME` into scratch.
> On a shared host, opt out with `SS_NO_SH_REPOINT=1` / `SS_NO_INSTALL=1`. Nothing
> ever touches a real `~/.smartSlurm`. See [CONTEXT.md](CONTEXT.md) for details.

## Layout

| Path | What it is |
|------|------------|
| `bootstrap.sh` | One command: resolve `SMARTSLURM_ROOT`, activate, run the baseline block. |
| `test/run_tests.sh` | The main test suite (~40 blocks; non-zero exit on any failure). |
| `test/checkrun_ui_tests.sh` | The `checkRun` / `checkRunGrid` UI suite, driven over a pseudo-terminal. |
| `test/lib/harness.sh` | Environment activation, the resource-model API, and assertion helpers. |
| `test/mockSlurm/` | The mock Slurm commands (`sbatch`, `sacct`, `srun`, …) + a `sacct`-format drift check. |
| `test/fixtures/` | Tiny annotated pipelines exercising specific behaviors. |
| `test/models/` | Resource-curve models for the estimator learning-loop tests. |
| `test/version.txt` | The SmartSlurm baseline this harness was last validated against. |
| `CONTEXT.md` | Setup, the harness API, failure injection, and gotchas. |
| `AGENTS.md` / `CLAUDE.md` | Model-neutral and Claude-specific onboarding. |
| `.claude/skills/` | The `/slurmless-harness` skill stub (plain Markdown; any agent can read it). |

## Learn more

- **[CONTEXT.md](CONTEXT.md)** — setup, the harness API, failure injection, and gotchas.
- **[test/README.md](test/README.md)** — the full block catalog, the resource-model
  layer, fidelity tiers, and how to add a fixture.
- **[test/mockSlurm/README.md](test/mockSlurm/README.md)** — the per-command mock contracts.
- **[AGENTS.md](AGENTS.md)** — the session-startup guide for coding agents.
