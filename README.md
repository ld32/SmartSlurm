# SmartSlurm

**SmartSlurm** automatically estimates and optimizes memory and run-time for Slurm jobs to increase resource efficiency and decrease the amount of time jobs spend in pending status. It is pure Bash — no daemon, no database, no special runtime. Setup is to clone the repo and put it on your `PATH`.

It has two parts you can use together or separately:

| Tool | What it is | Use it when |
|------|-----------|-------------|
| **`ssbatch`** | A drop-in `sbatch` wrapper that estimates memory/time from your past jobs, picks a partition, resubmits on OOM/OOT, and sends informative emails. | You submit individual jobs and want smart resource sizing. |
| **`runAsPipeline`** | A workflow runner built on `ssbatch`. You annotate an ordinary Bash script with `#@` markers; it wires up dependencies and submits each step as its own smart job. | You have a multi-step pipeline with dependencies. |

Two companion tools help you watch and debug runs:

| Tool | What it is |
|------|-----------|
| **`checkRun`** | An interactive monitor: a status **grid** of every step, with per-job log viewing, diagnosis, and reruns. |
| **`diagnose.sh`** | A one-shot "what happened / why / what to do" report for a single job — used inside `checkRun`, or on its own. |

> [!NOTE]
> Because it is plain Bash, installation is just `git clone`, and it slots into existing command-line tools and pipelines (Snakemake, Nextflow, Cromwell) without rewriting them.

---


## Features

- **Auto-sizes memory and run-time** from statistics of your earlier jobs.
- **Auto-selects the partition** that matches the requested run-time.
- **Auto-resubmits** jobs that die out-of-memory (OOM) or out-of-time (OOT), with doubled resources.
- **Optional checkpointing**: snapshot a long job before it hits its limit, then resume from the snapshot.
- **Consolidated emails**: instead of one email per job, SmartSlurm sends one per-sample alert (failed / completed) and one per-run efficiency summary — text-only, no inbox flood.
- **Conda-style activation**: a one-time `bash smartSlurmSetup.sh` adds a `smartslurm` command to your shell (it edits `~/.bashrc` for you — no manual `export`); `smartslurm activate` turns SmartSlurm on (aliasing `sbatch` → smart `ssbatch`, with a `{smartSlurm} ` prompt indicator) and `smartslurm deactivate` cleanly reverses it.
- **Resource efficiency report**: a per-run `efficiency.tsv` and an `htop`-style view in `checkRun` (press `e`) showing time/memory used vs allocated, per sample and per step.
- **(runAsPipeline) Dependency management**: steps wait for their prerequisites automatically.
- **(runAsPipeline) Explicit `#@ … #@end` blocks**: write a step's body verbatim to a per-job script (`awk`, quotes, heredocs, multi-line all survive), instead of flattening it into a single `--wrap` string.
- **(runAsPipeline) `--lint`**: validate a pipeline's block structure without submitting anything.
- **(runAsPipeline) Smart reruns**: an unchanged script is reused as-is; on re-issue you pick, from one menu, whether to redo everything, only failures, or from a chosen step.
- **(checkRun) Status grid + diagnosis**: a sample × step grid, addressable per cell, with inline "why did this fail" powered by `diagnose.sh`.

---


## Quick start

SmartSlurm is pure Bash — installation is a `git clone` plus a one-time setup
script that adds the `smartslurm` command to your shell:

```bash
bash smartSlurmSetup.sh      # one-time; edits ~/.bashrc for you (no manual export)
smartslurm activate          # turn SmartSlurm on (in a new shell); aliases sbatch -> ssbatch
```

Then submit a single smart job with `ssbatch`, or run an annotated multi-step
pipeline with `runAsPipeline`.

**&rarr; Install details:** [docs/installation.md](docs/installation.md)
&middot; **Full walkthrough:** [docs/quick-start.md](docs/quick-start.md)

---

## Documentation

Full guides live in [`docs/`](docs/):

**Getting started**
- [Installation & setup](docs/installation.md) — clone, the `smartslurm` command, activate/deactivate, optional Graphviz.
- [Quick start](docs/quick-start.md) — a first `ssbatch` job and a first `runAsPipeline` pipeline, end to end.

**Using SmartSlurm**
- [Core concepts](docs/core-concepts.md) — how a pipeline runs (convert → submit), the `#@` job block, scopes, loops, resource estimation, jobRecords, notifications.
- [Tools reference](docs/tools.md) — `runAsPipeline`, `ssbatch`, `reviewJobRecords.py`, `cancelAllJobs`, `sbatchAndTop`, and the under-the-hood helpers.
- [Monitoring & diagnosing runs](docs/monitoring.md) — the `checkRun` status grid, `diagnose.sh`, the resource-efficiency view, and rerunning.
- [Configuration](docs/configuration.md) — `config.txt`, key settings, personal overrides.

**Going further**
- [Using ssbatch with other pipeline managers](docs/pipeline-managers.md) — Snakemake, Nextflow, Cromwell.
- [Testing off-cluster (slurmless)](docs/testing-off-cluster.md) — the Slurm-free sandbox, vendored at `slurmless/`.
- [FAQ & troubleshooting](docs/faq.md) — common `ssbatch` and `runAsPipeline` questions.
- [Upgrading](docs/upgrading.md) · [Contributing](docs/contributing.md)

See [`docs/`](docs/README.md) for the full index.
