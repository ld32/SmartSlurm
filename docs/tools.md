[&larr; SmartSlurm README](../README.md) &middot; [Docs index](README.md)

# Tools

## Tools at a glance

| Tool | Kind | What it does | Details |
|------|------|--------------|---------|
| `runAsPipeline` | command | Convert and run an annotated Bash script as a dependency-aware pipeline. | [runAsPipeline](#runaspipeline) |
| `ssbatch` | command | Drop-in `sbatch` wrapper with smart sizing, partition selection, and OOM/OOT resubmit. | [ssbatch](#ssbatch) |
| `checkRun` | command | Interactive status grid + log/diagnosis/efficiency browser for a run. | [The status grid](monitoring.md#the-status-grid-checkrun) |
| `checkRunGrid` | command | Non-interactive render/resolve engine behind `checkRun` (scriptable). | [The status grid](monitoring.md#the-status-grid-checkrun) |
| `diagnose.sh` | command | "What happened / why / what to do" report for one job (human or `--json`). | [Diagnosing a job](monitoring.md#diagnosing-a-job-diagnosesh) |
| `reviewJobRecords.py` | command | Terminal tool to prune outlier records so estimates stay accurate. | [reviewJobRecords.py](#reviewjobrecordspy) |
| `cancelAllJobs` | command | Cancel running/pending jobs submitted from the current directory. | [cancelAllJobs](#cancelalljobs) |
| `sbatchAndTop` | command | Submit with `ssbatch`, then `scontrol top` the job. | [sbatchAndTop](#sbatchandtop) |
| `cleanUp.sh` | internal | Per-job EXIT-trap teardown: writes the jobRecord row, prunes/invalidates stats, decides OOM/OOT requeue, sends emails. | [Under-the-hood tools](#under-the-hood-tools) |
| `estimateResource.sh` | internal | Computes the memory/time estimate for a (pipeline, step, ref) group. | [Under-the-hood tools](#under-the-hood-tools) |
| `jobRecordSchema.sh` | internal | Single source of truth for the jobRecord schema; sourced by writers/readers. | [Under-the-hood tools](#under-the-hood-tools) |
| `jrRepair.sh` | command | One-time migration of a legacy (v1) jobRecord store to the current (v2) schema. | [Under-the-hood tools](#under-the-hood-tools) |

## runAsPipeline

`runAsPipeline` turns an ordinary Bash script into a dependency-aware pipeline. You annotate steps with [`#@` markers](core-concepts.md#the--job-block); each marked step is submitted through `ssbatch`, so it inherits smart sizing, resubmission, and emails.

**What it gives you**

- **Dependency management** — steps wait for the steps they depend on.
- **Per-step smart sizing** — every step is an `ssbatch` job.
- **Smart reruns** — an unchanged script is reused without reconversion; on re-issue you choose what to redo ([Rerunning a pipeline](monitoring.md#rerunning-a-pipeline)).
- **Failure containment** — if a step fails, its downstream steps are not run.
- **Multi-account support** — add `-A`/`--account=` and all jobs use that Slurm account.

**Usage**

```
runAsPipeline --script "SCRIPT [ARGS]" --tmp {useTmp|noTmp} [--pipeline NAME]
              [--sbatch-options "SBATCH_OPTIONS"] [--mode dryrun] [--lint]
              [--email noEmail|noSuccEmail] [--special checkpoint|excludeFailedNodes]
              [-v|--verbose | -q|--quiet]
```


Arguments are **named options**:

| Option | Description | Required | Default |
|--------|-------------|:--------:|---------|
| `--script "SCRIPT [ARGS]"` | Your annotated script and its arguments, quoted as one string. | Yes | n/a |
| `--tmp useTmp\|noTmp` | Copy each step's `reference` files to node-local `/tmp` (faster for big references) or not. | Yes | n/a |
| `--pipeline NAME` | Name stamped into the `pipeline` column of `jobRecord.txt`, which groups resource stats. Defaults to the pipeline-script basename; a generic basename triggers a warning suggesting you pass `--pipeline` so pipelines sharing a jobRecord store stay distinct. | No | script basename |
| `--sbatch-options "SBATCH_OPTIONS"` | Default sbatch options for steps that don't specify their own. | No | `"sbatch -p short -c 1 --mem 2G -t 50:0"` |
| `--mode dryrun` | Dry run only (builds the pipeline and logs planned jobs; no Slurm submission). | No | run (submit) |
| `--lint` | Validate block structure and exit — no conversion driver is run, nothing is submitted. | No | off |
| `--email noEmail\|noSuccEmail` | Silence all emails, or success emails only. | No | all emails |
| `--special checkpoint\|excludeFailedNodes` | Enable checkpointing, or exclude nodes where this job type failed before. | No | none |
| `-v\|--verbose` | Full firehose output (per-job "Submitted batch job N" + the jobs table). | No | off |
| `-q\|--quiet` | Errors only; silent on success. | No | off |

> [!NOTE]
> If `--sbatch-options` is omitted, or provided but doesn't start with `sbatch`, the default sbatch string above is used. That means every step must then get its resources either from that default or from its own `#@` line.

**Checking a pipeline without submitting: `--lint`.** `--lint` runs the conversion far enough to validate block structure, then stops — it never builds the run driver or submits anything.

```bash
runAsPipeline --script "myPipeline.sh args" --tmp noTmp --lint
```

It reports the block mode (explicit vs legacy), flags block-closure problems (orphan/duplicate `#@end`, a new `#@` before the previous block closed, an unclosed block at end-of-file), and — in explicit mode — notes any block still using a legacy blank-line close. Exit status is `0` when clean, non-zero on a structural error, so it drops into CI.

**Submission output and verbosity (`-v` / `-q`).** During Stage 2, `runAsPipeline` prints a compact **live** table — one row per sample, each step's glyph appearing as that job submits:

```text
Stage 2: Submitting jobs

  sample           steps
---------------------------------------------------------
  alpha            1✓ 2✓ 3✓ 4✓✓✓ 5✓ 6✓
  bravo            1✓ 2✓ 3✓ 4✓✓✓ 5✓ 7✓
  charlie          1✓ 2✓ 3✓ 4✓✓✓ 5✓ 6✓
Submitted 24, kept 0, error 0.  Details: less .smartSlurm.log  |  checkRun to monitor.
```

- `✓` = submitted, `✗` = submission error, `-` = kept (an already-successful step that was skipped).
- The step **number** prints when submission starts; its glyph is appended when Slurm returns.
- A fan-out step (multiple jobs at the same step) collapses under one number: `4✓✓✓` = step 4 submitted three chunk jobs.
- Rows are labelled by **sample** (the loop variable); an experiment-level step (no sample) is labelled by its job name.

Two flags change the volume:

| Flag | Output |
|------|--------|
| *(default)* | The compact live table above. |
| `-v` / `--verbose` | The old firehose — per-job "Working on inputs…", "Submitted batch job N", and the full jobs table. |
| `-q` / `--quiet` | Errors only; nothing on success. |

The full per-job detail is always written to `.smartSlurm.log` regardless of verbosity (`less .smartSlurm.log`).

**Example pipelines.** `example-scripts/` holds small, self-contained pipelines that double as templates and quick cluster smoke tests. Run each in its own directory.

| Pipeline | Shows |
|----------|-------|
| `demo_linear_fan_branch.sh` | the full progression — a linear chain, a fan-out/fan-in pair, and an `if/else` that submits one of two alternative steps. |
| `demo_phrase_pipeline.sh` | a minimal `#@end` linear pipeline with an experiment-level fan-in `summary` — the "start here" demo. |
| `bashScript-tutorial/` | the [Quick start](quick-start.md#quick-start) tutorial (`bashScriptV1`→`V2`, plus edge-case variants). |
| `resource-estimation/` | a size-scaled load (`sizedLoad.sh`) for exercising jobRecord fitting and OOM/OOT reruns. |

Each directory has its own README with the exact commands and expected output.

## ssbatch

For most programs, memory and run-time scale with input size. ssbatch records what resources prior jobs actually used, fits that relationship (see [Resource estimation](core-concepts.md#resource-estimation)), and uses it to allocate resources for subsequent jobs.

**Usage**

```
ssbatch [-P PROGRAM] [-R REF] [-I INPUTS] [-F FLAG] [-d DEPS] [SBATCH_OPTIONS] --wrap="COMMAND" [dryrun]
ssbatch [-P PROGRAM] [-R REF] [-I INPUTS] [-F FLAG] [-d DEPS] [SBATCH_OPTIONS] SCRIPT.sh [ARGS] [dryrun]
```

| Option | Description | Required |
|--------|-------------|----------|
| `-P PROGRAM` | Program name used to group resource statistics. If omitted, ssbatch derives it from the wrapped command or script name. | No |
| `-R REF` | Reference name (e.g. `hg19`) used in the stats grouping key. Without it, standalone jobs carry `ref=none`. | No |
| `-I INPUTS` | Input file(s)/dir(s), or an explicit size as `jobSize:N`. Drives resource estimation. | No |
| `-F FLAG` | Unique job identifier. If omitted, `program`+`input` (or a random suffix) is used. | No |
| `-d DEPS` | Slurm dependencies (used by `runAsPipeline`). | No |
| `--wrap="CMD"` | Command to run. | Yes* |
| `SCRIPT.sh [ARGS]` | A script to run instead of `--wrap`. Its first line must be a shebang. | Yes* |
| `SBATCH_OPTIONS` | Any standard sbatch option (`-c`, `--mem`, `-t`, `-p`, `-A`, `--mail-user=`, …). Defaults exist for memory and time. | No |
| `dryrun` | As the **last** argument: does not actually submit, only estimates and builds the job script. | No |

<sub>*Provide **either** `--wrap` **or** a script file.</sub>



> [!TIP]
> None of `-P`, `-I`, `-F` is mandatory, but supplying `-P` (and `-I` when a meaningful input exists) is what lets ssbatch build good per-program statistics. Without `-I`, estimation falls back to the distribution of that program's past run-times/memory rather than a size-based fit.

**Worked example**

```bash
# One-time setup (see "Installation and setup" above), then activate
git clone https://github.com/ld32/SmartSlurm.git $HOME/SmartSlurm
bash $HOME/SmartSlurm/smartSlurmSetup.sh
source ~/.bashrc          # load the new `smartslurm` command (or open a new shell)
smartslurm activate

# Make some test inputs
mkdir -p SmartSlurmTest && cd SmartSlurmTest
createNumberFiles.sh

# Run a few jobs so ssbatch can learn this program's memory/time profile.
# "findNumber" is just a label you choose with -P.
ssbatch -P findNumber -I numbers1.txt -F find1 --mem 4G -t 2:0:0 \
    --wrap="findNumber.sh 12345 numbers1.txt"

ssbatch -P findNumber -I numbers3.txt -F find3 --mem 4G -t 2:0:0 \
    --wrap="findNumber.sh 12345 numbers3.txt"

ssbatch -P findNumber -I numbers5.txt -F find5 --mem 4G -t 2:0:0 \
    --wrap="findNumber.sh 12345 numbers5.txt"
```

> [!IMPORTANT]
> Estimation needs **at least 3 completed records** for a given program/reference. Until then, jobs use the default memory and time. After 3 successful jobs, the next job is sized automatically — you can request `--mem 4G` and still see ssbatch submit it with, say, 21 M and a 13-minute limit.

```bash
# The 4th job is auto-sized from the first three:
ssbatch -P findNumber -I numbers2.txt -F find2 --mem 4G -t 2:0:0 \
    --wrap="findNumber.sh 12345 numbers2.txt"

# Multiple inputs are fine (their combined size is used):
ssbatch -P findNumber -I "numbers1.txt numbers2.txt" -F find12 --mem 4G -t 2:0:0 \
    --wrap="findNumber.sh 12345 numbers1.txt numbers2.txt"

# No -I? Estimation uses the 90th percentile of this program's past usage instead:
ssbatch -P findNumber -F find21 --mem 4G -t 2:0:0 \
    --wrap="findNumber.sh 12345 numbers2.txt"

# Check status
checkRun

# Cancel everything submitted from this directory
cancelAllJobs
```

> [!NOTE]
> Re-submitting a job with the **same program and input** that already succeeded will prompt you to confirm the rerun, so you don't silently repeat completed work.

## reviewJobRecords.py

`reviewJobRecords.py` is a terminal-only tool (no browser, X11, or extra Python packages) for pruning outliers from your records so estimates stay accurate. It groups records by `(pipeline, step, ref)` — the same key the estimator fits on.

```bash
reviewJobRecords.py                    # the default job record
reviewJobRecords.py path/to/jobRecord.txt
```

How it works:

1. Pick a program/step from the numbered list.
2. An ASCII scatter plot of **Input Size (G)** vs **Memory (G)** is drawn, each point labeled `0`–`9` then `a`–`z`, with a table of Job IDs and exact values below.
3. Delete outliers by index — single `2`, list `1,4,7`, range `3-6`, or mixed `0,2-4,8`.
4. The plot redraws immediately; remaining points **keep their original indices** and the **axis scale stays fixed**, so numbering doesn't shift under you.
5. `b` back to the program list · `s` save (a timestamped backup is written first) · `q` quit.


## cancelAllJobs

Cancels all running and pending jobs that `runAsPipeline` submitted from the current directory:

```bash
cancelAllJobs
```

## sbatchAndTop

Submit a job with `ssbatch` and immediately run `scontrol top` on it (with SmartSlurm active — see [Installation and setup](installation.md#installation-and-setup)):

```bash
sbatchAndTop -p short -c 1 -t 2:0:0 --mem 4G --wrap "my_application para1 para2"
# -p is optional — ssbatch picks the partition from the run-time.

# or with a script:
sbatchAndTop job.sh
```

## Under-the-hood tools

These are invoked by SmartSlurm itself (not typed by you day-to-day), but they are documented here because you may see them in logs or run them for maintenance. Descriptions are taken from each script's own header.

| Tool | Role |
|------|------|
| `cleanUp.sh` | Runs from each job's **EXIT trap** (on success and failure). It emits the `========== SmartSlurm job report ==========` section, writes the job's `efficiency.tsv` row (unconditionally) and its `jobRecord.txt` row, prunes/invalidates stale stats on OOM/OOT, decides whether to requeue with doubled resources, and sends the consolidated [emails](core-concepts.md#notifications-emails). |
| `estimateResource.sh` | The estimator. Called as `estimateResource.sh PROGRAM REF INPUTS …`; returns memory (in M) and time (in minutes) for a `(pipeline, step, ref)` group, fitting input-size→mem/time when ≥3 records exist, otherwise falling back to defaults. Writes fitted plots under `$smartSlurmJobRecordDir/stats/`. |
| `jobRecordSchema.sh` | **Single source of truth** for the jobRecord schema. It is *sourced* (defines variables and functions, no side effects) by the writers (`ssbatch`/`cleanUp.sh`), the readers (`estimateResource.sh`, `calculateMemTime.sh`), and `jrRepair.sh`, so the column layout is defined in exactly one place. It detects version by header signature (a v2 file's header begins with `pipeline,`; a legacy v1 header began with `1jobID,`). |
| `jrRepair.sh` | One-time, **manual** converter of a legacy (v1) jobRecord store to the v2 schema; run once per store after upgrading. It converts `jobRecord.txt` **and** `jobRecord.oom.oot.over.txt` in place after backing each up to `<base>.back.<timestamp>.txt`, and is idempotent (a v2 file is left untouched). Options: `--dry-run` (report, write nothing) and `--check` (health report only; exits non-zero if any store has unrecognized rows). Rows that are not clean 21-field v1 rows are copied to `<name>.jrRepair.skipped.<ts>` and counted, so nothing is silently lost. See [jobRecords](core-concepts.md#jobrecords). |

---


---

[&larr; SmartSlurm README](../README.md) &middot; [Docs index](README.md)
