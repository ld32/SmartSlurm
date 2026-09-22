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

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [ssbatch: smart sbatch](#ssbatch-smart-sbatch)
  - [Usage](#usage)
  - [Quick example](#quick-example)
  - [Utilities: unExport](#utilities-unexport)
  - [How ssbatch works](#how-ssbatch-works)
    - [jobRecord.txt](#jobrecordtxt)
    - [config.txt](#configtxt)
    - [Resource estimation](#resource-estimation)
    - [Partition selection](#partition-selection)
    - [OOM / OOT auto-resubmit](#oom--oot-auto-resubmit)
    - [Checkpointing](#checkpointing)
    - [Informative emails](#informative-emails)
  - [ssbatch FAQ](#ssbatch-faq)
- [runAsPipeline](#runaspipeline)
  - [Usage](#runaspipeline-usage)
  - [Writing a pipeline: the `#@` job block](#writing-a-pipeline-the--job-block)
  - [Closing a block: `#@end` vs blank line](#closing-a-block-end-vs-blank-line)
  - [Checking a pipeline without submitting: `--lint`](#checking-a-pipeline-without-submitting---lint)
  - [How runAsPipeline runs: the two phases](#how-runaspipeline-runs-the-two-phases)
  - [Passing files and variables between steps](#passing-files-and-variables-between-steps)
  - [Loops](#loops)
  - [Submission output and verbosity (`-v` / `-q`)](#submission-output-and-verbosity--v---q)
  - [Tutorial](#tutorial)
  - [Example pipelines](#example-pipelines)
  - [checkRun: monitor and debug](#checkrun-monitor-and-debug)
  - [diagnose.sh: explain one job](#diagnosesh-explain-one-job)
  - [Rerunning a pipeline](#rerunning-a-pipeline)
  - [cancelAllJobs](#cancelalljobs)
  - [runAsPipeline FAQ](#runaspipeline-faq)
- [Using ssbatch with other pipeline managers](#using-ssbatch-with-other-pipeline-managers)
  - [Snakemake](#snakemake)
  - [Nextflow](#nextflow)
  - [Cromwell](#cromwell)
- [Reviewing and cleaning job records](#reviewing-and-cleaning-job-records)
- [sbatchAndTop](#sbatchandtop)
- [Upgrading](#upgrading)

---

## Features

- **Auto-sizes memory and run-time** from statistics of your earlier jobs.
- **Auto-selects the partition** that matches the requested run-time.
- **Auto-resubmits** jobs that die out-of-memory (OOM) or out-of-time (OOT), with doubled resources.
- **Optional checkpointing**: snapshot a long job before it hits its limit, then resume from the snapshot.
- **Informative emails**: Slurm emails are just a subject line; SmartSlurm attaches the job script, the exact submit command, and the stdout/stderr logs.
- **(runAsPipeline) Dependency management**: steps wait for their prerequisites automatically.
- **(runAsPipeline) Explicit `#@ … #@end` blocks**: write a step's body verbatim to a per-job script (`awk`, quotes, heredocs, multi-line all survive), instead of flattening it into a single `--wrap` string.
- **(runAsPipeline) `--lint`**: validate a pipeline's block structure without submitting anything.
- **(runAsPipeline) Smart reruns**: an unchanged script is reused as-is; on re-issue you pick, from one menu, whether to redo everything, only failures, or from a chosen step.
- **(checkRun) Status grid + diagnosis**: a sample × step grid, addressable per cell, with inline "why did this fail" powered by `diagnose.sh`.

---

## Installation

```bash
# 1. Download SmartSlurm
git clone https://github.com/ld32/SmartSlurm.git $HOME/SmartSlurm

# 2. Put it on your PATH (add this line to ~/.bashrc to make it permanent)
export PATH=$HOME/SmartSlurm/bin:$PATH
```

**Optional** — only needed for the workflow-chart (`w`/`g`) option in `checkRun`:

```bash
module load conda/miniforge3/24.11.3-0
mamba create -n smartSlurmEnv -c conda-forge graphviz
```

> [!NOTE]
> The `module load conda/...` line is specific to an HPC that provides conda as an environment module. On other systems, activate conda however your site does. If your site uses a different conda module name, adjust it here and in `checkRun`.

---

# ssbatch: smart sbatch

For most programs, memory and run-time scale with input size. ssbatch records what resources prior jobs actually used, fits that relationship, and uses it to allocate resources for subsequent jobs.

<em><b>Figure 1</b> — Memory usage tracks input size, so input size is a good proxy for allocating memory.</em>
<div align="center">
<img src="https://github.com/ld32/SmartSlurm/blob/master/stats/back/findNumber.none.time.png" width="50%">
</div>

<em><b>Figure 2</b> — Early jobs run with the default memory; once enough have finished, ssbatch estimates memory for later jobs, sharply cutting wasted RAM.</em>
<div align="center">
<img src="https://github.com/ld32/SmartSlurm/blob/master/stats/back/barchartMem.png" width="50%">
</div>

<em><b>Figure 3</b> — The same for run-time: later jobs get much tighter time allocations.</em>
<div align="center">
<img src="https://github.com/ld32/SmartSlurm/blob/master/stats/back/barchartTime.png" width="50%">
</div>

## Usage

```
ssbatch [-P PROGRAM] [-I INPUTS] [-F FLAG] [SBATCH_OPTIONS] --wrap="COMMAND" [dryrun]
ssbatch [-P PROGRAM] [-I INPUTS] [-F FLAG] [SBATCH_OPTIONS] SCRIPT.sh [ARGS] [dryrun]
```

| Option | Description | Required |
|--------|-------------|----------|
| `-P PROGRAM` | Program name used to group resource statistics. If omitted, ssbatch derives it from the wrapped command or script name. | No |
| `-I INPUTS` | Input file(s)/dir(s), or an explicit size as `jobSize:N`. Drives resource estimation. | No |
| `-F FLAG` | Unique job identifier. If omitted, `program`+`input` (or a random suffix) is used. | No |
| `--wrap="CMD"` | Command to run. | Yes* |
| `SCRIPT.sh [ARGS]` | A script to run instead of `--wrap`. Its first line must be a shebang. | Yes* |
| `SBATCH_OPTIONS` | Any standard sbatch option (`-c`, `--mem`, `-t`, `-p`, `-A`, `--mail-user=`, …). Defaults exist for memory and time. | No |
| `dryrun` | As the **last** argument: does not actually submit, only estimates and builds the job script. | No |

<sub>*Provide **either** `--wrap` **or** a script file.</sub>

> [!TIP]
> None of `-P`, `-I`, `-F` is mandatory, but supplying `-P` (and `-I` when a meaningful input exists) is what lets ssbatch build good per-program statistics. Without `-I`, estimation falls back to the distribution of that program's past run-times/memory rather than a size-based fit.

## Quick example

```bash
# One-time setup
git clone https://github.com/ld32/SmartSlurm.git $HOME/SmartSlurm
export PATH=$HOME/SmartSlurm/bin:$PATH

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

## Utilities: unExport

`unExport` removes SmartSlurm from your `PATH` in the current shell, restoring the system `sbatch`:

```bash
source unExport; unExport
```

## How ssbatch works

### jobRecord.txt

A CSV of resource records for every successful job, one row per job. Its location is set by `smartSlurmJobRecordDir` in [`config.txt`](#configtxt) and defaults to `~/.smartSlurm/jobRecord.txt`.

> [!WARNING]
> The file has **18 columns** (created with the header below). Estimation matches on **program (col 12)** and **reference (col 13)** and reads **memory used (col 7)** and **time used (col 8)**. Input size (col 2) is used for the size-vs-resource fit.

```text
1jobID,2inputSize,3memDefault,4timeDefault,5memAllocated,6timeAllocated,7memUsed,8timeUsed,9jobStatus,10userID,11saccMem,12program,13reference,14flag,15core,16extraMem,17extraTime,18date
```

| Col | Name | Col | Name | Col | Name |
|----:|------|----:|------|----:|------|
| 1 | jobID | 7 | **memUsed** ⭐ | 13 | **reference** ⭐ |
| 2 | **inputSize** ⭐ | 8 | **timeUsed** ⭐ | 14 | flag |
| 3 | memDefault | 9 | jobStatus | 15 | core |
| 4 | timeDefault | 10 | userID | 16 | extraMem |
| 5 | memAllocated | 11 | saccMem | 17 | extraTime |
| 6 | timeAllocated | 12 | **program** ⭐ | 18 | date |

<sub>⭐ = used directly by resource estimation.</sub>

Example rows:

```text
46531,1465,4G,2:0:0,4G,0-2:0:0,3.52,1,COMPLETED,ld32,,findNumber,none,...
46535,2930,4G,2:0:0,4G,0-2:0:0,6.38,2,COMPLETED,ld32,,findNumber,none,...
46534,4395,4G,2:0:0,4G,0-2:0:0,9.24,4,COMPLETED,ld32,,findNumber,none,...
```

### config.txt

Ships at `SmartSlurm/config/config.txt`. It defines the partition time limits, the `adjustPartition` function, and defaults such as where records and logs live.

> [!TIP]
> Copy it to `~/.smartSlurm/config/config.txt` to override settings for just yourself. Your personal copy wins over the shared one — handy on a shared cluster. (SmartSlurm refuses to run if your personal copy is *older* than the shared one, to stop you using a stale config; see [Upgrading](#upgrading) to refresh it.)

Key settings:

```bash
export smartSlurmJobRecordDir=$HOME/.smartSlurm  # where jobRecord.txt + stats live
export smartSlurmLogDir=smartSlurmLog            # per-run log dir (relative to your cwd)

export firstBatchCount=5   # runAsPipeline: how many independent jobs run before the rest are held
export defaultMem=4096     # M   — used until estimation is available
export defaultTime=120     # min — used until estimation is available
export defaultExtraMem=500   # M   — safety margin added to estimates
export defaultExtraTime=10  # min — safety margin added to estimates

export partition1TimeLimit=720    # minutes: run-time  >0    and ≤12h  (720 min)
export partition2TimeLimit=7200   # minutes: run-time >12h   and ≤5 days (7200 min)
export partition3TimeLimit=43200  # minutes: run-time  >5d   and ≤30 days (43200 min)

adjustPartition() { ...; }
```

> [!IMPORTANT]
> **`firstBatchCount` (5) and the "3 records" rule are different things.**
> - **3 records** is the minimum ssbatch needs before it can *estimate* resources for a program.
> - **`firstBatchCount=5`** is used only by `runAsPipeline`: it lets the first 5 independent jobs of a step run immediately, and submits the rest **held** (`-H`) until early jobs finish and produce records, then releases them with estimated resources.

### Resource estimation

Performed by `estimateResource.sh`. Two modes:

- **Input size given (`-I`)** — ssbatch fits a straight line to *input size vs. memory* and *input size vs. time* using past records, and reads the estimate off the fit. Needs **≥3 records** to fit; below that, it uses the defaults.
- **No input size** — ssbatch takes the **90th percentile** of this program's past memory and time (so ~90% of jobs finish within the allocation). Also needs **≥3 records**; below that, defaults.

Fitted plots are written under `$smartSlurmJobRecordDir/stats/` for inspection:

<div align="center">
<img src="https://github.com/ld32/SmartSlurm/blob/master/stats/back/findNumber.none.mem.png" width="45%" style="display:inline-block; margin-right:2%;">
<img src="https://github.com/ld32/SmartSlurm/blob/master/stats/back/findNumber.none.time.png" width="45%" style="display:inline-block; margin-left:2%;">
</div>

### Partition selection

`ssbatch` (via `adjustPartition` in `config.txt`) chooses a partition from the requested run-time using the `partitionNTimeLimit` values above. If you pass `-p` yourself and it doesn't fit the run-time, it is adjusted for you — so you generally don't need to specify `-p` at all.

### OOM / OOT auto-resubmit

When a job ends, `cleanUp.sh` records its memory and time usage into `jobRecord.txt`. If the job **failed out-of-memory or out-of-time**, ssbatch resubmits it with **double** the memory or time, and clears the stored formula for that program/reference so later jobs re-learn from fresh data.

### Checkpointing

Optional. When enabled (the `checkpoint` mode in `runAsPipeline`, or a `.Checkpoint` program name), a long-running job is snapshotted before it would hit its memory/time limit and resubmitted to resume from the checkpoint, instead of restarting from scratch.

### Informative emails

Slurm's own email is only a subject line. `cleanUp.sh` instead emails you the **job script**, the **exact submit command used**, and the **stdout/stderr logs** — enough to diagnose a failure without logging in. Use `noSuccEmail` (failures only) or `noEmail` (none) to reduce volume.

## ssbatch FAQ

<details>
<summary><b>Do I have to wait for the first jobs to finish before later jobs get estimated resources?</b></summary>

For **standalone ssbatch**: there is no waiting — each job is submitted immediately. Estimation simply kicks in once ≥3 records exist. For **runAsPipeline**: yes, by design it runs the first `firstBatchCount` (5) independent jobs and holds the rest until estimates are available.
</details>

<details>
<summary><b>Are -P, -I, and -F optional?</b></summary>

Yes, all three.
- `-P` omitted → program name is taken from the wrapped command or the script name.
- `-I` omitted → estimation uses the program's past-usage distribution (90th percentile) instead of a size-based fit.
- `-F` omitted → the unique flag becomes `program`+`input`, or `program`+a random suffix.
</details>

<details>
<summary><b>Can -I take a size directly instead of a file?</b></summary>

Yes: `-I jobSize:12` — `12` is the input size (any integer). Useful when the meaningful "size" isn't a single file.
</details>

<details>
<summary><b>Can I pass -c and other sbatch options?</b></summary>

Yes. Any standard sbatch option works and is passed through.
</details>

<details>
<summary><b>How does ssbatch know a job already ran?</b></summary>

On success it creates `<flag>.success` in `smartSlurmLogDir`. That file's presence is how reruns are detected.
</details>

---

# runAsPipeline

`runAsPipeline` turns an ordinary Bash script into a dependency-aware pipeline. You add `#@` comment markers above the commands you want submitted as Slurm jobs; everything else stays a normal shell command. Each marked step is submitted through `ssbatch`, so it inherits smart sizing, resubmission, and emails.

**What it gives you**

- **Dependency management** — steps wait for the steps they depend on.
- **Per-step smart sizing** — every step is an `ssbatch` job.
- **Smart reruns** — an unchanged script is reused without reconversion; on re-issue you choose what to redo ([Rerunning a pipeline](#rerunning-a-pipeline)).
- **Failure containment** — if a step fails, its downstream steps are not run.
- **Multi-account support** — add `-A`/`--account=` and all jobs use that Slurm account.

## runAsPipeline Usage

```
runAsPipeline --script "SCRIPT [ARGS]" --tmp {useTmp|noTmp} [--sbatch-options "SBATCH_OPTIONS"]
              [--mode dryrun] [--lint] [--email noEmail|noSuccEmail]
              [--special checkpoint|excludeFailedNodes] [-v|--verbose | -q|--quiet]
```

Arguments are **named options**:

| Option | Description | Required | Default |
|--------|-------------|:--------:|---------|
| `--script "SCRIPT [ARGS]"` | Your annotated script and its arguments, quoted as one string. | Yes | n/a |
| `--tmp useTmp\|noTmp` | Copy each step's `reference` files to node-local `/tmp` (faster for big references) or not. | Yes | n/a |
| `--sbatch-options "SBATCH_OPTIONS"` | Default sbatch options for steps that don't specify their own. | No | `"sbatch -p short -c 1 --mem 2G -t 50:0"` |
| `--mode dryrun` | Dry run only (builds the pipeline and logs planned jobs; no Slurm submission). | No | run (submit) |
| `--lint` | Validate block structure and exit — no conversion driver is run, nothing is submitted. See [`--lint`](#checking-a-pipeline-without-submitting---lint). | No | off |
| `--email noEmail\|noSuccEmail` | Silence all emails, or success emails only. | No | all emails |
| `--special checkpoint\|excludeFailedNodes` | Enable checkpointing, or exclude nodes where this job type failed before. | No | none |
| `-v\|--verbose` | Full firehose output (per-job "Submitted batch job N" + the jobs table). | No | off |
| `-q\|--quiet` | Errors only; silent on success. | No | off |

> [!NOTE]
> If `--sbatch-options` is omitted, or provided but doesn't start with `sbatch`, the default sbatch string above is used. That means every step must then get its resources either from that default or from its own `#@` line.

## Writing a pipeline: the `#@` job block

A **job block** is one `#@` annotation line plus the command line(s) that follow it:

```
#@ stepID , dependIDs , name , reference , inputs , sbatchOptions
<the command(s) to run for this step>
```

| Field | Meaning | Example | Required |
|-------|---------|---------|:--------:|
| `stepID` | Unique integer identifying the step. | `1`, `2` | **Yes** |
| `dependIDs` | `0` for no dependency, or upstream step IDs joined by dots. | `0`, `1`, `1.3` | **Yes** |
| `name` | Program/label; groups resource stats and names the job. | `findNumber` | **Yes** |
| `reference` | Reference file(s)/dir(s), dot-joined; synced to `/tmp` under `useTmp`. | `genome.fa`, `db1.db2` | No |
| `inputs` | Input file(s), dot-joined; size drives resource estimation. | `sample.fq`, `in1.in2` | No |
| `sbatchOptions` | sbatch options for this step. Omit to use the command-line default. | `sbatch -c 4 -t 2:0:0` | No |

> [!IMPORTANT]
> **Only the first three fields (`stepID`, `dependIDs`, `name`) are required.** `reference`, `inputs`, and `sbatchOptions` may be left empty — just keep the commas as placeholders.

**Examples**

```bash
# Step 1, no dependency, program "findNumber", input $input, explicit resources:
#@1,0,findNumber,,input,sbatch -p short -c 1 --mem 4G -t 50:0
findNumber.sh 1234 $input > $number.$i.txt

# Step 2 depends on step 1; no reference, no input, no sbatch options
# (uses the command-line default sbatch string):
#@2,1,mergeNumber,,,
cat $number.*.txt > all$number.txt

# Step 4 depends on steps 1 AND 3; one reference to sync, one input, custom resources:
#@4,1.3,map,genome.fa,reads.fq,sbatch -p short -c 4 -t 2:0:0
map.sh genome.fa reads.fq > out.bam
```

> [!TIP]
> Every step's `stepID` must be unique, or the run aborts. In an `if/else` where a step exists in both branches, give the two branches **distinct** step IDs — reusing one ID across branches drops the second branch.

## Closing a block: `#@end` vs blank line

A block needs a clear end so the parser knows which lines belong to the job. There are **two ways**, and a pipeline uses one mode throughout (chosen automatically):

### Explicit `#@end` (recommended)

Close the block with a line containing exactly `#@end`. The body between the marker and `#@end` is written **verbatim** to a per-job script (`smartSlurmLog/jobs/<flag>/cmd.sh`) and that file is submitted:

```bash
#@3,2,measure,,in,sbatch -p short -c 1 --mem 100M -t 5:00
# any awk, quotes, multi-line, heredocs, and command substitution survive intact:
awk '{ print length(\$1)"\t"\$1 }' tokens.$s.txt | sort -k1,1nr -k2,2 > measured.$s.txt
#@end
```

Because the body isn't flattened into a `--wrap` string, this is the robust choice for anything beyond a trivial one-liner. **If a pipeline uses `#@end` anywhere, it is in explicit mode: blank lines are ordinary body, and *only* `#@end` closes a block.**

Two escaping rules apply to the body of an `#@end` block:

- A `$` that must reach the **tool** literally (an `awk` field like `$1`, a `perl` `$_`) must be written **`\$`** — otherwise it is read as the pipeline's own shell variable. (See the `\$1` in the example above.)
- A variable **assigned inside the block** expands on the node; a variable from **submit scope** (including the loop variable) is baked in at submission. See [Passing files and variables between steps](#passing-files-and-variables-between-steps).

> [!TIP]
> Avoid a `for`/`while` loop *inside* a block body — a loop variable assigned inside a block isn't recognized as node-scope. Use explicit names, or put the loop in submit scope wrapping the `#@` block (that's how fan-out works; see [Loops](#loops)).

### Legacy blank-line close

If a pipeline contains **no** `#@end`, blocks close the old way: a `#@` block collects every following line into one command and keeps going until it reaches a **blank line** (empty or all-whitespace). The block is flattened into a single `--wrap` string.

> [!IMPORTANT]
> In legacy mode a **blank line is the only thing that ends a block** — leave one before a `done`, the next `#@`, or any following plain command, or it gets swept into the previous job's command.

How each kind of comment is treated **in legacy mode**:

| Comment style | Ends the block? | In the job command? | Notes |
|---|:---:|---|---|
| Full-line `# comment` | No | No | Passed through to the converted script, but not part of the job's command. |
| Trailing `code # comment` | No | Only the `code` before ` #` | Everything from the first space-`#` to end of line is stripped. |
| Heredoc `: << EOF … EOF` | No | Mangled | **Not supported in legacy mode** — flattened into the command and breaks the block. Use `#@end` if you need a heredoc. |
| **Blank line** | **Yes** | n/a | The intended, and only, way to close a legacy block. |

> [!WARNING]
> Legacy inline-comment stripping is **textual, not shell-aware**: everything from the first space-`#` (` #`) to end of line is removed. Avoid a literal space-`#` inside a legacy command — *even inside quotes* — or it is truncated (e.g. `sed 's/ #/x/'` becomes `sed 's/`). A `#` with no space before it (`grep '#'`) is safe. The legacy `--wrap` flattener also **collapses quoted whitespace** (e.g. `tr -s ' '` loses its argument). Both problems disappear with `#@end`, which writes the body verbatim — prefer `#@end` for any non-trivial body.

## Checking a pipeline without submitting: `--lint`

`--lint` runs the conversion far enough to validate block structure, then stops — it never builds the run driver or submits anything.

```bash
runAsPipeline --script "myPipeline.sh args" --tmp noTmp --lint
```

It reports the block mode (explicit vs legacy), flags block-closure problems (orphan/duplicate `#@end`, a new `#@` before the previous block closed, an unclosed block at end-of-file), and — in explicit mode — notes any block still using a legacy blank-line close. Exit status is `0` when clean, non-zero on a structural error, so it drops into CI.

## How runAsPipeline runs: the two phases

Understanding **when** each line runs is critical when adapting scripts to use `runAsPipeline`. Scripts are parsed in two phases.

```
                 ┌─────────────────────── PHASE 1: CONVERT (once, on submit host) ──────────────────────┐
  your_script.sh │ read top → bottom:                                                                    │
                 │   • for/while/#loopStart  → remember the loop variable (becomes part of each job flag) │
                 │   • #@ marker + body      → turn into an ssbatch call (--wrap, or a per-job cmd.sh)    │
                 │   • any other line         → copy through unchanged                                    │
                 └──────────────────────────────────────────┬───────────────────────────────────────────┘
                                                             ▼
                       smartSlurmLog/slurmPipeLine.<md5>.sh  (the "converted script")
                                                             │
                 ┌───────────────────────── PHASE 2: SUBMIT (run the converted script) ──────────────────┐
                 │   plain lines  → run right now, on the submit host                                     │
                 │   ssbatch call → SUBMIT a job, capture its ID; the command runs LATER, on a node,      │
                 │                  once dependencies (-d afterok:…) are satisfied.                        │
                 │                  Independent jobs beyond firstBatchCount(5) are submitted HELD (-H),    │
                 │                  then released after early jobs produce ≥3 records to estimate from.    │
                 └───────────────────────────────────────────────────────────────────────────────────────┘
```

**Phase 1 — Convert.** `runAsPipeline` reads your script once and writes a converted script named `smartSlurmLog/slurmPipeLine.<md5>.sh`. The `<md5>` is a hash of your script's contents, so re-running with an **unchanged** script skips reconversion and reuses the existing one.

**Phase 2 — Submit.** The converted script is executed top to bottom. Plain shell lines run **immediately, on the submit host**. Each `ssbatch` call **submits** a job and captures its ID — it does **not** run your command inline. Your command runs later, on a compute node, after its dependencies complete.

> [!WARNING]
> Because plain lines run at **submission time**, any shell logic that inspects a job's **output** runs *before that output exists*. The check below always prints `process`, never `skip`:
> ```bash
> for sample in sample1 sample2; do
>     #@1,0,process,,sample
>     process_sample.sh ${sample}.fq > ${sample}.result
>
>     # Runs during Phase 2 on the submit host — BEFORE job #@1 executes on a node.
>     # ${sample}.result does not exist yet, so this is always "process".
>     [ -f ${sample}.result ] && status=skip || status=process
>     echo "$sample: $status" >> summary.txt
> done
> ```
> **Fix:** move any logic that depends on a step's results into a *later* `#@` step, so it runs on a node after the upstream job finishes. The same rule governs an `if/else` that chooses which step to submit — the **condition must be evaluable at submission time** (from the input/manifest), not from a job's output.

> [!TIP]
> Functions you define in the plain part of the script are **not** automatically available inside `#@` job blocks (those run in separate jobs). Export them first: `export -f myfunction`.

## Passing files and variables between steps

Every runAsPipeline script has **two scopes**. Knowing which one a variable lives in tells you exactly *when* it is evaluated and *who* can see it.

| Scope | What lives here | When it runs | Who can see it |
|-------|-----------------|--------------|----------------|
| **Submit scope** | Every plain line, including variable assignments *outside* any `#@` block | Once, on the submit host, during submission ([Phase 2](#how-runaspipeline-runs-the-two-phases)) | Later plain lines, the `reference`/`inputs`/`sbatch` fields of markers, and the command text of *every* `#@` block (substituted at submission) |
| **Node scope** | Code *inside* a `#@` block | Later, on a compute node, when the job runs | Only that one job |

The rule for when a `$var` is expanded:

> [!IMPORTANT]
> - A `$var` **used in a `#@` block but assigned in submit scope** is substituted with its value **at submission time** — the literal value is baked into the job's command before it goes to a node.
> - A `$var` **both assigned and used inside the same `#@` block** is left alone and expands **on the node at run time** (runAsPipeline sees the in-block assignment and defers it).
> - A variable assigned inside a `#@` block **exists nowhere else** — not in submit scope, not in other jobs.

### Data crosses between steps through files, not variables

Because each job runs in its own node scope, you cannot hand a value from one job to the next through a shell variable. Steps communicate through **files on shared storage**:

1. Name the output path once, as a **plain variable in submit scope**.
2. The **producer** step writes to that path; the **consumer** step reads from it.
3. Declare the dependency in the consumer's marker so it waits for the producer.

The same submit-scope variable is substituted into both commands at submission, so both refer to the identical literal path — while the data itself flows through the filesystem, gated by the dependency:

```bash
mapInputR1=$outDir/fastq/$sampleName.1.noadap.fastq   # submit scope: name the file once
mapInputR2=$outDir/fastq/$sampleName.2.noadap.fastq

#@1,0,cutadaptSeqtk,,fileR1.fileR2,sbatch -c 4 -p short -t 3:0:0 --mem 12G
... ; seqtk trimfq -e 1 ...trim.paired.fastq > $mapInputR1; cp ... $mapInputR2   # producer writes the files

bowtieOut=$outDir/mapping/${refPrefix}_$genomeSpike.bam                          # submit scope: name step 3's output
refIndex=$bwtPath/${genomeRef}_$genomeSpike/${genomeRef}_$genomeSpike

#@3,1.2,bowtie2,refIndex,mapInputR1.mapInputR2,sbatch -c 6 -p short -t 2:00:0 --mem 24G
... bowtie2 ... -1 $mapInputR1 -2 $mapInputR2 ... | samtools sort ... -o $bowtieOut ...   # consumer reads the files
```

At submission `$mapInputR1` expands to the same path in step 1's and step 3's commands. Step 3 depends on steps 1 and 2 (`#@3,1.2,...`), so it is held until they finish — only then does the file exist, get measured, and step 3 get sized and released.

> [!NOTE]
> Listing a **produced** file as a downstream step's `inputs` is safe even though it doesn't exist at submission time. Dependent steps are submitted **held**, and their resource estimation is deferred until the upstream step completes and the file exists.

### Recommended convention: name paths in submit scope, next to the step

Assign each step's input/output paths as plain variables **immediately above the `#@` marker that first uses them**. This gives one source of truth per path, keeps producer and consumer in agreement automatically, and reads top-to-bottom.

> [!WARNING]
> **Submit-scope variables persist across steps, branches, and loop iterations.** A value set for one step is still set when a later step is submitted, so a *missing* assignment silently reuses a stale value. Two habits prevent this:
> - Assign the variable on **every branch** that leads to a step using it.
> - Inside a **per-sample loop**, (re)assign paths **within the loop**, next to the step — not above the loop — or every iteration will submit jobs pointing at the first iteration's paths.

### Anti-patterns

```bash
# ✗ Using a variable that was set inside a job block, from outside that block
#@1,0,stepA,,in
out=$outDir/result.bam; process.sh > $out     # 'out' lives only on step A's node
#@2,1,stepB,,in
summarize.sh $out                             # empty here: 'out' was never in submit scope
```
Fix: define `out` in submit scope (above `#@1`) so both steps see the same path.

```bash
# ✗ Reading a job's output during submission
producedByStep1=$outDir/a.bam
nReads=$(samtools view -c $producedByStep1)   # runs at SUBMISSION, before step 1 has made the file
#@2,1,stepB,,in
do_something.sh $nReads                        # $nReads is empty/garbage
```
Fix: move the `samtools view -c` **inside** step 2's block, where it runs on a node after step 1 completes.

## Loops

`runAsPipeline` preserves loop structure but extracts the `#@` blocks inside. **Every enclosing loop variable becomes part of each job's flag**, so per-iteration jobs get distinct names (e.g. `1.0.findNumber.1`, `1.0.findNumber.2`, …), and a **nested** loop (fan-out) stays collision-free (`4.3.countchunk.$sample.$chunk`).

**`for` loops** work directly — the variable right after `for` is detected automatically:

```bash
for file in `ls someFolder`; do
    #@1,0,process,,file
    process.sh $file
done
# each iteration submits a job flagged ...process.$file
```

**Fan-out / fan-in.** A nested `for` loop around a `#@` block submits one job per iteration (fan-out); a later step that depends on it waits for all of them (fan-in):

```bash
for c in 1 2 3; do
    #@4,3,countchunk,,in,sbatch -p short -c 1 --mem 100M -t 5:00
    wc -w < chunk.$s.$c > cc.$s.$c.txt      # 3 jobs: ...countchunk.$s.1/.2/.3
    #@end
done
#@5,4,merge,,in,sbatch -p short -c 1 --mem 100M -t 5:00
cat cc.$s.*.txt > merged.$s.txt             # depends on step 4 -> waits for all three
#@end
```

**`while` loops need a hint.** A `while` header doesn't name its loop variable in a position the parser can read, so declare it with `#loopStart:VAR` on the line above:

```bash
#loopStart:f1
while read -r f1 f2 f3 f4; do
    #@1,0,process,,f1
    process.sh "$f1"
done < samples.txt
```

## Submission output and verbosity (`-v` / `-q`)

During Stage 2, `runAsPipeline` prints a compact **live** table — one row per sample, each step's glyph appearing as that job submits:

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

## Tutorial

The tutorial scripts (`bashScriptV1.sh` … `bashScriptV3.sh`) live in `example-scripts/bashScript-tutorial/`. Run the commands below from that directory, or copy the scripts into your working directory. They rely on `createNumberFiles.sh` and `findNumber.sh` (both in `bin/`, on your PATH).

Start from a plain script, `bashScriptV1.sh`:

```bash
#!/bin/sh
number=$1
[ -z "$number" ] && echo -e "Error: number is missing.\nUsage: bashScript <number>" && exit 1

for i in {1..5}; do
    input=numbers$i.txt
    findNumber.sh 1234 $input > $number.$i.txt
done

cat $number.*.txt > all$number.txt
```

It searches for a number in `numbers1.txt … numbers5.txt`, then merges the results. To run the search step and the merge step as Slurm jobs, add `#@` markers — this is `bashScriptV2.sh`:

```bash
#!/bin/bash
number=$1
[ -z "$number" ] && echo -e "Error: number is missing.\nUsage: bashScript <number>" && exit 1

for i in {1..5}; do
    input=numbers$i.txt
    #@1,0,findNumber,,input,sbatch -p short -c 1 --mem 4G -t 50:0
    findNumber.sh 1234 $input > $number.$i.txt
done

#@2,1,mergeNumber,,,sbatch -p short -c 1 --mem 4G -t 50:0
cat $number.*.txt > all$number.txt
```

**Reading the markers**

- `#@1,0,findNumber,,input,sbatch …` — step **1**, depends on **nothing** (`0`), program **findNumber**, **no** reference, input is `$input`. Inside the `for` loop, it submits five jobs, one per `$i`.
- `#@2,1,mergeNumber,,,sbatch …` — step **2**, depends on **step 1**. Slurm holds it until all five step-1 jobs finish.

**Dry run** (adds `--mode dryrun`, so nothing is submitted — you see the plan and fake IDs):

```bash
runAsPipeline --script "bashScriptV2.sh 123" --sbatch-options "sbatch -p short -t 10:0 -c 1" --tmp useTmp --mode dryrun
```

**Real run** (default mode, no `--mode` needed):

```bash
runAsPipeline --script "bashScriptV2.sh 1234" --sbatch-options "sbatch -p short -t 10:0 -c 1" --tmp useTmp
```

Its Stage-2 output is the compact table described in [Submission output](#submission-output-and-verbosity--v---q) — one row per `$i`, five `findNumber` submissions plus the `mergeNumber` row. Add `-v` to see the full per-job "Submitted batch job N" firehose and the jobs table:

```text
Stage 2: Submitting jobs
...
step: 1, depends on: 0, job name: findNumber, flag: 1.0.findNumber.1
Submitted batch job 69308
...
step: 2, depends on: 1, job name: mergeNumber, flag: 2.1.mergeNumber
Submitted batch job 69313

All submitted jobs:
job_id       depend_on                      job_flag          program     reference  inputs
69308        null                           1.0.findNumber.1  findNumber  none       numbers1.txt
...
69313        69308:69309:69310:69311:69312  2.1.mergeNumber   mergeNumber none       none
```

After it runs:

```bash
ls -l smartSlurmLog   # per-step .sh (job scripts), .out (logs), .success / .failed flags,
                      # and (for #@end steps) a jobs/<flag>/cmd.sh
checkRun              # interactive status grid + log/diagnosis browser
cancelAllJobs         # cancel running/pending jobs from this directory
```

## Example pipelines

`example-scripts/` holds small, self-contained pipelines that double as templates and quick cluster smoke tests. Run each in its own directory.

| Pipeline | Shows |
|----------|-------|
| `demo_linear_fan_branch.sh` | the full progression — a linear chain, a fan-out/fan-in pair, and an `if/else` that submits one of two alternative steps. |
| `phrase_pipeline_end.sh` | a minimal `#@end` pipeline with an experiment-level fan-in `summary`. |
| `phrase_pipeline_legacy.sh` | the legacy blank-line (`--wrap`) style and its constraints. |
| `bashScript-tutorial/` | the tutorial above (`bashScriptV1`→`V2`, plus edge-case variants). |
| `resource-estimation/` | a size-scaled load (`sizedLoad.sh`) for exercising jobRecord fitting and OOM/OOT reruns. |

Each directory has its own README with the exact commands and expected output.

## checkRun: monitor and debug

`checkRun` is an interactive monitor. **Run it from the directory where you launched `runAsPipeline`** (it reads `.smartSlurm.log` and `smartSlurmLog/` there):

```bash
checkRun
```

It opens on a **status grid** — samples down the side, steps across the top, one glyph per cell:

```text
Status:  ✓ 12 done   ✗ 0 error   ▸ 0 running   · 12 pending

   row  sample      normali tokeniz  split  countch  merge  detaile  brief
                       1       2       3       4       5       6       7
   a    alpha      │   +   │   +   │   +   │   @   │   +   │   +   │   ·   │
   b    bravo      │   +   │   +   │   +   │   @   │   +   │   ·   │   +   │
   ...
```

| Glyph | Meaning |
|:-----:|---------|
| `+` | done (`.success`) |
| `X` | error (`.failed`) |
| `x` | orphaned (an upstream dependency didn't succeed) |
| `>` | running |
| `-` | pending |
| `?` | unknown (no flag and not in the queue — often killed / node failure) |
| `@` | multi-job cell (a fan-out step; resolving its address lists all the jobs) |
| `·` | not applicable (this sample never ran this step — e.g. the untaken `if/else` branch) |

**Addressing a cell.** A job address is a **row letter + step number** — `b2` is row b, step 2; `e3` is experiment-level step 3. Both orders parse (`b2` or `2b`). At the prompt:

| Type | Action |
|------|--------|
| `b2` | view that job's log (the `.out`) |
| `d b2` | **diagnose** — why it failed / what to do (runs [`diagnose.sh`](#diagnosesh-explain-one-job)) |
| `a b2` | ask the AI helper about the job |
| `r b2` | view the job's raw files (`.out` / `.err` / `.sh` …) |
| `o` | browse other SmartSlurm log folders in this directory |
| `l` | list view (the classic per-job browser; all its capabilities are preserved) |
| `g` | workflow chart (DAG image; needs graphviz — see [Installation](#installation)) |
| `p` | show/hide pending jobs |
| `Enter` | refresh (re-query `squeue`, redraw) |
| `q` / `qq` | back up a level / quit |

**Non-interactive verbs.** The same grid and diagnosis are available without entering the UI:

```bash
checkRun grid              # print the status grid and exit
checkRun why b2            # diagnose the job at address b2 and exit
checkRun log b2            # print the path to that job's .out
checkRun --json            # a diagnose.sh JSON record for every job (for scripts)
```

`checkRunGrid` is the underlying render/resolve engine; `checkRun` calls it.

**Typical debugging flow:** open `checkRun`, find a cell showing `X` (or `x`/`?`), type `d <addr>` to get the cause and suggested action, then `<addr>` to read the log.

## diagnose.sh: explain one job

`diagnose.sh` answers "what happened, why, and what to do" for a single job. It is authoritative about *whether* a job succeeded (`.success` + sacct) and heuristic about *why* (an ordered pattern table over the tool's output region), and it never asserts a cause without an evidence line.

```bash
diagnose.sh <flag> [--log DIR] [--json]
# e.g.
diagnose.sh 1.0.normalize.alpha
diagnose.sh 1.0.normalize.alpha --json     # machine-readable record
```

The human report gives a one-line state/action header (e.g. `[OK]`, `[RETRYABLE]`, `[DOOMED]`, `[INSPECT]`), a plain-language summary, the evidence line it keyed on (with file:line), and the log path. The `--json` record carries fields `state, category, action, retryable, evidence, summary, log` — consume it by field name.

## Rerunning a pipeline

Re-issue the **same** `runAsPipeline` command to rerun. The unchanged script reuses its converted driver (same md5), so failed and never-run steps rerun automatically; already-successful steps prompt once, up front, with a single menu:

| Choice | Effect |
|--------|--------|
| `all` | re-run everything, including completed steps. |
| `failed` *(default)* | re-run only failed/incomplete steps; keep successful ones. |
| `step` → *N* | re-run from a chosen step onward (offers only steps that actually completed). |
| `cancel` | do nothing. |

For scripts/CI, set the choice non-interactively with the (ephemeral) `smartSlurmRerun` environment variable:

```bash
smartSlurmRerun=failed        runAsPipeline --script "..." ...   # default
smartSlurmRerun=all           runAsPipeline --script "..." ...
smartSlurmRerun=fromStep:3    runAsPipeline --script "..." ...
```

Rerunning a step re-runs its downstream steps too (their dependency is satisfied by a fresh upstream run).

## cancelAllJobs

Cancels all running and pending jobs that `runAsPipeline` submitted from the current directory:

```bash
cancelAllJobs
```

## runAsPipeline FAQ

<details>
<summary><b>Do later jobs wait for the first jobs before getting estimated resources?</b></summary>

If a step's jobs are independent, `runAsPipeline` submits them all at once, lets the first `firstBatchCount` (5) run, and holds the rest. Once early jobs finish and produce ≥3 records, the held jobs are released with estimated resources.
</details>

<details>
<summary><b>Should I use `#@end` or the blank-line close?</b></summary>

Prefer **`#@end`** for anything beyond a trivial one-liner: the body is written verbatim to a per-job script, so `awk`, quotes, heredocs, multi-line bodies, and command substitution all work. The legacy blank-line close flattens the body into a single `--wrap` string, which mangles quoted whitespace and multi-line bodies. A pipeline is all-one-mode: if it uses `#@end` anywhere, every block must close with `#@end`.
</details>

<details>
<summary><b>Can inputs be given as a size instead of a file?</b></summary>

Yes. Assign a variable and reference it in the marker:
```bash
someVar=jobSize:12   # 12 is the input size (any integer)
#@1,0,runshard,,someVar,sbatch -p short -c 4 -t 0-12:00 --mem 8G
```
</details>

<details>
<summary><b>Multiple inputs?</b></summary>

Yes — dot-join them in the marker (`#@2,1,find,,input1.input2,...`) or use a shell variable holding a space-separated list.
</details>

<details>
<summary><b>Fewer or no emails?</b></summary>

Add `noSuccEmail` (failures only) or `noEmail` (none) with `--email`.
</details>

<details>
<summary><b>Can I drop the command-line sbatch options?</b></summary>

Yes, **if every step sets its own** `sbatchOptions`. Then omit `--sbatch-options`.
</details>

<details>
<summary><b>Does runAsPipeline run my script's commands in their original order?</b></summary>

Plain commands run top to bottom during submission. Commands under a `#@` marker do **not** run then — they are submitted as Slurm jobs and run later on compute nodes (respecting dependencies). See [the two phases](#how-runaspipeline-runs-the-two-phases).
</details>

<details>
<summary><b>How does runAsPipeline detect an already-completed step?</b></summary>

By the `<flag>.success` file in `smartSlurmLogDir` (set in `config.txt`).
</details>

---

# Using ssbatch with other pipeline managers

Because `ssbatch` accepts standard sbatch syntax, you can point Snakemake, Nextflow, or Cromwell at it as their submit command and get smart sizing for free.

## Snakemake

```bash
git clone https://github.com/ld32/SmartSlurm.git $HOME/SmartSlurm

mkdir -p SmartSlurmTest && cd SmartSlurmTest
export PATH=$HOME/SmartSlurm/bin:$PATH

cp $HOME/SmartSlurm/bin/Snakefile .
cp $HOME/SmartSlurm/config/config.yaml .

module load conda/miniforge3/24.11.3-0
mamba env create --name snakemakeEnv --file $HOME/SmartSlurm/config/snakemakeEnv.yaml
conda activate snakemakeEnv

# Use ssbatch as the cluster submit command:
snakemake -p -j 999 --latency-wait=80 --cluster "ssbatch -t 100 --mem 1G -p short"

# With a specific Slurm account:
snakemake -p -j 999 --latency-wait=80 --cluster "ssbatch -A mySlurmAccount -t 100 --mem 1G"

checkRun
```

## Nextflow

```bash
git clone https://github.com/ld32/SmartSlurm.git $HOME/SmartSlurm

mkdir -p SmartSlurmTest && cd SmartSlurmTest
module load conda/miniforge3/24.11.3-0
mamba create -n nextflowEnv -c bioconda -y nextflow

export PATH=$HOME/SmartSlurm/bin:$HOME/SmartSlurm/sbatchBin:$PATH
conda activate nextflowEnv

cp $HOME/SmartSlurm/bin/nextflow.nf .
cp $HOME/SmartSlurm/config/nextflow.config .

# For a specific Slurm account, edit nextflow.config and uncomment/set:
#   process.clusterOptions = '--account=mySlurmAcc'

nextflow run nextflow.nf -profile slurm
checkRun

# When finished, restore the system sbatch:
export PATH="${PATH/:$HOME\/SmartSlurm\/sbatchBin/}"
```

## Cromwell

> [!NOTE]
> Cromwell support is planned and not yet documented.

---

# Reviewing and cleaning job records

`reviewJobRecords.py` is a terminal-only tool (no browser, X11, or extra Python packages) for pruning outliers from your records so estimates stay accurate.

```bash
reviewJobRecords.py                    # the default job record
reviewJobRecords.py path/to/jobRecord.txt
```

How it works:

1. Pick a program from the numbered list.
2. An ASCII scatter plot of **Input Size (G)** vs **Memory (G)** is drawn, each point labeled `0`–`9` then `a`–`z`, with a table of Job IDs and exact values below.
3. Delete outliers by index — single `2`, list `1,4,7`, range `3-6`, or mixed `0,2-4,8`.
4. The plot redraws immediately; remaining points **keep their original indices** and the **axis scale stays fixed**, so numbering doesn't shift under you.
5. `b` back to the program list · `s` save (a timestamped backup is written first) · `q` quit.

---

# sbatchAndTop

Submit a job with `ssbatch` and immediately run `scontrol top` on it:

```bash
export PATH=$HOME/SmartSlurm/bin:$PATH

sbatchAndTop -p short -c 1 -t 2:0:0 --mem 4G --wrap "my_application para1 para2"
# -p is optional — ssbatch picks the partition from the run-time.

# or with a script:
sbatchAndTop job.sh
```

---

# Upgrading

Back up your config before upgrading, then pull:

```bash
cd ~/SmartSlurm
cp ~/.smartSlurm/config/config.txt ~/.smartSlurm/config/config.txt.backup
git pull

# then re-apply any personal changes to the new config.txt manually
```
