# SmartSlurm

**SmartSlurm** automatically estimates and optimizes memory and run-time for Slurm jobs to increase resource efficiency and decrease the amount of time jobs spend in pending status. It is pure Bash — no daemon, no database, no special runtime. Setup is to clone the repo and put it on your `PATH`.

It has two parts you can use together or separately:

| Tool | What it is | Use it when |
|------|-----------|-------------|
| **`ssbatch`** | A drop-in `sbatch` wrapper that estimates memory/time from your past jobs, picks a partition, resubmits on OOM/OOT, and sends informative emails. | You submit individual jobs and want smart resource sizing. |
| **`runAsPipeline`** | A workflow runner built on `ssbatch`. You annotate an ordinary Bash script with `#@` markers; it wires up dependencies and submits each step as its own smart job. | You have a multi-step pipeline with dependencies. |

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
  - [How runAsPipeline runs: the two phases](#how-runaspipeline-runs-the-two-phases)
  - [Loops](#loops)
  - [Tutorial](#tutorial)
  - [checkRun: monitor and debug](#checkrun-monitor-and-debug)
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
- **(runAsPipeline) Smart reruns**: an unchanged script is reused as-is, and already-successful steps are skipped unless you ask to rerun them.

---

## Installation

```bash
# 1. Download SmartSlurm
git clone https://github.com/ld32/SmartSlurm.git $HOME/SmartSlurm

# 2. Put it on your PATH (add this line to ~/.bashrc to make it permanent)
export PATH=$HOME/SmartSlurm/bin:$PATH
```

**Optional** — only needed for the workflow-chart (`w`) option in `checkRun`:

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
| `dryrun` | As the **last** argument: not actually submit, only estimates and builds the job script but does not submit. | No |

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
> Copy it to `~/.smartSlurm/config/config.txt` to override settings for just yourself. Your personal copy wins over the shared one — handy on a shared cluster. (SmartSlurm refuses to run if your personal copy is *older* than the shared one, to stop you using a stale config; run [`upgrade.sh`](#upgrading) to refresh it.)

Key settings:

```bash
export smartSlurmJobRecordDir=$HOME/.smartSlurm  # where jobRecord.txt + stats live
export smartSlurmLogDir=smartSlurmLog            # per-run log dir (relative to your cwd)

export firstBatchCount=5   # runAsPipeline: how many independent jobs run before the rest are held
export defaultMem=4096     # M   — used until estimation is available
export defaultTime=120     # min — used until estimation is available
export defaultExtraMem=500   # M   — safety margin added to estimates
export defaultExtraTime=10  # min — safety margin added to estimates

export partition1TimeLimit=720   # hours: run-time  >0h  and ≤12h
export partition2TimeLimit=7200  # hours: run-time >12h  and ≤5 days
export partition3TimeLimit=43200  # hours: run-time  >5d  and ≤30 days

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
- **Smart reruns** — an unchanged script is reused without reconversion; already-successful steps are skipped unless you choose to rerun.
- **Failure containment** — if a step fails, its downstream steps are not run.
- **Multi-account support** — add `-A`/`--account=` and all jobs use that Slurm account.

## runAsPipeline Usage

```
runAsPipeline --script "SCRIPT [ARGS]" --tmp {useTmp|noTmp} [--sbatch-options "SBATCH_OPTIONS"] [--mode dryrun] [--email noEmail|noSuccEmail] [--special checkpoint|excludeFailedNodes]
```

Arguments are **named options**:

| Option | Description | Required | Default |
|--------|-------------|:--------:|---------|
| `--script "SCRIPT [ARGS]"` | Your annotated script and its arguments, quoted as one string. | Yes | n/a |
| `--tmp useTmp\|noTmp` | Copy each step's `reference` files to node-local `/tmp` (faster for big references) or not. | Yes | n/a |
| `--sbatch-options "SBATCH_OPTIONS"` | Default sbatch options for steps that don't specify their own. | No | `"sbatch -p short -c 1 --mem 2G -t 50:0"` |
| `--mode dryrun` | Dry run only (builds the pipeline and logs planned jobs; no Slurm submission). | No | run (submit) |
| `--email noEmail\|noSuccEmail` | Silence all emails, or success emails only. | No | all emails |
| `--special checkpoint\|excludeFailedNodes` | Enable checkpointing, or exclude nodes where this job type failed before. | No | none |

> [!NOTE]
> If `--sbatch-options` is omitted, or provided but doesn't start with `sbatch`, the default sbatch string above is used. That means every step must then get its resources either from that default or from its own `#@` line.

## Writing a pipeline: the `#@` job block

A **job block** is one `#@` annotation line plus the command line(s) directly beneath it:

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

# Step 3 depends on steps 1 AND 2; two references, no input, default sbatch:
#@3,1.2,align,db1.db2,,
align.sh $db1 $db2
```

> [!TIP]
> The multi-line command under a `#@` marker can be split with trailing backslashes. Every step's `stepID` must be unique, or the run aborts.

### Comments, and where a job block ends

> [!IMPORTANT]
> A `#@` block collects **every following line into one command** and keeps going until it reaches a **blank line** (empty or all-whitespace). **A blank line is the only thing that ends a block.** This is why you must leave a blank line before a `done`, the next `#@`, or any following plain command — otherwise it gets swept into the previous job's command.

How each kind of comment is treated:

| Comment style | Ends the block? | In the job command? | Notes |
|---|:---:|---|---|
| Full-line `# comment` | No | No | Passed through to the converted script, but not part of the job's command. |
| Trailing `code # comment` | No | Only the `code` before ` #` | Everything from the first space-`#` to end of line is stripped. |
| Heredoc `: << EOF … EOF` | No | Mangled | **Not supported** — flattened into the command and breaks the block. |
| **Blank line** | **Yes** | n/a | The intended, and only, way to close a block. |

**Example 1 — full-line comments between two commands**

```bash
#@1,0,job1
script1.sh
# comment A
# comment B
script2.sh
              # ← blank line ends the block
```

> **Does `script2.sh` run in job1? Yes.** Full-line comments don't end the block, so both commands are joined into a single job command (`script1.sh; script2.sh`) and run in job1. Comments A and B are copied into the converted script but are not part of the job.

**Example 2 — trailing comments on the command lines**

```bash
#@1,0,job1
script1.sh # comment A
script2.sh # comment B
              # ← blank line ends the block
```

> **Do both scripts run in job1? Yes.** Each ` # comment` is stripped from its line; the code before it (`script1.sh`, then `script2.sh`) is kept, so both run in job1.

**Example 3 — a heredoc used as a comment**

```bash
#@1,0,job1
script1.sh
: << EOF
a heredoc comment
EOF
script2.sh
```

> **Does `script2.sh` run in job1? No — the block is broken.** The parser is line-based and has no concept of heredocs. It joins every line with `;` and collapses whitespace into one line, producing roughly:
> ```
> script1.sh; : << EOF; a heredoc comment; EOF; script2.sh;
> ```
> With no real newlines, the `<< EOF` heredoc swallows everything after it — including `script2.sh` — as its body, so `script2.sh` never executes. **Don't use heredocs (or heredoc-style comments) inside a `#@` block; use `#` comments instead.**

> [!WARNING]
> Inline-comment stripping is **textual, not shell-aware**: everything from the first space-`#` (` #`) to end of line is removed. So avoid a literal space-`#` inside your command — *even inside quotes* — or it will be truncated. For example `sed 's/ #/x/'` gets cut down to `sed 's/`. A `#` with no space before it (e.g. `grep '#'`) is safe.


## How runAsPipeline runs: the two phases

Understanding **when** each line runs is critical when adapting scripts to use `runAsPipeline`. Scripts are parsed in two phases.

```
                 ┌─────────────────────── PHASE 1: CONVERT (once, on submit host) ──────────────────────┐
  your_script.sh │ read top → bottom:                                                                    │
                 │   • for/while/#loopStart  → remember the loop variable (becomes part of each job flag) │
                 │   • #@ marker + command    → turn into an  ssbatch --wrap "command"  call              │
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
> **Fix:** move any logic that depends on a step's results into a *later* `#@` step, so it runs on a node after the upstream job finishes.

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

The same submit-scope variable is substituted into both commands at submission, so both refer to the identical literal path — while the data itself flows through the filesystem, gated by the dependency. This is exactly what the proseq example does:

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

At submission `$mapInputR1` expands to the same path in step 1's and step 3's commands. Step 3's marker resolves `mapInputR1.mapInputR2` (its inputs) and `refIndex` (its reference) the same way. Step 3 depends on steps 1 and 2 (`#@3,1.2,...`), so it is held until they finish — only then does the file exist, get measured, and step 3 get sized and released.

> [!NOTE]
> Listing a **produced** file as a downstream step's `inputs` is safe even though it doesn't exist at submission time. Dependent steps are submitted **held**, and their resource estimation is deferred until the upstream step completes and the file exists.

### Recommended convention: name paths in submit scope, next to the step

Assign each step's input/output paths as plain variables **immediately above the `#@` marker that first uses them**, as the example does. This gives one source of truth per path, keeps producer and consumer in agreement automatically, and reads top-to-bottom as *"here is the file, here is the job that makes it, here is the job that uses it."*

Hoisting every variable to the very top of the script also works for *static* paths, but it reads worse and is fragile inside loops — prefer define-near-use.

> [!WARNING]
> **Submit-scope variables persist across steps, branches, and loop iterations.** A value set for one step is still set when a later step is submitted, so a *missing* assignment silently reuses a stale value. Two habits prevent this:
> - Assign the variable on **every branch** that leads to a step using it. (The proseq `if/else` sets `mapInputR1`/`mapInputR2` in *both* branches — do the same.)
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

`runAsPipeline` preserves loop structure but extracts the `#@` blocks inside. The **loop variable becomes part of each job's flag**, so per-iteration jobs get distinct names (e.g. `1.0.findNumber.1`, `1.0.findNumber.2`, …).

**`for` loops** work directly — the variable right after `for` is detected automatically:

```bash
for file in `ls someFolder`; do
    #@1,0,process,,file
    process.sh $file
done
# each iteration submits a job flagged ...process.$file
```

**`while` loops need a hint.** A `while` header doesn't name its loop variable in a position the parser can read, so declare it with `#loopStart:VAR` on the line above:

```bash
#loopStart:f1
while read -r f1 f2 f3 f4; do
    #@1,0,process,,f1
    process.sh "$f1"
done < samples.txt
```

## Tutorial

> **Note:** these tutorial scripts now live in `example-scripts/bashScript-tutorial/`. Run the commands below from that directory (or copy the scripts into your working directory).

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

- `#@1,0,findNumber,,input,sbatch …` — step **1**, depends on **nothing** (`0`), program **findNumber**, **no** reference, input is `$input`, with the given sbatch options. Because it's inside the `for` loop, it submits five jobs, one per `$i`.
- `#@2,1,mergeNumber,,,sbatch …` — step **2**, depends on **step 1**, program **mergeNumber**, no reference, no input. Slurm holds it until all five step-1 jobs finish.

**Dry run** (adds `--mode dryrun`, so nothing is submitted — you just see the plan and fake IDs):

```bash
runAsPipeline --script "bashScriptV2.sh 123" --sbatch-options "sbatch -p short -t 10:0 -c 1" --tmp useTmp --mode dryrun
```

**Real run** (default mode, no `--mode` needed):

```bash
runAsPipeline --script "bashScriptV2.sh 1234" --sbatch-options "sbatch -p short -t 10:0 -c 1" --tmp useTmp
```

Abbreviated output:

```text
runAsPipeline run date: 2024-04-28_16-03-36
Running: .../bin/runAsPipeline --script "bashScriptV2.sh 1234" --sbatch-options "sbatch -p short -t 10:0 -c 1" --tmp useTmp
===========
Stage 1: Processing Pipeline
    Converting pipeline to execution script (.../slurmPipeLine.<md5>.sh)
==========
Stage 2: Submitting jobs
---------------------------------------------------------
step: 1, depends on: 0, job name: findNumber, flag: 1.0.findNumber.1
Submitted batch job 69308
step: 1, depends on: 0, job name: findNumber, flag: 1.0.findNumber.2
Submitted batch job 69309
...
step: 2, depends on: 1, job name: mergeNumber, flag: 2.1.mergeNumber
Submitted batch job 69313

All submitted jobs:
job_id       depend_on                      job_flag          program     reference  inputs
69308        null                           1.0.findNumber.1  findNumber  none       numbers1.txt
69309        null                           1.0.findNumber.2  findNumber  none       numbers2.txt
...
69313        69308:69309:69310:69311:69312  2.1.mergeNumber   mergeNumber none       none
---------------------------------------------------------
```

> [!NOTE]
> Steps without their own sbatch options use the command-line default (here `-t 10:0`). In the script above both steps set `-t 50:0` themselves, so they override the default. This is how you mix a global default with per-step overrides.

After it runs:

```bash
ls -l smartSlurmLog   # per-step .sh (job scripts), .out (logs), .success / .failed flags
checkRun              # interactive status + log browser
cancelAllJobs         # cancel running/pending jobs from this directory
```

## checkRun: monitor and debug

`checkRun` is an interactive, three-level browser for a pipeline's status and logs. **Run it from the directory where you launched `runAsPipeline`** (it needs the `.smartSlurm.log` file there).

```bash
checkRun
```

### Level 1 — pick a run

Lists each log folder/file with its job count; dry runs are flagged.

| Key | Action |
|-----|--------|
| *number* | Open that run |
| `r` | Reload (refresh the list and job states) |
| `h` | Show/hide runs that submitted no jobs |
| `q` | Quit |

### Level 2 — job status table

Shows every job from `allJobs.txt`, each with a colored status label:

| Label | Meaning |
|-------|---------|
| `Done` | finished successfully (`.success` exists) |
| `Fail` | finished with failure (`.failed` exists) |
| `Runn` | currently running |
| `Pend` | pending in the queue |
| `Requ` | was requeued (e.g. after OOM/OOT) |
| `Unkn` | no success/failure flag and not in the queue — often killed or a node failure |

| Key | Action |
|-----|--------|
| *number* | Open that job's log files (Level 3) |
| `w` | Render the dependency DAG as an image (needs graphviz; see below) |
| `p` | Show/hide pending jobs |
| `q` | Back to Level 1 |
| `qq` | Quit |

### Level 3 — pick a log file

Lists the files for the selected job, tagged in plain language:

| Tag | File | Contents |
|-----|------|----------|
| `out` | `<flag>.out` | **stdout/stderr — start here to see the actual error** |
| `sh` | `<flag>.sh` | the generated Slurm script for the step |
| `adjust` | `<flag>.adjust` | resource-adjustment log |
| `success` / `failed` | flag files | status only, no contents |

| Key | Action |
|-----|--------|
| *number* | Open the file in `less` |
| `q` | Back to Level 2 |
| `qq` | Quit |

**Typical debugging flow**

```text
checkRun
  → [Level 1] pick your run
  → [Level 2] find the row labeled  Fail  or  Unkn
  → [Level 3] open its  out  log to read the error
```

**Workflow chart (`w`)** requires graphviz in the `smartSlurmEnv` conda environment (see [Installation](#installation)):

```bash
module load conda/miniforge3/24.11.3-0
conda activate smartSlurmEnv
```

`checkRun` then generates and displays a DAG of the pipeline's jobs.

> [!NOTE]
> To keep things fast, `checkRun` caches the `squeue` result for ~2 minutes and a generated DAG for ~10 minutes. Use `r` at Level 1 to force a refresh.

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

Add `noSuccEmail` (failures only) or `noEmail` (none) with `--email`:
```bash
runAsPipeline --script "bashScriptV2.sh 123" --sbatch-options "sbatch -p short -t 10:0 -c 1" --tmp useTmp --email noSuccEmail
```
</details>

<details>
<summary><b>Can I drop the command-line sbatch options?</b></summary>

Yes, **if every step sets its own** `sbatchOptions`. Then omit `--sbatch-options`:
```bash
runAsPipeline --script "bashScriptV2.sh 123" --tmp useTmp
```
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

`backup your config.txt before upgrade`:

```bash
cd ~/SmartSlurm
cp ~/.smartSlurm/config/config.txt ~/.smartSlurm/config/config.txt.backup
git pull

# the modify the new config.txt manually
```
