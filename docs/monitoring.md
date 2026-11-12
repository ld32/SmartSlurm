[&larr; SmartSlurm README](../README.md) &middot; [Docs index](README.md)

# Monitoring and diagnosing runs

## The status grid (checkRun)

`checkRun` is an interactive monitor. **Run it from the directory where you launched `runAsPipeline`** (it reads `.smartSlurm.log` and `smartSlurmLog/` there):

```bash
checkRun
```

It opens on a **status grid** — samples down the side, steps across the top, one glyph per cell:

```text
Found 24 jobs

Status:  ✓ 12 done   ✗ 0 error   ▸ 0 running   · 12 pending

════════════════════════════════════════════════════════════
Showing run: /n/scratch/users/.../myProject
════════════════════════════════════════════════════════════

   row  sample      normali tokeniz  split  countch  merge  detaile  brief
                       1       2       3       4       5       6       7
   a    alpha      │   +   │   +   │   +   │   @   │   +   │   +   │   ·   │
   b    bravo      │   +   │   +   │   +   │   @   │   +   │   ·   │   +   │
   ...
```

A band at the top names the run currently on screen (its full path). If other SmartSlurm runs exist in this directory, it also shows `--> N other run(s) available — press o to switch`.

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

**Addressing a cell.** A job address is a **row letter + step number** — `b2` is row b, step 2; experiment-level steps are shown with a capital `E` (e.g. `E3`). Both orders parse (`b2` or `2b`, `E3` or `3E`). At the prompt:

| Type | Action |
|------|--------|
| `b2` | view that job's log (the `.out`) |
| `d b2` | **diagnose** — why it failed / what to do (runs [`diagnose.sh`](#diagnosing-a-job-diagnosesh)) |
| `a b2` | ask the AI helper about the job |
| `r b2` | view the job's raw files (`.out` / `.err` / `.sh` …) |
| `o` | browse other SmartSlurm log folders in this directory |
| `l` | list view (the classic per-job browser; all its capabilities are preserved) |
| `g` | workflow chart (DAG image; needs graphviz — see [Installation and setup](installation.md#optional-graphviz-for-workflow-charts)) |
| `e` | **resource efficiency** — per-sample time/mem used vs allocated, with drill-down (see [Resource efficiency](#resource-efficiency)) |
| `p` | show/hide pending jobs |
| `Enter` | refresh (re-query `squeue`, redraw) |
| `q` / `qq` | back up a level / quit |

> [!NOTE]
> `checkRun` is **interactive only** — it takes no arguments (other than `-h`).
> There is no `checkRun grid` / `checkRun why` / `checkRun --json`.

**For scripting**, call the engines `checkRun` drives from its UI directly:

```bash
diagnose.sh <flag> --log smartSlurmLog --json   # one job's machine-readable record
checkRunGrid --render  --log smartSlurmLog       # print the status grid and exit
checkRunGrid --resolve b2 --log smartSlurmLog    # address -> job flag(s), one per line
checkRunGrid b2 why    --log smartSlurmLog        # diagnose the job(s) at address b2
```

`checkRunGrid` is the render/resolve engine and [`diagnose.sh`](#diagnosing-a-job-diagnosesh) the per-job analyzer; `checkRun` drives both from its keyboard UI.

**Typical debugging flow:** open `checkRun`, find a cell showing `X` (or `x`/`?`), type `d <addr>` to get the cause and suggested action, then `<addr>` to read the log.

## Diagnosing a job (diagnose.sh)

`diagnose.sh` answers "what happened, why, and what to do" for a single job. It is authoritative about *whether* a job succeeded (`.success` + sacct) and heuristic about *why* (an ordered pattern table over the tool's output region), and it never asserts a cause without an evidence line.

```bash
diagnose.sh <flag> [--log DIR] [--json]
# e.g.
diagnose.sh 1.0.normalize.alpha
diagnose.sh 1.0.normalize.alpha --json     # machine-readable record
```

The human report's first line is `<flag>  [TAG]  state=<STATE>  action=<ACTION>`, where `TAG` is `OK`, `RETRYABLE`, `DOOMED`, or the action name (`INSPECT`, `CHECK`, `CHECK_UPSTREAM`, `CHECK_CLUSTER`, `CHECK_INPUT`, …). Below it come a plain-language summary, the evidence line it keyed on (`evidence (<log>:<lineno>): …`, when one was found), and the full-log path (for failed jobs).

The `--json` record carries these fields — consume it by name:

| Field | Meaning |
|-------|---------|
| `flag` | the job flag diagnosed |
| `failed` | `true`/`false` — authoritative pass/fail verdict |
| `state` | `COMPLETED`, `OOM`, `OOT`, `OOMOOT`, `NODE_FAIL`, `CANCELLED`, `FAILED`, `COMPLETED_NOSUCCESS`, `ORPHANED`, `UNCONFIRMED` |
| `category` | cause class, e.g. `resource:oom`, `environment:missing-command`, `tool:crash`, `orphaned:upstream-failed`, `infrastructure:uncertain` |
| `action` | advice: `NONE`, `RETRYABLE`, `RETRYABLE_EXHAUSTED`, `DOOMED`, `INSPECT`, `CHECK`, `CHECK_UPSTREAM`, `CHECK_CLUSTER`, `CHECK_INPUT` |
| `retryable` | `true`/`false` |
| `requeued` | `true`/`false` — runAsPipeline already auto-bumped resources and resubmitted |
| `evidence` | the exact log line the cause keyed on (may be empty) |
| `evidence_line` | its line number in the `.out` (or `null`) |
| `region` | where the scan looked: `tool-region` or `whole-log(legacy)` |
| `log` | path to the job's `.out` |
| `summary` | the one-line plain-language explanation |

## Resource efficiency

Every job writes one row to a per-run **`smartSlurmLog/efficiency.tsv`** (flag, sample, step, status, and the memory/time it was **allocated vs actually used**). Press **`e`** in `checkRun` to explore it as `htop`-style bars — so you can see where you over- or under-reserved and tighten future runs.

Bars are **percent-based** (fixed width, so absolute scale never distorts the picture); the `eff | used | alloc` columns follow each bar. **Time bars are cyan, memory bars magenta** (distinct and deliberately not green, so the time/mem colors never read as a red/green pass–fail scale; plain under `--no-color` or when piped). A row whose efficiency is at or below `effLowPct` (default **25%**, set in [`config.txt`](configuration.md#configtxt)) turns its bar **red** — a waste alert for allocations that went mostly unused, worth trimming. For a sample, **time is summed across its steps** (step-time is additive) and **memory is the peak step** (steps run sequentially, so a sample never holds the sum at once).

`e` opens the **overview**: an `Overall` total (also shown as bars), then one row per sample, then — if the pipeline has any **experiment-level steps** (jobs that span all samples, e.g. a merge or a MultiQC summary, so they belong to no single sample) — an `Experiment-level (across all samples)` section with one row each, then the sample/step menus. Steps are listed in ascending numeric order; sample drill-in letters are **UPPERCASE** (lowercase keys are commands).

```text
Resource efficiency -- overview
                                                         eff    used   alloc
Overall        (time) [|||||||||||||||||||||.........]   69%   32.9h   48.0h
               (mem)  [||||||||||||||||||||||||......]   82%   78.3G   96.0G

liver (A)      (time) [||||||||||||||||||||||||......]   79%   12.6h   16.0h
               (mem)  [|||||||||||||||||||||||||||||.]   95%   30.5G   32.0G
brain (B)      (time) [|||||||||||||||...............]   50%    8.0h   16.0h
               (mem)  [||||||||||||||||||............]   61%   19.5G   32.0G
heart (C)      (time) [|||||||||||||||||||||||.......]   77%   12.3h   16.0h
               (mem)  [|||||||||||||||||||||||||||...]   89%   28.3G   32.0G

Samples:  liver (A),  brain (B),  heart (C)
Steps:    cutadapt (1),  bowtie (2),  count (3)
```

From there you drill down, and each view opens instantly (re-rendered from the same `efficiency.tsv`):

| Key | View |
|-----|------|
| a **sample letter** (`A`) | that sample's **per-step** breakdown (total row + one row per step) |
| a **step number** (`2`) | that step **across all samples** (total row + one row per sample) |
| `t` / `m` | toggle the time / memory rows |
| `o` | back to the overview |
| `q` | leave the efficiency view (from the overview, back to `checkRun`) |

```text
# after typing 2  (bowtie across all samples)
Step view: showing bowtie (2)
                                                         eff    used   alloc
bowtie         (time)                                    73%   22.0h   30.0h
               (mem)                                     82%   78.3G   96.0G

liver (A)      (time) [|||||||||||||||||||||||||||...]   90%    9.0h   10.0h
               (mem)  [|||||||||||||||||||||||||||||.]   95%   30.5G   32.0G
brain (B)      (time) [|||||||||||||||...............]   50%    5.0h   10.0h
               (mem)  [||||||||||||||||||............]   61%   19.5G   32.0G
heart (C)      (time) [||||||||||||||||||||||||......]   80%    8.0h   10.0h
               (mem)  [|||||||||||||||||||||||||||...]   89%   28.3G   32.0G
```

> [!NOTE]
> In aggregate rows (`Overall`, and the sample/step totals) **time** is the true sum and **memory** is the **sum of the per-sample (or per-step) peaks** — a reservation-efficiency ratio, not a simultaneous peak.

For scripting, the engine renders each view non-interactively:

```bash
checkRunGrid --efficiency                      --log smartSlurmLog   # overview
checkRunGrid --efficiency --eff-view sample --eff-select A --log smartSlurmLog
checkRunGrid --efficiency --eff-view step   --eff-select 2 --log smartSlurmLog
checkRunGrid --efficiency --eff-show time   --log smartSlurmLog      # one dimension
```

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

---


---

[&larr; SmartSlurm README](../README.md) &middot; [Docs index](README.md)
