[&larr; SmartSlurm README](../README.md) &middot; [Docs index](README.md)

# AI assistance (optional)

SmartSlurm can call a language model to explain failures and sanity-check its own estimates. **This is entirely optional.** Every AI call is guarded, so if AI is not set up — or is turned off, or the call fails or times out — each script behaves exactly as it did before: same emails, same estimates, same menus.

## Setup

AI calls go through a small setuid helper binary, `ai_call`, which holds the API token so it never appears in a script, a process list, or a log.

```bash
# Built and installed by an admin, typically at:
~/.smartSlurm/bin/ai_call          # setuid (mode 4755), compiled from ai_call.c

# Access is granted per user by the helper's owner:
echo yourUserName >> ~/.smartSlurm/allowed_users.txt
```

`ai_call` is found automatically at `~/.smartSlurm/bin/ai_call`, then next to the SmartSlurm scripts. Set `smartSlurmAiCall` in [`config.txt`](configuration.md) to point somewhere else.

If it is missing, or your username is not in `allowed_users.txt`, SmartSlurm says so once and carries on without AI.

## What AI does, and when it runs

| Where | When it runs | What you get |
|---|---|---|
| `ssbatch` | automatically, when `-P`, `-R` or `-I` is missing — **once per job shape**, then cached | Program name, reference and input files for the estimate. See [Automatic program and input detection](#automatic-program-and-input-detection) |
| `cleanUp.sh` (every job) | automatically, when a job ends **not** COMPLETED (Canceled excluded) | Failure analysis at the top of the job email, also saved to `smartSlurmLog/ai.<jobid>.txt` |
| `estimateResource.sh` | automatically, **only when a resource curve is rebuilt** and there are ≥ `smartSlurmAiMinRecords` records | Review of the fit: outliers, non-linearity, regime changes. Saved to `<recordDir>/stats/<program>.<ref>.aiReview.txt` and copied into the step's `.out` log |
| [`checkRun`](monitoring.md) | on demand: `1a`, `2a`, … or `?` | Per-job analysis (tailored to whether the job is done, failed, running or pending), or a summary of the whole run |
| [`reviewJobRecords.py`](tools.md#reviewjobrecordspy) | on demand: `?` | Outliers and trends in a program's records, ending with a `SUGGESTED DELETIONS:` line |

`calculateMemTime.sh` runs on **every** submission, so it deliberately makes **no** AI call — it only prints a review that `estimateResource.sh` produced earlier.

> [!NOTE]
> AI never sets memory or time itself, deletes a record, cancels a job, or blocks a submission. The one thing it decides is *which* program, reference and input files a job has (only when you did not say), and the estimate is then computed from your job history as usual. Everything else it does is text for you to read.

## Automatic program and input detection

Resource estimates are keyed on a program name and scale with input size, so `ssbatch` needs to know both. You can say so with `-P program -R reference -I input1,input2`. When you leave any of them out, `ssbatch` asks AI to read the `--wrap` command or the job script and fill in the missing ones:

```bash
ssbatch -p short -t 2:0:0 --mem 8G -c 4 --wrap "bwa mem -t 4 ref/hg38.fa A_R1.fq.gz A_R2.fq.gz > A.sam"
#AI detect (new AI query, saved as 1ce5f821df15798c): program=bwa-mem-4core ref=hg38 inputs=A_R1.fq.gz,A_R2.fq.gz
```

- **AI names, code checks.** AI names only the executable; the program name is built from the command itself (tool + subcommand such as `bwa-mem`, or the script name for `python`/`Rscript`/`java -jar` jobs, plus the core count). A tool that is not in the job is rejected, and only input paths that exist are kept. Input size is still measured with `du`, and memory/time still come from your job records.
- **One AI call per job shape.** The answer is cached in `~/.smartSlurm/aiCache/<key>.json`. The key ignores file names, option values and numbers, so every sample of the same step reuses it with no AI call, and finds its own input files locally. The log line says `(cache …, no AI call)` or `(new AI query …)`.
- **Your flags always win.** Anything given with `-P`, `-R` or `-I` is used as is.
- **Fixing a wrong guess.** Edit that step's JSON file in `aiCache/` (the log line names it); every later job of that shape uses your fix. Delete the file to ask AI again.
- **Not found:** inputs that only exist at run time (built from a sample sheet, a loop or `$(…)`, or written by an earlier step), inputs given as a prefix (`plink --bfile geno`), and files a program opens without naming them. Such jobs fall back to the per-program estimate; pass `-I` if you know the input.
- **Preview first.** `ssbatch -dryrun …` shows the detection and the estimate without submitting.

### Evaluating detection and estimates after test runs

**1. Was the detection right?** Each job's log in `smartSlurmLog/` keeps the detection line. From the run folder:

```bash
grep -h "^#AI detect" smartSlurmLog/*.out | sort | uniq -c
```

Check that:
- the **program name** is the same for every sample of one step, and different for different steps (`bwa-mem-4core` vs `samtools-sort-4core`);
- the **inputs** are the files whose size drives the cost, not references, indexes, configs or outputs;
- **`inputs=none`** appears only for jobs that really have no input.

See what AI decided for each job shape, and how often it was asked:

```bash
ls ~/.smartSlurm/aiCache/ | wc -l       # AI calls ever made (one per job shape)
grep -h '"example"\|"program"\|"input_flags"' ~/.smartSlurm/aiCache/*.json
grep -c "AI detect (cache" .smartSlurm.log; grep -c "AI detect (new" .smartSlurm.log   # cache hits vs AI calls
```

The first sample of each step should be `new` and the rest `cache`. A step that keeps showing `new` has a per-sample value the cache key does not mask; pass `-P/-I` for it, or report the two commands. To correct a guess, edit its JSON file or pass `-P/-R/-I`.

**2. Did the estimates improve?** A step needs at least 3 COMPLETED jobs under the same name before it is fitted, so run 4–5 samples per step. Then summarise [`jobRecord.txt`](core-concepts.md#jobrecords). The commands below find each column by its header name, so they keep working if columns are added or reordered:

```bash
f=~/.smartSlurm/jobRecord.txt
awk -F, 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; next}
  { notes = ("notes" in c) ? $c["notes"] : ""
    by = (notes ~ /^AI-found/) ? "ai" : "user"            # who named program/inputs
    k = $c["step"] " " $c["ref"] " " by
    n[k]++; if ($c["status"] != "COMPLETED") bad[k]++
    if ($c["inputSize"] + 0 == 0) zero[k]++
    if ($c["memAllocated"]  > 0) mu[k] += $c["memUsed"]  / $c["memAllocated"]
    if ($c["timeAllocated"] > 0) tu[k] += $c["timeUsed"] / $c["timeAllocated"] }
  END { printf "%-40s %5s %5s %6s %8s %8s\n", "step ref by", "jobs", "fail", "size0", "mem-use", "time-use"
        for (k in n) printf "%-40s %5d %5d %6d %7.0f%% %7.0f%%\n", k, n[k], bad[k], zero[k], 100*mu[k]/n[k], 100*tu[k]/n[k] }' $f | { IFS= read -r h; echo "$h"; sort; }
```

| Column | Good | Bad sign |
|---|---|---|
| `by` | `ai` rows and `user` rows of the same step look alike | `ai` rows clearly worse: detection named the step or its inputs differently from what you would have |
| `fail` | near 0 | out-of-memory/time jobs: each was resubmitted with more resources, costing a retry |
| `size0` (jobs with input size 0) | 0 for steps that read input | equal to `jobs`: inputs were not found, so the estimate cannot scale with size |
| `mem-use`, `time-use` (share of the request actually used) | about 60–90% | ~10%: still over-requesting, usually too few records or a poor fit |

`by` comes from the [`notes`](core-concepts.md#jobrecords) column: a row is `ai` when AI detection filled in its program, reference or inputs. To list just those rows with what AI found:

```bash
awk -F, 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; next}
  ("notes" in c) && $c["notes"] ~ /^AI-found/ {print $c["jobID"], $c["step"], $c["ref"], $c["sample"], $c["inputSize"]"K", $c["status"], $c["notes"]}' ~/.smartSlurm/jobRecord.txt
```

Look closer with [`reviewJobRecords.py`](tools.md#reviewjobrecordspy) (records per step, with an optional AI review), [`checkRun`](monitoring.md) → `e` (efficiency of this run), and the fitted plots `~/.smartSlurm/stats/<pipeline>.<step>.<ref>.mem.png` / `.time.png`.

For a fair comparison, run the same few samples twice, once with `-P/-R/-I` and once without, then compare the `user` and `ai` lines of the summary. If AI found the same inputs you would have given, the step names match and the estimates are identical.

## Settings

All in [`config.txt`](configuration.md) (your `~/.smartSlurm/config/config.txt` copy wins):

```bash
export smartSlurmAiEnabled=yes        # no = turn AI off completely
export smartSlurmAiCall=              # empty = auto-detect ai_call
export smartSlurmAiTimeout=90         # s, default cap for any AI call
export smartSlurmAiDetect=yes         # no = never ask AI for program/reference/inputs

export smartSlurmAiMinRecords=10      # records needed before a fitted curve is reviewed
export smartSlurmAiReviewTimeout=30   # s, curve review — short, a submission is waiting

export smartSlurmAiRecordReviewTimeout=180  # s, reviewJobRecords.py — you asked for it, so it may take longer
```

> [!TIP]
> Because `config.txt` is *sourced* (its `export`s override your environment), turn AI off by editing the config file — not with `smartSlurmAiEnabled=no ssbatch ...` on the command line.

Every timeout exists so AI can never delay a job or a submission: when it expires, the call is abandoned and the script continues.

## For developers

Two small libraries provide the guarded interface; use them rather than calling `ai_call` directly:

```bash
# bin/aiHelper.sh — shell
. "$(dirname $0)/aiHelper.sh"
aiEnabled || return            # false when AI is off, missing, or not permitted
echo "$prompt" | aiAsk         # prints the answer; on any failure prints nothing and returns non-zero
aiWhyNot                       # one-line reason, for a log or a menu
```

```python
# bin/aiHelper.py — python
import aiHelper
if aiHelper.ai_enabled():
    answer = aiHelper.ai_ask(prompt)   # str, or None on any failure
```

`aiAsk`/`ai_ask` never write a partial answer to stdout and never raise, so a caller can treat failure as simply "no AI today". Prompts are capped at 48 KB (`ai_call` rejects anything over 64 KB).

---
