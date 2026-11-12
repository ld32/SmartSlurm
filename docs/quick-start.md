[&larr; SmartSlurm README](../README.md) &middot; [Docs index](README.md)

# Quick start

A minimal `runAsPipeline` example. The tutorial scripts (`bashScriptV1.sh` …
`bashScriptV3.sh`) live in `example-scripts/bashScript-tutorial/`. Run the
commands below from that directory, or copy the scripts into your working
directory. They rely on `createNumberFiles.sh` and `findNumber.sh` (both in
`bin/`, on your `PATH`).

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

**Reading the markers** (full field reference: [The `#@` job block](core-concepts.md#the--job-block)):

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

Its Stage-2 output is the compact table described in [runAsPipeline](tools.md#runaspipeline) — one row per `$i`, five `findNumber` submissions plus the `mergeNumber` row. Add `-v` to see the full per-job "Submitted batch job N" firehose and the jobs table:

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

More worked examples live under `example-scripts/` — see
[Tools at a glance](tools.md#tools-at-a-glance) and the per-directory READMEs.

---


---

[&larr; SmartSlurm README](../README.md) &middot; [Docs index](README.md)
