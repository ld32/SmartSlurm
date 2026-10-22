# bashScript tutorial — the canonical "annotate your pipeline" example

These are SmartSlurm's original teaching pipelines, migrated here from `SmartSlurm/scripts/`.
They show the core move: take a plain bash pipeline and add `#@` markers so rAP submits its
steps as dependency-aware cluster jobs. `bashScriptV2.sh` is the example referenced in
`runAsPipeline`'s own `--help`/usage text.

## The teaching arc

| Script | What it is |
|---|---|
| `bashScriptV1.sh` | The **plain** pipeline, *before* any annotation — a loop that runs `findNumber.sh` over `numbers1..5.txt`, then `cat`s the results. This is "your pipeline as you'd normally write it." |
| `bashScriptV2.sh` | The **same pipeline, annotated** — `#@1,0,findNumber,...` on the per-sample step and `#@2,1,mergeNumber,...` on the merge. The canonical before→after comparison. |
| `bashScriptV2x.sh` | A variant with a **nested loop** and a `sleep`, plus commented-out alternative markers — a scratchpad for experimenting with loop shapes. Less a clean example than a mechanics probe. |
| `bashScriptV3.sh` | Deliberately declares an **input that does not exist at submit time** (`all$number.txt`, produced only by the merge) to show how rAP handles a not-yet-present input. |

## Dependencies (already in `SmartSlurm/bin/`, on PATH)

These call two core helpers that ship with SmartSlurm — they are **not** duplicated here:

- `findNumber.sh` — the toy "tool": greps a number from a file. (The `bin/` version also
  allocates memory and sleeps in proportion to input size — see the note below.)
- `createNumberFiles.sh` — generates the graded inputs `numbers1..10.txt`, where
  `numbers1` is the base size and `numbers2..5` are 2×…5× larger. **Run it first.**

## Run it

```bash
export PATH=/path/to/SmartSlurm/bin:$PATH
createNumberFiles.sh                      # makes numbers1..10.txt (varied sizes)
runAsPipeline --script "bashScriptV2.sh 777" \
              --sbatch-options "sbatch -p short -c 1 --mem 2G -t 50:0" --tmp noTmp
```

## Note: these double as a resource-estimation testbed

Because `createNumberFiles.sh` produces **differently-sized** inputs and the `bin/`
`findNumber.sh` consumes memory and wall-time **in proportion to input size**, running
`bashScriptV2.sh` across `numbers1..5.txt` already generates jobRecords at several input
sizes — enough for rAP to fit a mem/time curve. This is the existing basis for the
resource-estimation / auto-rerun testing we were about to build from scratch; see the
top-level README's "Pending" note.
