# SmartSlurm example scripts

Small, self-contained rAP pipelines that double as **best-practice templates** for
converting a pipeline to SmartSlurm and as **quick smoke tests** on the cluster. Every
step uses only common shell tools (`tr`, `sed`, `awk`, `sort`, `wc`, `cut`, `printf`),
asks for trivial resources, and produces output a human can verify by eye.

## How to run any of them

```bash
# add the SmartSlurm bin/ to PATH first (so rAP finds its sibling tools)
export PATH=/path/to/SmartSlurm/bin:$PATH

# preview the conversion without submitting anything:
runAsPipeline --script "demo_linear_fan_branch.sh demo_manifest.txt" --tmp noTmp --mode dryrun
# check block structure only (no submit, exit 1 on a problem):
runAsPipeline --script "demo_linear_fan_branch.sh demo_manifest.txt" --tmp noTmp --lint
# real run:
runAsPipeline --script "demo_linear_fan_branch.sh demo_manifest.txt" \
              --sbatch-options "sbatch -p short -c 1 --mem 100M -t 5:00" --tmp noTmp
```

Run each example in its **own directory** — pipelines share step names, so running two in
one directory would mix their `smartSlurmLog/` job-flag namespaces.

## What's here

| Script | Style | Patterns it teaches | Expected result |
|---|---|---|---|
| `demo_linear_fan_branch.sh` | `#@end` (scripts-on-disk) | the full progression: a **linear** chain (normalize→tokenize), a **fan-out/fan-in** pair (split into 3 chunks → count each → merge), and an **if/else** that submits one of two alternative terminal steps | per-sample report; `total_words` equals the token count; long inputs → `DETAILED`, short → `brief` |
| `phrase_pipeline_end.sh` | `#@end` | a minimal linear pipeline plus an experiment-level **fan-in** `summary` — the smallest "start here" `#@end` example | `summary.tsv` with `longest=brown/dozen/vexingly/quickly` |
| `phrase_pipeline_legacy.sh` | legacy blank-line (`--wrap`) | the older single-line-body style, and its constraint (avoid quoted-whitespace args — the `--wrap` flattener collapses them) | `all_counts.txt` (24 lines) |
| `bashScript-tutorial/` | mixed | the original "**annotate your pipeline**" tutorial: `V1` plain → `V2` annotated (the example named in `runAsPipeline` usage), plus edge-case variants. See its own README. | — |
| `resource-estimation/` | `#@end` | **jobRecords / resource estimation**: a `load` step whose memory + wall-time scale with input size, for testing curve-fitting and OOM/OOT auto-rerun. See its own README. | `load.<s>.txt` per sample; jobRecords across 4 sizes |

`demo_*` and `phrase_*` each use their own manifest. The `bashScript-tutorial/` scripts use
`createNumberFiles.sh` + `findNumber.sh` from `SmartSlurm/bin/`.

## `#@end` vs legacy

Prefer **`#@end`** blocks for new pipelines: the body is written to a per-job script
verbatim, so `awk`, quotes, multi-line bodies, and command substitution all survive. The
legacy blank-line style flattens each block into a single `--wrap` string, which mangles
quoted whitespace and multi-line bodies — `phrase_pipeline_legacy.sh` is included to show
that style and its limits, not as the recommended default.

## Conventions these examples follow (and why)

- **Marker format:** `#@ step , dependsOn , jobName , reference , inputs , sbatchOpts`.
- **Loop var is submit-scope:** `$s` from `for s in $samples` is baked into each job as it
  is submitted; a variable *assigned inside* a block expands on the node.
- **Escape `$` meant for a tool:** awk field refs must be written `\$1`, or the pipeline's
  own shell argument is substituted instead.
- **Fan-out flags are unique:** rAP includes *every* enclosing `for`-loop variable in the
  job flag (`…countchunk.$s.$c`), so nested chunk jobs never collide.
- **Branch conditions must be knowable at submit time** (from the input/manifest), never
  from a job's output — the driver submits every job before any of them runs.
- **if/else branches use distinct step numbers** — reusing a step number across branches
  makes rAP drop the second block.
- **No `for` loops inside a block body** — a loop variable assigned inside a block is not
  recognized as node-scope and would bake empty; use explicit filenames.

## `resource-estimation/`

A tidied, self-contained replacement for the old `useMemTime*` scripts: `sizedLoad.sh`
burns memory and wall-time in proportion to input rows (clean `awk`-array allocation,
env-tunable slopes/noise, and a deterministic overshoot knob), driven by
`jobrecords_pipeline.sh` over graded input sizes. It seeds the jobRecord fit and, with
`FAIL_PROB=1`, forces an OOM/OOT to exercise resource auto-increment and rerun. See
`resource-estimation/README.md` for the protocol.
