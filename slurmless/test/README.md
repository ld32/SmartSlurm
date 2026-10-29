# SmartSlurm test harness (Slurm-free)

This directory lets you run **unmodified** `runAsPipeline` / `ssbatch` on a
machine with **no Slurm** (laptop, CI, web sandbox), while still executing real
job payloads — so pipeline conversion, dependency wiring, and failure attribution
can be tested for real. It is the regression net for the UX-overhaul work
(Module **M0** in the plan): "did we break rAP?"

```
test/
├── run_tests.sh          # dependency-free test runner (CI exit code)
├── lib/
│   └── harness.sh        # environment activation + assertion helpers
├── mockSlurm/            # stand-ins for sbatch/sacct/srun/… (see its README)
└── fixtures/             # tiny pipelines exercising specific behaviors
    ├── clean_2step/                 # baseline — MUST pass
    ├── fail_exit0_empty/            # "symptom ≠ cause" — north-star scenario
    └── known_broken_swallowed_done/ # the blank-line block foot-gun (M4)
```

## Quick start

```bash
# from the repo root
./test/run_tests.sh                 # pipeline + resource-model net (M0/M1)
./test/run_tests.sh clean           # only fixtures whose name matches "clean"

./test/checkrun_ui_tests.sh         # interactive-monitor UI net (M3: checkRun/checkRunGrid)
./test/checkrun_ui_tests.sh grid    # only UI blocks whose name matches "grid"

# interactive experimentation
source test/lib/harness.sh
ss_activate
ss_new_workspace /tmp/play
cp test/fixtures/clean_2step/* .
ss_run_pipeline "./pipeline.sh 777" "sbatch -p short -c 1 --mem 4G -t 50:0"
cat smartSlurmLog/allJobs.txt
```

Scratch state (job records, mock scheduler state, run logs) goes under
`$SS_SCRATCH` (default `/tmp/ss-scratch`); nothing touches a real `~/.smartSlurm`.

### UI regression net — `checkrun_ui_tests.sh`

A second suite covering the interactive monitor (M3): the `checkRunGrid` engine
(framed grid render, base-26 row labels, address/number resolve, gap/multi-dep
step handling, sample-label trimming + `--samples`) and the `checkRun` controller
(grid-as-home, `<verb> <addr>` dispatch, list/raw/switcher views, in-place redraw),
driven over a pseudo-terminal (needs `python3`). It shares `lib/harness.sh` and the
`clean_2step` fixture; run it from within a SmartSlurm checkout like `run_tests.sh`.

It pins the behavior of the checkRun/checkRunGrid UI rework (e.g. `--samples`,
base-26 labels). That rework has landed in `ld32/smartSlurm` `master`, so this
suite passes on current `master`; on an older SmartSlurm without it, the
UI-specific checks fail until it is present.

## How it works

`runAsPipeline` converts an annotated script into a driver that submits jobs via
`ssbatch`, which in turn calls the real `sbatch`. The harness interposes **shims**
on `$PATH` (ahead of the real Slurm binaries and the repo's `bin/`) so those calls
resolve to our fakes. This is the same `$PATH`-interposition trick SmartSlurm
itself uses for `sbatchBin`/`unExport` — not aliases (which don't survive into
child scripts), not passthroughs (there is no real Slurm here).

The mock `sbatch` **actually executes the generated job script**, so the real
`cleanUp.sh` trap and the real payload run, producing genuine exit codes,
`.out` logs, output files, and `.success`/`.failed` flags. The mock `sacct`
replays each job's recorded state so `cleanUp.sh`'s OOM/TIMEOUT/FAILED detection
behaves. See `mockSlurm/README.md` for the per-command contracts.

Because jobs run **synchronously in submission order**, `afterok` dependencies are
always satisfied by the time a dependent job runs, so `-d`/`-H` need no scheduler.

**Why the shims are staged to a neutral path.** `ssbatch` finds the real `sbatch`
with `which -a sbatch | grep -v SmartSlurm` — it deliberately skips any `sbatch`
whose **path contains the substring "SmartSlurm"** so it never picks up its own
wrapper. The mocks are committed under `test/mockSlurm/` (inside the repo, whose
path contains "SmartSlurm"), so at activation `ss_activate` copies them to
`$SS_SCRATCH/slurm-shims/` and puts *that* on `$PATH`. If you invoke a mock
`sbatch` directly for debugging, run it from the staged copy, not from
`test/mockSlurm/`, or `ssbatch` will filter it out and `realSbatch` will be empty
(symptom: `ssbatch: line …: --mail-type=FAIL: command not found`).

## Important finding: `#!/bin/sh` must be bash

SmartSlurm's scripts start with `#!/bin/sh` but use bash-only features
(`[[ … ]]`, arrays, `${x//…}`, `set -o pipefail`, …). On the **production RHEL
cluster `/bin/sh` is bash**, so this is fine there. On Debian/Ubuntu (including
this sandbox) `/bin/sh` is **dash**, and the scripts break in subtle ways — e.g.
`estimateResource.sh` silently returns nothing, so `ssbatch` aborts with
"did not find time limit"; and the generated driver dies with "Bad substitution".

`ss_activate` therefore repoints `/bin/sh` → `bash` when it can (it needs root;
the sandbox is disposable). **This is faithful emulation of production, not a
hack.** If you cannot repoint `/bin/sh` (no root), run on a system where `sh` is
already bash, or the `#!/bin/sh` scripts will misbehave.

## The fixtures

**`clean_2step` — baseline, must always pass.** Two parallel `findNumber` jobs
feed one `mergeNumber` job. Proves the whole path works: conversion → valid bash
driver → real payload execution → `.success` flags → correct `afterok` wiring in
`allJobs.txt`. If this breaks, rAP is broken.

**`fail_exit0_empty` — the scenario the overhaul exists to fix.** The producer
cannot read its input, writes an **empty** file, and **exits 0**; the consumer
then fails on the empty input. The harness asserts the *current* behavior, which
is wrong on purpose: the producer is marked `.success` and the **consumer** (the
symptom, not the cause) is marked `.failed`. After Module **M1** (producer
output-check) + **M2** (DAG-aware root-cause tracer), this fixture should be
updated to assert the producer is flagged as the root cause. It is the executable
definition of the problem.

**`known_broken_swallowed_done` — the M4 foot-gun.** A block whose body is
immediately followed by the loop's `done` (no blank line) makes the blank-line
terminated parser swallow `done` into the wrapped command, leaving the outer loop
unterminated — the converter emits a driver that is **invalid bash**. (The
shipped `scripts/bashScriptV2x.sh` trips exactly this.) The fixture asserts the
driver fails `bash -n`. Module **M4**'s explicit `#@ … #@end` closure fixes it;
until then this is the regression anchor.

## M1: resource capture, curve fitting, and jobRecords

The base harness runs tiny payloads, so the real `memCpuMonitor.sh`/`sacct`
capture records a near-zero peak and ~1 min — useless for exercising
`estimateResource.sh`'s curve fitting. The **resource-model layer** fixes this:
it synthesizes the mem/time a job of a given input size would *measure* on a real
cluster, feeds that to the **unmodified** `cleanUp.sh` → `jobRecord.txt` →
`estimateResource.sh` → `calculateMemTime.sh` loop, and lets us test the learning
loop end-to-end, Slurm-free.

**Nothing in `bin/` is modified.** The injection lives entirely in three shims
(`sbatch`, `sacct`, `memCpuMonitor.sh`) plus a model evaluator
(`mockSlurm/resourceModel.sh`). When no model is active the shims behave exactly
as before (functional blocks still pass, `memCpuMonitor.sh` execs the real tool).

### How it works

A **model** (`test/models/*.model`, `KEY=VALUE`) declares a ground-truth curve of
mem and time as a function of input size (KB): `linear` (`a*x+b`), `exp`
(`b*exp(a*x)`), `poly` (`a*x²+b*x+c`), or `const`, with optional reproducible
gaussian noise (a self-contained LCG — mawk's `srand` is not deterministic). At
submit time the mock `sbatch` reads the job's program/ref/inputSize and the
*allocated* mem/time from the generated `cleanUp.sh` trap line, evaluates the
model, and caches `<jobid>.resource`. Then:

- mock **`memCpuMonitor.sh`** writes a synthetic `job_<id>.memCPU.txt` whose line-1
  `START` is back-dated so `cleanUp.sh` computes the modeled minutes, and whose
  sample peak is the modeled memory;
- mock **`sacct`** reports the modeled `MaxRSS`;
- **emergent errors** fall out for free: if the (bad-fit) *allocation* is smaller
  than the modeled need, the job genuinely OOMs/OOTs — no explicit injection.

### API (in `lib/harness.sh`)

```bash
ss_use_model NAME            # activate test/models/NAME.model (or a path)
ss_clear_model
ss_run_sweep NAME "S1 S2 …" "sbatchOpts"   # submit the 1-step pipeline once per size (KB)
ss_fit PROGRAM REF           # force a fresh gnuplot fit; echo the mem.stat keys
ss_seed_record P R SIZE MEM TIME STATUS [n] [saccMem]   # append crafted rows
assert_near EXPECTED GOT TOL_FRACTION MSG
stat_val KEY FILE
```

`ss_activate` installs **gnuplot** if missing (the fit needs it); curve-fit blocks
skip cleanly when it is unavailable. Input sizes are real: sweeps make
`truncate -s <N>K` sparse files so the genuine `du --apparent-size` path runs.

**HOME is redirected into `$SS_SCRATCH`.** `config.txt` sets
`smartSlurmJobRecordDir=$HOME/.smartSlurm` unconditionally and `ssbatch`/`cleanUp`
re-source it, so the only way to isolate the learning DB is to point `HOME` at the
scratch tree — which is exactly how production resolves the path.

### Test blocks

| Block | What it protects |
|---|---|
| `scaling_linear` | Linear mem/time vs size; real gnuplot fit recovers slope/intercept; prediction tracks the line. |
| `noise_fit` | ~8% noise; fit still recovers the slope; `STDFIT>0` reported. |
| `scaling_exp` | Exponential truth vs SmartSlurm's **linear** fit — the model mismatch under-provisions at large sizes. |
| `contaminate_exit0` | An exit-0-empty / mis-sized COMPLETED row drags the fitted slope down. |
| `contaminate_exclude` | OOM/failed rows are **not** fit directly (only `COMPLETED` is). |
| `contaminate_prune` | An OOM storm trips fail-pruning; asserts the evolved contract — the failing program's **failed + failure-disproven** COMPLETED rows (bogus under-reporters) go while its **validated COMPLETED history is kept**, and unrelated programs are fully spared. Forces the triggering OOM (`MOCK_FAIL_MODE=oom`) since improved estimation would otherwise prevent it. Slower (~45 s: real requeue). |
| `negslope_robust` | Numeric edge cases (negative/sub-1 line results, empty STDFIT) always yield a valid floored integer — never the old `bc "-"` garbage (issue A3). |
| `constant_model` | Constant mem (size-independent) + scaling time: mem uses the constant (percentile) model, not a negative-slope collapse to the floor; time still scales (A1/A2). |
| `failure_floor` | Estimates never drop to/below a proven-insufficient allocation; a genuinely smaller input is unaffected (B2). |
| `manual_rerun_floor` | End-to-end: after a real emergent-OOM run, a fresh (manual) rerun estimate exceeds the proven-insufficient allocation (B2/B3). |
| `firstpass_rate` | Measured first-pass success tracks `targetFirstPassPct` and responds monotonically to the knob (B1/B4). |

### Resource-estimation quality (issues A and B)

The estimator (`estimateResource.sh` + `calculateMemTime.sh`) was hardened so the
initial estimate aims for ~90–95% first-pass completion, auto-rerun handles the
rest, and **no submission is ever at or near a proven-insufficient allocation**:

- **A — graceful slope handling.** Robust numerics (no `bc "-"`/empty-STDFIT
  breakage); per-dimension model selection — a positive slope uses the fitted line,
  a non-positive/flat slope (e.g. constant memory while only time scales, or
  scheduler-timing noise) falls back to a size-independent constant (percentile)
  model instead of collapsing to the floor.
- **B1/B4 — service level.** The line margin is the `targetFirstPassPct`-th
  percentile of *relative* residuals (scale-invariant; no normality assumption),
  config knob `targetFirstPassPct` (default 93).
- **B2/B3 — failure floor.** `estimateResource.sh` raises every estimate above any
  proven-insufficient OOM/OOT allocation for a comparable-or-smaller input, escalated
  with the same `nextMemAfterFail`/`nextTimeAfterFail` helpers (config.txt) that
  `cleanUp.sh`'s auto-requeue uses — so initial submit, auto-rerun, and a fresh
  manual rerun stay consistent and never regress to a known-too-small value.

### When do bad/error data contaminate jobRecords? (findings)

1. **Model mismatch.** `estimateResource.sh` fits a straight line. For
   super-linear (e.g. exponential) growth it under-predicts at large sizes →
   systematic OOM once inputs exceed the fitted range (`scaling_exp`).
2. **exit-0-empty / wrong inputSize.** A job that "succeeds" having produced
   nothing records a near-zero mem at its input size; a mis-measured `du` records
   the right mem at the wrong size. Either enters the fit as a `COMPLETED` row and
   skews the slope (`contaminate_exit0`).
3. **Fail-pruning over-deletes (FIXED).** On >`numberOfFailedJobsLastHourToRecalculateStats`
   recent fails, `cleanUp.sh` kept rows matching `$12 != prog && $13 != ref`. To
   drop one program's records it must keep `$12!=prog || $13!=ref`; the `&&`
   also deleted every *other* program sharing the ref — and nearly all jobs use
   `ref=none`, so one program's OOM storm could wipe the whole DB. Fixed in three
   steps: `||` stops the cross-program wipe; a follow-up keeps COMPLETED history
   by pruning only the offending program's failed rows; and a final step also
   drops COMPLETED rows the failures *disprove* — a record `(size s, mem m)` is
   removed when an OOM at size `fs<=s` had a proven-insufficient allocation
   `fm>=m` (same rule on time for OOT). So bogus under-reporters are purged while
   validated history and unrelated programs are kept, and the curve is recomputed
   higher instead of wiped. `contaminate_prune` asserts all of this.
4. **gnuplot-6 scraping drift (FIXED).** `estimateResource.sh` recorded `RSquare=0`
   always under gnuplot 6 (the `Correlation: r = …` line format changed, so
   `cut -d' ' -f7` missed it). Fixed to parse the last field; `noise_fit` asserts
   `RSquare>0.5`. (`STDFIT` is still empty for perfectly collinear noiseless
   history — a zero-residual degenerate case, not seen with real noisy data;
   prediction math is unaffected.)

Both fixes are part of the estimator/jobRecords hardening patch set submitted
upstream to `ld32/smartSlurm` as a pull request; the two blocks above assert the
corrected behavior, so they pass once that PR is in `master` and are red on
unpatched `master` — by design, flagging the gap.

## Adding a fixture

1. `mkdir test/fixtures/my_case` with a `pipeline.sh` (annotated `#@` script) and
   any step tools it calls by bare name (they'll be on `$PATH` in the workspace).
2. Add a block to `run_tests.sh`: `ss_new_workspace`, copy the fixture, run,
   and assert with the helpers (`assert_file`, `assert_nonempty`, `assert_empty`,
   `assert_contains`, `assert_valid_bash`, `assert_invalid_bash`, `assert_rc`).

## Failure injection

Force pathologies without special tools, via env read by the mocks:

```bash
MOCK_FAIL_JOB=<jobName-substring>   # only affect matching jobs (by -J name)
MOCK_FAIL_MODE=nonzero              # payload fails non-zero → FAILED
MOCK_FAIL_MODE=oom                  # job killed → sacct OUT_OF_MEMORY, no .success
MOCK_FAIL_MODE=timeout             # job killed → sacct TIMEOUT, no .success
```

Injected jobs are killed **before** touching `.success` (the injection is applied
in mock `srun`), so `cleanUp.sh` falls through to the `sacct` state exactly as on
a real cluster — an injected OOM job ends with no `.success` and `jobStatus: OOM`.
The `oom_injection` fixture in `run_tests.sh` locks this in. (The exit-0-empty
case is modeled by a fixture's own tool, not injection.)

## Fidelity tiers

- **Tier 1 (this harness):** shims, fast inner loop, no Slurm. Great for rAP
  conversion, the log/data model, diagnosis, and block-syntax work.
- **Tier 2 (recommended for CI before merge):** stock **slurm-in-docker at your
  production Slurm version (25.11.7)** — real scheduler, real `sacct`. Catches
  format/behavior drift the shims can't. *Not* the UBCCR `slurm_simulator` fork,
  which simulates scheduling in virtual time and never executes payloads.
- **Tier 3:** real-cluster canaries for partitions, email, and scale.

> **`sacct`/`squeue` format is pinned to production Slurm 25.11.7.** The 25.11
> release changed nothing in `sacct`'s tabular output, and the columns are fixed
> by `cleanUp.sh`'s explicit `--format=` string, so the mock matches the cluster.
> If you move to a Slurm version that changes `sacct` output, capture
> `sacct --version` and a sample row and update `mockSlurm/sacct`. See
> `mockSlurm/README.md`.
