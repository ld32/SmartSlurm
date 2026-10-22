# resource-estimation testbed

A tidied, self-contained replacement for the old `useMemTime*` / `findNumber.sh` resource
scripts. It exercises SmartSlurm's jobRecord machinery: recording mem/time per input size,
fitting a curve once there are ≥3 records, estimating resources for new inputs, and
incrementing resources on an OOM/OOT **auto-rerun** (or **manual rerun** if that's
exhausted).

## Files

- **`sizedLoad.sh`** — the load "tool". Given an input file, it holds MEMORY and burns
  WALL-TIME in proportion to the file's row count, deterministically. Cleaner than the old
  approach: memory is an `awk` array of distinct ~1 KB entries (linear, predictable RSS)
  instead of bash string concatenation / `eval array$i=…`.
- **`jobrecords_pipeline.sh`** — a one-step rAP pipeline that runs `sizedLoad.sh` once per
  sample. The `load` step is the thing whose resources get recorded and fitted.
- **`jobrecords_manifest.txt`** — samples at graded sizes (`s1 200 … s4 1600`), so a single
  run produces four records across a size range — enough to fit.

## `sizedLoad.sh` knobs (env vars)

| Var | Default | Meaning |
|---|---|---|
| `MEM_KB_PER_ROW` | 20 | memory held per input row (KB) |
| `SEC_PER_ROW` | 0.02 | seconds slept per input row |
| `NOISE_PCT` | 5 | ± jitter on both totals, so records aren't perfectly collinear |
| `FAIL_PROB` | 0 | chance [0..1] of a deliberate overshoot |
| `FAIL_MODE` | both | `oom` \| `oot` \| `both` — what to overshoot |
| `OVERSHOOT_X` | 4 | multiplier applied on overshoot |

The pipeline forwards these from your submit shell, baking them into each job.

Measured scaling (single host): 200/400/800 rows → ≈10/14/23 MB peak RSS and ≈4/8/16 s —
linear in rows (atop a fixed ~6 MB `awk` baseline). Tune the knobs / manifest sizes to land
where you want.

## Test protocol

Copy all three files into a fresh run directory, with `SmartSlurm/bin` on `PATH`.

**1. Clean run — seed the fit (all jobs succeed):**
```bash
runAsPipeline --script "jobrecords_pipeline.sh jobrecords_manifest.txt" \
              --sbatch-options "sbatch -p short -c 1 --mem 64M -t 2:00" --tmp noTmp
```
Four `load` jobs run at four sizes → four jobRecords. With ≥3 records, the next run fits a
mem/time curve instead of using defaults (watch for the "fit" vs "less than 3 records"
message from `estimateResource.sh` in `.smartSlurm.log`).

**2. Rerun test — force an overshoot past the request:**
```bash
FAIL_PROB=1 FAIL_MODE=oom OVERSHOOT_X=4 \
  runAsPipeline --script "jobrecords_pipeline.sh jobrecords_manifest.txt" \
                --sbatch-options "sbatch -p short -c 1 --mem 64M -t 2:00" --tmp noTmp
```
The large sample's `load` job exceeds `--mem` → OOM → SmartSlurm requeues it with more
memory (auto-rerun); if the increment is exhausted, the manual rerun menu takes over. Use
`FAIL_MODE=oot` (with a small `-t`, e.g. `-t 1:00`) to test the time path instead.

## Caveats

- **Real OOM/OOT fires only on real Slurm.** In the mock harness, failures are *simulated*
  via `MOCK_FAIL_MODE`; the actual cgroup OOM / wall-clock TIMEOUT that drives requeue only
  happens on the cluster.
- Memory has a fixed ~6 MB `awk` baseline; size your `--mem` and `MEM_KB_PER_ROW` so the
  *scaling* part dominates and a clean run sits comfortably under the request.
- Keep sizes modest (the defaults finish in seconds) so a full accumulate-and-fit cycle is
  quick.
