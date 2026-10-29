# jobRecords schema revision — design document

Status: IMPLEMENTED on `feat/jobrecords-v2` (off ld32 master 71e9d9c). This is the as-built design record for the jobRecords v1→v2 revision. (fitQC, a separate tool, was paused on this and resumes after.)

### Progress
- [x] **P1 SSOT** — `bin/jobRecordSchema.sh` (commit 62ae051). Additive, nothing sources it yet.
- [x] **P2 repair tool** — `bin/jrRepair.sh` (commit 62ae051). Tested on synthetic v1 stores (mapping, sample-derivation per-sample/experiment/nested, idempotent, dry-run writes nothing, malformed-row diversion, backup). **TODO: add the matching harness block in the standalone slurmless repo.**
- [x] **P3 writer** — commit 7cc4d93. cleanUp record=→v2; ssbatch/ssbatchSimple parse `--pipeline=`, write v2 header via SSOT, quote excludeFailedNodes in the trap; runAsPipeline `--pipeline` (default=script basename) injected as `--pipeline=` into every ssbatch call. Pipeline flag is **`--pipeline=NAME`** (both `-N` and `-L` were taken by Slurm --nodes/--licenses). Verified: clean_2step writes a correct v2 row.
- [~] **P4 readers** — IN PROGRESS.
  - [x] **P4a** (commit d0aa066): compat shim `jrEmitV2` + robust `jrSchemaVersion` (handles headerless v1); jrRepair refactored onto the shim; estimateResource extracts + fail-analysis → shim+name-based; cleanUp 7-day prune version-branched (v1 unchanged, v2 name-based, header-preserving) + records read via shim. VERIFIED: history persists; estimation works on v1+v2. Compat strategy: read-only readers use jrEmitV2 (no file change); write-backs version-branch (v1 byte-identical, no surprise migration).
  - [x] **P4b** (commit 7a94aff): cleanUp fail-prune block → v2. `jrFailPrune()` shared helper (version-branched); `recentFails` via shim. Verified in isolation.
  - [x] **P4c** (commit af8d3ad): fit per (pipeline, step, refGroup). `jrExtract()` with FALLBACK LADDER (pipeline-scoped → pool step+refGroup when <3 pts → default); statKey names stats; pipeline threaded ssbatch→estimateResource→calculateMemTime; ref size-class (jrRefClass/jrRefMatchRe, default pools mm*/hg*/GRC*). Verified on scaling_linear: pipeline-keyed stats, exact fit (0.5x+200,R²=1), fallback pooled correctly.
  - [x] **P4d** (commit 3382fcb): jrExtractBoth() builds mem+time in ONE jrEmitV2 pass; gnuplot emits STDFIT from stdout (print FIT_STDFIT) with fit logfile silenced, replacing the race-prone `cat fit.log` scrape. Re-validated: estimator 16/17 (only extrapolate_lowdata).
  - NOTE: harness assertions that hard-code v1 column positions (e.g. scaling_linear:457/460 read $9/$12) FAIL on v2 output — these are slurmless-side updates (with P5/P6), not code bugs.
- [~] **P6 harness updates** — DONE for run_tests.sh/harness.sh (slurmless commit c97104b, pushed). ss_seed_record→v2 rows; ss_fit→v2-keyed stats (unknown.PROG.REF); JR_REF_CLASS="" for determinism; all block jobRecord reads re-indexed to v2 cols; stat refs re-keyed. **38/40 blocks green.** The 2 reds are unmerged-feature gaps, NOT v2: extrapolate_lowdata (thin-data fork) and the .out.<jobid> marker check (retire-unused PR not in this base). Both go green once those branches reach the base. Remaining P6: docs/release note ("run jrRepair.sh once"); commit this design doc.
- [x] **P5** (commit ad3488f): reviewJobRecords.py → csv.DictReader + (pipeline,step,ref) grouping; v1 refused with jrRepair pointer. Verified it separates rnaseq/bwa/hg38 from wgs/bwa/hg38.


## DECISIONS (signed off)
1. Name-based access + v1 compat shim — **YES**.
2. Optional columns — **add `sample` and `node`; drop `memReqDefault`/`timeReqDefault`**.
3. Retire set — **drop all**: `saccMem`, `ratioM`, `ratioT`, old `date`, `extraMem`, `extraTime` (and the req-defaults from #2).
4. Versioning — **by header signature, NOT a `#comment` line** (we read the header anyway for name access; v1 = `1jobID,2inputSize,…`, v2 = `pipeline,step,ref,…`). No reader needs `/^#/` skipping.
5. Delivery — **one coordinated PR (P1–P6)**.

Final v2 schema (17 cols, from 21):
`1 pipeline · 2 step · 3 ref · 4 sample · 5 jobID · 6 inputSize · 7 memUsed · 8 memAllocated · 9 timeUsed · 10 timeAllocated · 11 cores · 12 status · 13 user · 14 node · 15 epoch · 16 isoDate · 17 flag`
Required (11): pipeline, step, ref, jobID, inputSize, memUsed, timeUsed, memAllocated, timeAllocated, status, epoch.

### Guiding principle (overrides style ties)
**Robust and clearly coded first; snappy second.** No clever one-liner that obscures intent. Every reader/writer change is covered by a harness block. Perf wins only if they don't cost clarity or robustness.

### Logistics
- Branch: **`feat/jobrecords-v2`** off **ntc-agnt's fork of ld32 `master`** (one coordinated PR → ld32).
- Harness: push to the **standalone slurmless repo** only. Merge into SmartSlurm's vendored copy is **on hold** (no `git subtree` sync until dev settles).
- Rollback: `jrRepair.sh` keeps timestamped backups of the pre-conversion file(s) (`jobRecord.back.<ts>.txt`, and the `.oom.oot.over` variant) so a user can restore v1 if needed. Never delete the backup.

## 0. Objective & constraints

Revise the SmartSlurm jobRecord schema so it is honest, minimal, name-addressed, and
ready for per-(pipeline, step, ref) fitting — **without breaking existing deployments**.

Hard constraints:
- No silent breakage of existing `~/.smartSlurm/jobRecord.txt` (and shared group stores).
- Conversion of existing files is a **manual, one-time, user-run** repair utility (not auto-on-access).
- The offline python review tool (`reviewJobRecords.py`) is revised **in tandem**.
- Don't modify tool behavior beyond the schema/consumers; CLAUDE.md "change a tool only when the task is to" — this IS that task.

## 1. Verified producer/consumer survey (exhaustive)

**Writers**
- Rows: `cleanUp.sh` ONLY — `record=` at 202/280/569; appended at 235/282 (`jobRecord.txt`), 570 (`jobRecord.oom.oot.over.txt`).
- Header: `ssbatch:90` + `ssbatchSimple:85` — two identical copies of the stale header.

**Readers of jobRecord.txt — all POSITIONAL (awk `$N`)**
- `estimateResource.sh`: extracts (44 `$7,$8`; 115 `$2,$7`; 117 `$2,$8`); fail analysis (193 `$9,$21`; 194 `$5,$2,$9`; 195 `$6,$2,$9`); keys `$12,$13`.
- `cleanUp.sh`: prune (227 `$21,$9`; 230 `$2,$7`; 252 `$21`), fail-prune (295/582 `$5,$6,$7,$8,$2,$9`); keys `$12,$13`.
- `reviewJobRecords.py`: `row[0]`=jobID, `row[1]`=inputSize, `row[6]`=memUsed, `row[11]`=program (0-indexed, positional).

**Readers of 2-col extracts** `.mem`/`.time`/`.memTime.noInput` (POSITIONAL, `cut -f1/-f2`): `calculateMemTime.sh`, `estimateResource.sh`.

**Reader of `.stat`** (NAME-based, sourced): `calculateMemTime.sh` → `Finala Finalb STDFIT RSquare SCount Maximum`. ← the pattern to emulate.

**TRUE consumed jobRecord fields** (awk-internal, shell-arg pollution removed):
`2 inputSize, 5 memAllocated, 6 timeAllocated, 7 memUsed, 8 timeUsed, 9 status, 12 program/step, 13 ref, 21 epoch` (+ `1 jobID` by review tool).
**Not consumed by anything**: 3 memDefault, 4 timeDefault, 10 user, 11 saccMem, 14 flag, 15 core, 16 extraMem, 17 extraTime, 18 ratioM, 19 ratioT, 20 humanDate.

## 1b. Pipeline propagation (how the new `pipeline` value reaches the row)

`cleanUp.sh` writes the row but doesn't know the pipeline. Path to wire (P3):
`runAsPipeline` knows it (`--script`, or a new one-time `--pipeline NAME`, default = script basename with a generic-name warning) → stamps it on each generated `ssbatch` call (new flag, e.g. `-N <pipeline>`) → `ssbatch`/`ssbatchSimple` carry it into the per-job env/cmd dir → `cleanUp.sh` reads it into `$record`. Old rows (pre-flag) and direct `ssbatch` use with no `-N` → `pipeline=unknown`.

## 2. Q1 — Minimum field set that satisfies all consumers (incl. fitQC + pipeline)

Eleven fields:

| field | why required |
|-------|--------------|
| `pipeline`   | NEW — fit key component (fitQC) |
| `step`       | fit key (was `program`/`-P`) |
| `ref`        | fit key (may be `none`) |
| `jobID`      | record identity — fitQC/review select a row to remove |
| `inputSize`  | regression X; fail-analysis size filter |
| `memUsed`    | regression Y (mem) |
| `timeUsed`   | regression Y (time) |
| `memAllocated`  | OOM fail-analysis (proven-insufficient floor) |
| `timeAllocated` | OOT fail-analysis |
| `status`     | COMPLETED filter; OOM/OOT detection |
| `epoch`      | time-window filters (recent-fails, prune) + date math |

Nothing else is required for correct estimation, pruning, or fitQC.

## 3. Q2 — Additional fields: useful + cheap (readability / future dev)

All are available at write time in `cleanUp.sh`/`ssbatch` with no extra queries → "cheap".

| field | value | keep/add rationale |
|-------|-------|--------------------|
| `user`      | `$USER` | per-user disambiguation in shared/group stores; readability ("who ran it"). Already captured — keep, rename from `useID`. |
| `flag`      | rAP `step.deps.name.sample` | sample-level identity — lets fitQC say *which sample* is an outlier; future per-sample analysis. Keep. |
| `cores`     | `-c` | multi-population diagnostic; future per-core normalization. Keep. |
| `isoDate`   | ISO-8601 of epoch | human readability (the one the user asked for); replaces the free-form `date`. |
| `node`      | `$SLURMD_NODENAME`/host | ADD — future node-anomaly diagnosis (a slow node → apparent outlier). Cheap. |
| `memReqDefault` / `timeReqDefault` | pre-estimate request | OPTIONAL — readability of the request→allocate→use chain; future "did estimation beat the default" analysis. |

Retire (unused AND superseded/redundant): `saccMem` (second mem measure, unused), `ratioM`/`ratioT` (superseded by efficiency.tsv + checkRun), old free-form `date` (replaced by `isoDate`), `extraMem`/`extraTime` (unused), `timeDefault` (unused; keep only if we adopt the optional req-default pair).

## 4. Proposed v2 schema

Column order (identity+keys → request/alloc/use chain → outcome → provenance). Required = ✱.

```
1 pipeline ✱  2 step ✱   3 ref ✱     4 sample    5 jobID ✱
6 inputSize ✱ 7 memUsed ✱ 8 memAllocated ✱ 9 timeUsed ✱ 10 timeAllocated ✱
11 cores     12 status ✱ 13 user     14 node     15 epoch ✱
16 isoDate   17 flag
```
- `sample` (4): parsed from `flag` at write time — explicit column for readability/fitQC.
- `node` (14): `$SLURMD_NODENAME`/host at write time — future node-anomaly diagnosis.
- Non-required columns (sample, cores, user, node, isoDate, flag) are kept for readability/future; no consumer depends on them, so they can be reordered/added later safely once access is name-based.
- Dropped vs v1: memDefault, timeDefault, saccMem, extraMem, extraTime, ratioM, ratioT, old free-form date.

**Versioning = header signature** (decision #4): no comment/marker line. Readers already parse the header for name access; the header's names ARE the version — v1 = `1jobID,2inputSize,…`, v2 = `pipeline,step,ref,…`. `jrRepair.sh` and the compat shim branch on that. No reader needs `/^#/` skipping. Future v2→v3 detected the same way (new/renamed columns in the header).

## 5. Access model — name-based (header-driven)

Replace positional `$N` with name resolution from the header. Shared awk prelude (one definition, sourced or pasted):
```
NR==1 { for(i=1;i<=NF;i++) c[$i]=i; next }      # c["memUsed"] → column index
```
- `.stat` stays name-based (already is).
- Bash callers that need a column: a `jrcol <name>` helper or the prelude.
- `reviewJobRecords.py`: `csv.DictReader` (names, not `row[N]`).

**v1 compatibility shim**: when a reader sees a v1 file (no v2 marker / old header), it applies the known v1→v2 column map internally so estimation never degrades before the user runs the repair. The shim is removable in a later release once migration is widespread.

## 6. Repair / update utility (manual, one-time) — `jrRepair.sh`

```
jrRepair.sh [jobRecordDir]        # default $smartSlurmJobRecordDir
jrRepair.sh --dry-run             # report what it would change, write nothing
```
Behavior:
1. Detect format by header signature: header already v2 (names `pipeline,step,…`) → report "already current", exit 0 (idempotent). Old `1jobID,…` header → convert.
2. Back up: `jobRecord.back.<ts>.txt` (existing convention).
3. Map v1 (21 positional fields) → v2 named schema:
   - carry: jobID, inputSize, memUsed, memAllocated, timeUsed, timeAllocated, status, step(=program), ref, cores, user(=useID), flag, epoch.
   - synthesize: `isoDate` from epoch; `sample` from flag; `pipeline` = `legacy` (unrecoverable from old rows — stated honestly); `node` = empty.
   - drop: saccMem, ratioM, ratioT, old date, extraMem, extraTime, memDefault, timeDefault.
4. Write new file with the v2 header (signature versioning; no marker line).
5. Also offer to migrate `jobRecord.oom.oot.over.txt` and leave `.back.*` archives untouched.

Manual by design: no surprise rewrites mid-run; the user runs it once at upgrade. Until they do, the compat shim (§5) keeps things working.

## 7. reviewJobRecords.py revision (in tandem)

- Read via `csv.DictReader` (name-based) — robust to column order; drop all `row[N]`.
- Group by the new key **pipeline + step + ref** (ref via size-class map, §8); `ref==none` drops ref.
- Keep its three functions (assess ASCII plot / identify / delete-range), now operating on the new schema and grouping.
- Detect v2 by header signature; if handed a v1 file, print "run jrRepair.sh first" (python tool is offline, so a hard prompt is fine there).
- Longer term it converges with fitQC, but for this work it must simply not break and must use names.

## 8. ref size-class merge (config option)

- A config map groups similarly-sized genomes: e.g. `REF_SIZE_CLASS="mm9 mm10 mm39 hg19 hg38 => mammal"` (mouse+human pooled per your call).
- Applied at grouping/extract time: the key uses `classOf(ref)` when enabled, raw `ref` otherwise. Default off (hardcoded table in config for now; full UI later).
- `ref==none` → ref dropped from the key entirely (ref-less tools).

## 8c. Performance & robustness review (read/write + gnuplot)

Priorities in order: **robust → clearly coded → snappy**. Keep gnuplot as the fitter (do NOT reimplement the fit in awk).

**Read hotspot (the real cost).** `estimateResource.sh` rebuilds the `.mem` and `.time` extracts by scanning the whole — shared, unbounded-growth — `jobRecord.txt` with `grep COMPLETED | awk`, **twice** (once per dimension), potentially on every submit. In a group store that's O(N) per job with N climbing forever.
- In-scope fix (clear + safe): derive **both** extracts in **one** pass (single awk emitting mem+time), keyed by `(pipeline,step,ref)`; keep the existing mtime/row-count cache. Fewer full scans, same logic, still legible.
- Deferred (flagged, not now): incremental per-group extract files appended by `cleanUp` on completion (no full rescans). More moving parts → only if the single-pass isn't enough. Robustness > speed.

**Write path.** `cleanUp` appends one row per job — fine; leave it. v2 is 17 cols vs 21, so rows are marginally smaller/faster to scan (incidental win).

**gnuplot data passing (robustness, keep the tool).** Today it passes data via the clean extract file (good) but then **scrapes stdout** with `grep/sed/awk` for `Finala/Finalb` and re-reads `fit.log` for `STDFIT` — fragile, and it already has to guard empty/collinear cases. Clean-ups: stabilize the output parse, handle the degenerate cases explicitly (noted: STDFIT empty for perfectly collinear history), and avoid re-reading the same extract more than once per dimension. One gnuplot invocation per dimension is acceptable; don't multiply calls.

**Name-based access cost.** One header parse per read — negligible; the robustness win dominates.

Net: schema shrink is perf-neutral-to-slightly-positive; the meaningful speed win is single-pass extract derivation; the meaningful robustness win is name-based access + a clean gnuplot output parse. All covered by harness blocks.

## 9. Phased work plan (each step shippable, harness-verified)

- **P1 — Schema SSOT**: one authoritative column definition (names+order+marker) in a shared location; dedupe the two header-writers to read it. No behavior change yet.
- **P2 — Name-access helper + `jrRepair.sh`, harness-first**: build the prelude/helper, the v1→v2 map, and the repair utility; prove on slurmless fixtures (v1 file, v2 file, ref-less, multi-pipeline, idempotency, dry-run) BEFORE touching production readers.
- **P3 — Writer**: `cleanUp.sh` `record=` emits v2; single shared header with marker.
- **P4 — Readers → name-based, one at a time** (each harness-verified): estimateResource extracts + fail-analysis, cleanUp prune, calculateMemTime extract key. Add the v1 compat shim. **Fold in the §8c perf/gnuplot clean-ups here**: single-pass `(pipeline,step,ref)` extract derivation + robust gnuplot output parse.
- **P5 — reviewJobRecords.py** to DictReader + new grouping.
- **P6 — ship** with release note: "run `jrRepair.sh` once to convert existing history."

fitQC resumes after P4 (its grain = v2 key), reusing the same name-access helper.

## 10. Test plan (slurmless harness)

Fixtures: a v1 jobRecord.txt (old 21-field), a v2 file, a ref-less group, a two-pipeline collision group.
Assert: repair maps v1→v2 correctly, is idempotent, backs up, `--dry-run` writes nothing; readers produce identical fits on v2 as the old code did on v1 (regression parity); compat shim reads v1 unchanged; grouping splits by pipeline+step+ref; ref-class pools mm/hg.

## 11. Open decisions (need sign-off)
1. Adopt name-based access + v1 compat shim? (vs positional forever) — recommended yes.
2. Optional columns: include `sample`, `node`, `memReqDefault`/`timeReqDefault`? (readability vs width)
3. Retire set: confirm `saccMem`, `ratioM`, `ratioT`, old `date`, `extraMem`, `extraTime` go.
4. Schema-version marker as a `#comment` first line — acceptable to all parsers? (grep-based readers must skip it.)
5. One coordinated PR (P1–P6) vs land P1 (honest SSOT header) first, then the rest.
