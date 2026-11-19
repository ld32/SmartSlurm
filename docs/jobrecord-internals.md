[&larr; SmartSlurm README](../README.md) &middot; [Docs index](README.md)

# jobRecord.txt internals (for developers)

This is the reference for code that reads or writes [`jobRecord.txt`](core-concepts.md#jobrecords). Everything below lives in one file, `bin/jobRecordSchema.sh`; source it rather than re-implementing any of it.

## Rules for readers and writers

- **Read columns by name, never by position.** Use the awk prelude `JOB_RECORDS_AWK_COLMAP` (`NR==1{for(i=1;i<=NF;i++)col[$i]=i; next}`, then `$(col["memUsed"])`), the bash helper `jobRecordsCol <name> <header>`, or `csv.DictReader` in Python (as `reviewJobRecords.py` does). This is what lets a column such as `notes` be added without breaking anything.
- **Read through `jobRecordsEmitV2 <file>`.** It returns the store in the current format whatever is on disk (old 21-column rows, 17-column rows, mixed files), so a reader never has to care which version a user still has.
- **Write only through `jobRecordsRow` / `jobRecordsAppendRow`** (one row, fields given as `name=value`) and **rewrite only through `jobRecordsRewrite`** (prunes). These place fields in the right column order, replace commas, spaces, tabs and newlines in values with `_` so every row has exactly the header's number of fields, take a lock (a `mkdir` lock, safe on NFS; a stale lock is reclaimed after `JOB_RECORDS_LOCK_STALE` seconds), and convert an old-format store before writing to it. `cleanUp.sh` is the only place that writes rows.

## Columns: why each exists

Eleven columns are **required** (`JOB_RECORDS_REQUIRED_V2`); estimation, pruning and failure analysis do not work without them:

| Column | Used for |
|---|---|
| `pipeline`, `step`, `ref` | the fit key: one memory/time fit per (pipeline, step, reference group) |
| `jobID` | record identity, so a review tool can name a row to delete |
| `inputSize` | x of the regression; size filter in failure analysis |
| `memUsed`, `timeUsed` | y of the memory and time regressions |
| `memAllocated`, `timeAllocated` | out-of-memory / out-of-time analysis: a proven-insufficient floor |
| `status` | keep only COMPLETED history; detect OOM/OOT |
| `epoch` | time windows (recent failures, 7-day prune) |

The rest are cheap to record and kept for people and future analysis; no code may depend on their presence or position: `sample` (from the flag, so a review can name the sample), `cores`, `user` (who ran it, in shared stores), `node` (to spot a slow node), `isoDate` (readable time), `flag` (the full job id) and `notes` (how the row was produced, e.g. what AI detection filled in).

## Version detection

The **header is the version**; there is no marker line. `jobRecordsSchemaVersion <file>` scans every row, not just the first, because stores can be mixed (an old node appending old-format rows into a converted store):

| On disk | Detected as | Action on next write |
|---|---|---|
| header `pipeline,step,…,notes` first, all rows 18 fields | current | none |
| any 20/21-field row with a numeric first field (old format), any 17-field row or header (before `notes`), or no header first | needs conversion | convert in place (backup first) |
| nothing recognisable | unknown | left untouched, with a warning |

Conversion is checked before it is committed: the new file must be in the current format and keep every recognisable row, or the store is left untouched and an error is printed. `jobRecordsRepair.sh --check` prints the same report for a human.

## Mapping from the old 21-column format

Old order: `jobID inputSize memDefault timeDefault memAllocated timeAllocated memUsed timeUsed status user saccMem program ref flag cores extraMem extraTime ratioM ratioT date epoch`.

- **Carried over:** jobID, inputSize, memUsed, memAllocated, timeUsed, timeAllocated, status, `program` → `step`, ref, cores, user, flag, epoch (taken from the last field only if it looks like a Unix time, because the old free-text date could contain a comma).
- **Derived:** `sample` from the flag (pipeline flags only, `step.deps.name.sample`); `isoDate` from the old date; `pipeline` = `legacy` (old rows never recorded it); `node` and `notes` empty.
- **Dropped** (nothing read them): memDefault, timeDefault, saccMem, extraMem, extraTime, ratioM, ratioT.

## Reference pooling

`JOB_RECORDS_REF_CLASS` maps references of similar size to one group so their history is fitted together. The default pools mouse and human (`mm9 mm10 mm39 hg19 hg38 GRCm38 GRCm39 GRCh37 GRCh38` → `mammal`). Set it to `""` to fit every reference separately. `none` (tools without a reference) is never pooled.

## Fitting performance

`estimateResource.sh` builds the memory and time extracts for a fit in **one** pass over the store (`jobRecordsExtractBoth`), and only when the cached fit is stale, because the store is shared and keeps growing. gnuplot does the fit; its results (`Finala`, `Finalb`, `STDFIT`, …) are printed by gnuplot itself rather than scraped from `fit.log`, which avoids a race between concurrent jobs.

---
