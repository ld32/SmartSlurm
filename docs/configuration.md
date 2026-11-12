[&larr; SmartSlurm README](../README.md) &middot; [Docs index](README.md)

# Configuration

## config.txt

Ships at `SmartSlurm/config/config.txt`. It defines the partition time limits, the `adjustPartition` function, and defaults such as where records and logs live.

## Key settings

```bash
export smartSlurmJobRecordDir=$HOME/.smartSlurm  # where jobRecord.txt + stats live
export smartSlurmLogDir=smartSlurmLog            # per-run log dir (relative to your cwd)

export firstBatchCount=5   # runAsPipeline: how many independent jobs run before the rest are held
export defaultMem=4096     # M   — used until estimation is available
export defaultTime=120     # min — used until estimation is available
export defaultExtraMem=500   # M   — safety margin added to estimates
export defaultExtraTime=10  # min — safety margin added to estimates

export partition1TimeLimit=720    # minutes: run-time  >0    and ≤12h  (720 min)
export partition2TimeLimit=7200   # minutes: run-time >12h   and ≤5 days (7200 min)
export partition3TimeLimit=43200  # minutes: run-time  >5d   and ≤30 days (43200 min)

adjustPartition() { ...; }
```

Other settings referenced elsewhere in this document:

- `effLowPct` (default **25**) — efficiency at or below this percent turns a [resource-efficiency](monitoring.md#resource-efficiency) bar red.
- `JR_REF_CLASS` — controls [ref pooling](core-concepts.md#jobrecords) (default pools `mm*`/`hg*`/`GRC*`; set to `""` to disable).

> [!IMPORTANT]
> **`firstBatchCount` (5) and the "3 records" rule are different things.**
> - **3 records** is the minimum ssbatch needs before it can *estimate* resources for a program.
> - **`firstBatchCount=5`** is used only by `runAsPipeline`: it lets the first 5 independent jobs of a step run immediately, and submits the rest **held** (`-H`) until early jobs finish and produce records, then releases them with estimated resources.

## Personal overrides

> [!TIP]
> Copy `config.txt` to `~/.smartSlurm/config/config.txt` to override settings for just yourself. Your personal copy wins over the shared one — handy on a shared cluster. (SmartSlurm refuses to run if your personal copy is *older* than the shared one, to stop you using a stale config; see [Upgrading](upgrading.md#upgrading) to refresh it.)

---


---

[&larr; SmartSlurm README](../README.md) &middot; [Docs index](README.md)
