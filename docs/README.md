# SmartSlurm documentation

Detailed guides for SmartSlurm. For the overview, see the
[main README](../README.md).

### Getting started
- [Installation & setup](installation.md) — clone, the `smartslurm` command, activate/deactivate, optional Graphviz.
- [Quick start](quick-start.md) — a first `ssbatch` job and a first `runAsPipeline` pipeline, end to end.

### Using SmartSlurm
- [Core concepts](core-concepts.md) — how a pipeline runs (convert &rarr; submit), the `#@` job block, scopes, loops, resource estimation, jobRecords, notifications.
- [Tools reference](tools.md) — `runAsPipeline`, `ssbatch`, `reviewJobRecords.py`, `cancelAllJobs`, `sbatchAndTop`, and the under-the-hood helpers.
- [Monitoring & diagnosing runs](monitoring.md) — the `checkRun` status grid, `diagnose.sh`, the resource-efficiency view, and rerunning.
- [Configuration](configuration.md) — `config.txt`, key settings, personal overrides.
- [AI assistance (optional)](ai-assistance.md) — automatic program/input detection for `ssbatch`, setup and settings.

### Going further
- [Using ssbatch with other pipeline managers](pipeline-managers.md) — Snakemake, Nextflow, Cromwell.
- [Testing off-cluster (slurmless)](testing-off-cluster.md) — the Slurm-free sandbox, vendored at `slurmless/`.
- [FAQ & troubleshooting](faq.md) — common `ssbatch` and `runAsPipeline` questions.
- [Upgrading](upgrading.md)
- [Contributing](contributing.md)

### Design notes
- [jobRecord.txt internals](jobrecord-internals.md) — rules for code that reads or writes the store, columns, version detection, v1&rarr;v2 mapping.

---

[&larr; SmartSlurm README](../README.md)
