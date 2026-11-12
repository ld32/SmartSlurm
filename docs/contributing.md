[&larr; SmartSlurm README](../README.md) &middot; [Docs index](README.md)

# Contributing

- **Test off-cluster first.** Develop and validate changes with the
  [slurmless](testing-off-cluster.md#testing-off-cluster-slurmless) harness
  (`./bootstrap.sh`, then `./test/run_tests.sh <block>`) before running them on a
  real cluster.
- **Smoke-test with the examples.** The pipelines under `example-scripts/` double
  as templates and quick cluster smoke tests; each directory has its own README
  with the exact commands and expected output.
- **Migrating an old store?** Convert a legacy jobRecord store once with
  [`jrRepair.sh`](tools.md#under-the-hood-tools) (keep the timestamped backup it writes).

---

[&larr; SmartSlurm README](../README.md) &middot; [Docs index](README.md)
