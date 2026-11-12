[&larr; SmartSlurm README](../README.md) &middot; [Docs index](README.md)

# Testing off-cluster (slurmless)

You can develop and test SmartSlurm **without a Slurm scheduler** using the
**slurmless** harness — a Slurm-free sandbox that stands in a mock scheduler so
`runAsPipeline`, `ssbatch`, `cleanUp.sh`, `checkRun`, `checkRunGrid`,
`diagnose.sh`, and resource estimation can be exercised on an ordinary machine.
SmartSlurm itself is **not modified**; the mocks actually execute each job's
payload, so the real cleanup trap and `.success`/`.failed` flags are produced
exactly as on a cluster.

slurmless is vendored into this repository at [`slurmless/`](../slurmless/), so
from a SmartSlurm checkout:

```bash
cd SmartSlurm
./slurmless/bootstrap.sh                 # prove the baseline first (clean_2step)
./slurmless/test/run_tests.sh --list     # see every test block and what it checks
./slurmless/test/run_tests.sh <block>    # run one block (the full suite is slow by design)
```

The harness resolves the enclosing SmartSlurm checkout automatically, so no path
setup is needed when it runs from `slurmless/`.

slurmless is also developed as a standalone repo,
[`ntc-agnt/slurmless_smartSlurm_testing`](https://github.com/ntc-agnt/slurmless_smartSlurm_testing#readme).
See [`slurmless/README.md`](../slurmless/README.md) for the full workflow: adding
a test fixture, injecting OOM/timeout/failure, and converting a pipeline to
`runAsPipeline` and validating it off-cluster.

---

[&larr; SmartSlurm README](../README.md) &middot; [Docs index](README.md)
