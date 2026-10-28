# slurmless — Context & Setup

**slurmless** is a **Slurm-free test harness** for SmartSlurm. It runs the
**unmodified** SmartSlurm tools (`runAsPipeline` / `ssbatch` / `cleanUp.sh` /
`checkRun` / `checkRunGrid` / `diagnose.sh`) on a machine with **no Slurm** — a
laptop, CI, or an AI-agent sandbox — while still executing real job payloads, so
pipeline conversion, dependency wiring, failure attribution, the status grid, the
rerun menu, and the resource-estimation learning loop can all be tested for real.

It is **transparent to SmartSlurm**: nothing in SmartSlurm's `bin/` is modified.
The harness interposes mock Slurm commands (`sbatch`, `sacct`, `srun`, …) on
`$PATH` ahead of the real binaries, exactly the way SmartSlurm itself locates
`sbatch`.

---

## Where this lives

slurmless is a **component of SmartSlurm** — ordinary files living at
`SmartSlurm/slurmless/`, developed **in-tree**. You edit `slurmless/…` on a
branch and open a pull request to SmartSlurm, exactly like any other part of the
repo; there is no separate slurmless repository to keep in sync.

```
SmartSlurm/                      # github.com/ld32/smartSlurm
├── bin/                         # the tools under test
├── config/
└── slurmless/                   # <- this directory (a component of the package)
    ├── CONTEXT.md               # this file
    ├── README.md
    └── test/
        ├── run_tests.sh
        ├── checkrun_ui_tests.sh
        ├── lib/harness.sh
        ├── mockSlurm/
        ├── fixtures/
        ├── models/
        └── version.txt
```

- **Project:** `github.com/ld32/smartSlurm`. The `ux-update` branch has been
  merged into `master` and deleted; **`master` is the target.**
- **One copy, versioned with the package.** Because slurmless is not a separate
  repo, a given SmartSlurm commit contains exactly one slurmless — so the harness
  and the tools it tests always share a commit and **cannot drift apart**. Just
  run the harness before pushing a tool change (`./slurmless/bootstrap.sh`, then
  the relevant blocks). `test/version.txt` records the baseline the harness was
  last validated against, for anyone running it against an *external* SmartSlurm
  checkout.

> **Resource-estimation patches:** a set of estimator / jobRecords hardening
> patches (robust numerics, per-dimension model selection, failure floor,
> service-level margin, conservative fail-pruning) has been submitted upstream to
> `ld32/smartSlurm` as a pull request. Several test blocks below assert that
> hardened behavior, so they pass once those patches are in `master` and are
> red until then — by design, flagging the gap. The patches are **not** kept
> here; they belong to SmartSlurm's `bin/`.

---

## Run it

From a SmartSlurm checkout, `SMARTSLURM_ROOT` auto-resolves:

```bash
cd SmartSlurm
./slurmless/test/run_tests.sh                 # every block (SLOW — see note)
./slurmless/test/run_tests.sh clean_2step     # one block by name (substring match)
./slurmless/test/run_tests.sh --list          # list blocks + what each verifies
```

Exit status is non-zero if any assertion fails (CI-friendly). Each block prints
`PASS - …` / `FAIL - …`, and the run ends with a grouped report and
`Blocks: …  Checks: …  RESULT: PASS|FAIL`.

> **Run one block at a time in a time-boxed environment.** The mock scheduler
> sleeps through SmartSlurm's cleanup counters, so the full suite is slow; a
> single block is typically well under ~4 minutes. The gnuplot curve-fit blocks
> (`scaling_*`, `noise_fit`, `contaminate_*`) need `gnuplot` and are slower.

### Running against a different SmartSlurm checkout

To test the harness against a SmartSlurm clone that is *not* its parent directory
(e.g. an older release, or a separate working copy), point it at the checkout
explicitly:

```bash
export SMARTSLURM_ROOT=/path/to/SmartSlurm
/path/to/slurmless/test/run_tests.sh clean_2step
```

`SMARTSLURM_ROOT` resolution: an explicit value wins; otherwise the harness walks
up from `test/` to the first directory containing `bin/runAsPipeline`. That makes
it correct whether slurmless sits at `SmartSlurm/slurmless/test/` (its home) or
`SmartSlurm/test/` (legacy), at any depth.

---

## Interactive use (experimentation / debugging)

```bash
cd SmartSlurm
source slurmless/test/lib/harness.sh
ss_activate                                    # prepares the env (see below)
ss_new_workspace /tmp/play                      # clean scratch workspace, cd's in
cp slurmless/test/fixtures/clean_2step/* .
ss_run_pipeline "./pipeline.sh 777" "sbatch -p short -c 1 --mem 4G -t 50:0"
cat smartSlurmLog/allJobs.txt                   # inspect the result
checkRunGrid --log smartSlurmLog                # try the grid, addressing, why, etc.
```

### Harness API (in `test/lib/harness.sh`)
- `ss_activate` — one-time environment setup (below).
- `ss_new_workspace DIR` — fresh scratch dir, `cd`s in.
- `ss_run_pipeline "SCRIPT [args]" "SBATCH_OPTS" [--dryrun]` — runs
  `runAsPipeline`; full output captured to `$LAST_RUN_LOG`.
- Resource-model API: `ss_use_model`, `ss_run_sweep`, `ss_fit`,
  `ss_seed_record`, `assert_near` (see `test/README.md` → "M1").
- Assertions: `assert_file`, `assert_no_file`, `assert_contains`, `assert_rc`,
  `assert_nonempty`, `assert_empty`, `assert_valid_bash`, and `_pass`/`_fail`.

### What `ss_activate` does (and its side effects)
1. **Installs `gnuplot`** if missing (needed by `estimateResource.sh` for curve
   fitting). Non-fatal: curve-fit blocks skip cleanly when it is unavailable.
2. **Repoints `/bin/sh` → `bash`** if `/bin/sh` is `dash`. SmartSlurm's scripts
   declare `#!/bin/sh` but use bash features; on the **production RHEL cluster
   `/bin/sh` is bash**, so this is faithful emulation, not a hack. It needs root;
   where it cannot (locked-down host), it warns, and the `#!/bin/sh` scripts may
   misbehave under dash — run on a box where `sh` is already bash.
3. Stages the mock commands into a **neutral path** (`$SS_SCRATCH/slurm-shims`,
   no "SmartSlurm" substring) and puts it first on `$PATH`, then
   `$SMARTSLURM_ROOT/bin`. (`ssbatch` finds the real `sbatch` with
   `which -a sbatch | grep -v SmartSlurm`, so a mock under a path containing
   "SmartSlurm" would be filtered out.)
4. Redirects `HOME` into `$SS_SCRATCH` so the learning DB, stats, and job records
   live under scratch — it **never** touches a real `~/.smartSlurm` or a real
   Slurm.

> These are environment-mutating side effects, appropriate for a disposable
> sandbox or CI runner. On a shared workstation, prefer a container or a box
> where `/bin/sh` is already bash.

---

## Simulating failures

The mock `sbatch`/`srun` inject failures so attribution can be tested:

```bash
MOCK_FAIL_MODE=oom      MOCK_FAIL_JOB=tokenize ss_run_pipeline "./pipeline.sh" "$opts"
MOCK_FAIL_MODE=timeout  MOCK_FAIL_JOB=align    ss_run_pipeline "./pipeline.sh" "$opts"
MOCK_FAIL_MODE=nonzero  MOCK_FAIL_JOB=count    ss_run_pipeline "./pipeline.sh" "$opts"
```

- `MOCK_FAIL_MODE` = `oom` | `timeout` | `nonzero` — forces the payload result /
  reported state.
- `MOCK_FAIL_JOB` = substring of the job name (`-J`); omit to affect all jobs.

Injected jobs are killed **before** `.success` is touched, so `cleanUp.sh` falls
through to the `sacct` state exactly as on a real cluster. See
`test/mockSlurm/README.md` for each command's contract and
`test/README.md` for the resource-model layer and the full block list.

---

## Gotchas

- **Slow by design** — run one block at a time in a time-boxed environment.
- **Scratch** — everything transient lives under `$SS_SCRATCH` (default
  `/tmp/ss-scratch`); delete it to reset. Nothing writes to a real Slurm or
  `~/.smartSlurm`.
- **`sh` vs `bash`** — the `sh`→`bash` repoint (above) is what makes the
  `#!/bin/sh` tools behave like production. `checkRunGrid`/`diagnose.sh` are
  `#!/bin/bash`.
- **`sacct` format is pinned to production Slurm (25.11.7).** If you move to a
  Slurm major version that changes `sacct`'s tabular output, update
  `test/mockSlurm/sacct` (see its README) — it is the most version-sensitive file.
- **Moving target** — SmartSlurm's `bin/` is under active development. When a tool
  changes, re-run the affected block(s) and update fixtures/asserts to the new
  intended behavior before relying on them as a regression net.
