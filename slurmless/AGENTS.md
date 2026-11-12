# AGENTS.md — working in slurmless

Model-neutral session-startup guide for any coding agent (Anthropic, OpenAI,
Cursor, …) operating in this directory. Humans: this doubles as a quick operator
guide; see `CONTEXT.md` and `test/README.md` for depth.

## What this is

**slurmless** is a **Slurm-free test harness** for
[SmartSlurm](https://github.com/ld32/smartSlurm). It runs SmartSlurm's
**unmodified** tools (`runAsPipeline`, `ssbatch`, `cleanUp.sh`, `checkRun`,
`checkRunGrid`, `diagnose.sh`) with **no Slurm scheduler**, by interposing mock
Slurm commands (`sbatch`, `sacct`, `srun`, …) on `$PATH`. The mocks actually run
each job's payload, so real `.success`/`.failed` flags and logs are produced —
failure attribution behaves as on a cluster.

It is vendored into SmartSlurm at `SmartSlurm/slurmless/` as a plain subdirectory
(a point-in-time snapshot of this standalone repo, **not** a `git subtree`). Until
the jobRecords-v2 migration lands, develop here and re-snapshot into SmartSlurm;
afterward, development moves into `SmartSlurm/slurmless/` and this repo is archived.

## Session startup (do this first)

1. **Confirm the SmartSlurm checkout.** Vendored (`SmartSlurm/slurmless/`):
   `SMARTSLURM_ROOT` auto-resolves. Standalone: `export SMARTSLURM_ROOT=/path/to/SmartSlurm`.
2. **Do not modify SmartSlurm's `bin/` to make a test pass.** The harness tests
   the *unmodified* tools. If a test reveals a tool bug, report it; change a tool
   only when the task explicitly is to change that tool (then update the affected
   fixtures/asserts).
3. **Prove the baseline before anything else:** `./bootstrap.sh` (runs the
   `clean_2step` block). If that is red, stop and diagnose — the environment is
   not set up.
4. **Run one block at a time.** The mock scheduler sleeps through SmartSlurm's
   cleanup counters, so the full suite is slow; a single block is usually quick.
   Curve-fit blocks (`scaling_*`, `noise_fit`, `contaminate_*`) need `gnuplot`.

## Golden-path commands

```bash
./bootstrap.sh                         # one-command: resolve SmartSlurm + run baseline
./test/run_tests.sh --list             # every block and what it verifies
./test/run_tests.sh <block>            # run one block (substring match)
./test/checkrun_ui_tests.sh            # the interactive-monitor UI suite (needs python3)
./test/mockSlurm/check_sacct_format.sh # sacct-format drift check (CI-friendly)
```

Interactive:

```bash
source test/lib/harness.sh && ss_activate
ss_new_workspace /tmp/play
cp test/fixtures/clean_2step/* .
ss_run_pipeline "./pipeline.sh 777" "sbatch -p short -c 1 --mem 4G -t 50:0"
```

## What `ss_activate` changes (and how to opt out)

It is built for a disposable sandbox / CI runner. It will: install `gnuplot` if
missing, repoint `/bin/sh → bash` (SmartSlurm's `#!/bin/sh` scripts need bash;
production RHEL `/bin/sh` *is* bash), stage the mocks on `$PATH`, and redirect
`HOME` into scratch (so nothing touches a real `~/.smartSlurm`). On a shared
host, opt out: `SS_NO_INSTALL=1`, `SS_NO_SH_REPOINT=1`. It prints a one-line
summary of what it touched.

## Adding a test

1. `mkdir test/fixtures/<name>/` with a `pipeline.sh` (annotated `#@` script) and
   any step tools it calls by bare name.
2. Add a block to `test/run_tests.sh`: `_begin <name>`, `ss_new_workspace`, copy
   the fixture, `ss_run_pipeline`, assert with the helpers (`assert_file`,
   `assert_nonempty`, `assert_contains`, `assert_valid_bash`, `assert_rc`, …).
3. Add `<name>` to `ALL_BLOCKS` and give it an `ss_purpose` line. Keep payloads
   tiny and deterministic.

## Simulating failures

```bash
MOCK_FAIL_MODE=oom|timeout|nonzero  MOCK_FAIL_JOB=<jobName-substring> ss_run_pipeline ...
```
Injected jobs are killed before `.success` is touched, so `cleanUp.sh` falls
through to the mock `sacct` state exactly as on a real cluster.

## User-workflow mode (Migration Assistant)

Besides testing SmartSlurm's own tools, the harness can run an **unmodified**,
hand-rolled SLURM workflow (stock `sbatch` + dependencies) invisibly — the basis
of the Migration Assistant. Enable it with `ss_wf_activate`; the mock `sbatch`
then REGISTERS each job instead of running it, so a script's `jid=$(sbatch ...)`
capture + `--dependency=...` build a real DAG that `ss_wf_drain` executes
topologically.

```bash
source test/lib/harness.sh && ss_activate
ss_new_workspace /tmp/play && cp test/fixtures/wf_pipeline/* .
ss_wf_activate
bash ./submit.sh      # registers jobs (no payload runs yet)
ss_wf_drain           # executes the DAG in dependency order
# or:  ss_run_workflow "./submit.sh"   (activate + run + drain)
```

Supported: the full `--dependency` grammar (after / afterany / afterok /
afternotok / afterburstbuffer / aftercorr / singleton, `,`=AND / `?`=OR, `+time`,
`--kill-on-invalid-dep`); job arrays (`--array`, with element-wise `aftercorr`);
`sbatch --wait`; and `scancel` / `squeue` live state. Everything is gated behind
`SS_WF_MODE` (set by `ss_wf_activate`), so default mode is unchanged.

Knobs:
- `MOCK_FAIL_JOB=<name-substr> MOCK_FAIL_MODE=oom|timeout|nonzero` — force a
  workflow step to fail (resource-sim), to exercise failure handling.
- `SS_WF_STUB="<name-glob> ..."` — run matching steps as success stubs (tiered
  execution) when the real tool/input isn't available; `SS_WF_STUB_MISSING=1`
  auto-stubs a step whose leading command is missing. A stub validates wiring,
  not data (it produces no outputs).
- `SS_WF_KILL_INVALID_DEP=yes|no` — default disposition of a job whose dependency
  can never be satisfied (kill vs hold); per-job override via
  `sbatch --kill-on-invalid-dep=`.

Blocks: `wf_afterok wf_depgrammar wf_arrays wf_wait wf_scancel wf_inject
wf_tiered wf_pipeline` — run them all with the `workflow` group
(`./test/run_tests.sh -j8 workflow`).

## Using the skills with any agent

Two task playbooks live under `skills/` as **plain Markdown** — usable by any agent
(or human), not just Claude Code:

- `skills/slurmless-harness/SKILL.md` — run/develop with this harness.
- `skills/migration-assistant/SKILL.md` — convert a hand-rolled SLURM workflow or
  loose scripts to SmartSlurm's `runAsPipeline`, validated here before a real cluster.

Claude Code auto-discovers them (and installs this repo as the `slurmless` plugin —
see `PACKAGING.md`). **Any other agent** (OpenAI / Codex, Gemini, Cursor, Aider, …)
or a human uses the exact same files, with no plugin:

1. Make this repo (or `SmartSlurm/slurmless/` when vendored) the working directory
   and set `$SL` to its root. Frameworks that read `AGENTS.md` / conventions files
   pick this guide up automatically; otherwise load the relevant `SKILL.md` as your
   system / developer instructions.
2. Follow that `SKILL.md` body — it is model-neutral prose. Ignore only the
   Claude-specific framing: the YAML front-matter, the `/migration-assistant`
   invocation, and `${CLAUDE_PLUGIN_ROOT}` (for you, `$SL` is just this repo root).
3. The tools it calls are plain bash:
   `skills/migration-assistant/tools/{extract_dag,emit_rap,compare_outputs}.sh`.
4. Verify end-to-end in one command (converts a tiny SLURM workflow and
   differentially tests the port; exits 0 on success):
   `bash skills/migration-assistant/examples/ma_dryrun.sh`.

The harness, tools, and `SKILL.md` procedures are the single source of truth; the
Claude plugin and this neutral entry are thin wrappers over them, so there is nothing
to keep in sync between them.

## Drift — the thing to protect

- **Tools vs harness:** vendored, they share one commit (can't drift); run the
  harness before pushing a tool change (`./bootstrap.sh`, then the relevant
  blocks). `test/version.txt` records the validated baseline for standalone runs.
- **sacct format** is pinned to production Slurm (25.11.7); it is the most
  version-sensitive file. If a Slurm upgrade changes `sacct` output, run
  `test/mockSlurm/check_sacct_format.sh` on the cluster and update
  `test/mockSlurm/sacct`.
- **Trap-line coupling:** the resource model reads the generated `cleanUp.sh`
  trap line positionally; the mock `sbatch` now warns loudly if that shape
  changes. Heed the warning rather than trusting a modeled result.

## Hard rules

- Never edit SmartSlurm `bin/` to force a green test (see startup #2).
- Never point the harness at a real `~/.smartSlurm` or a real cluster; scratch
  only (`$SS_SCRATCH`, default `/tmp/ss-scratch`).
- The vendored copy in SmartSlurm is a plain-subdirectory snapshot of this repo
  (not a `git subtree`): until the v2 migration lands, develop here and re-snapshot
  into `SmartSlurm/slurmless/`; afterward, develop there directly.
