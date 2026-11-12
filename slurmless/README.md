# slurmless

**A Slurm-free test harness for [SmartSlurm](https://github.com/ld32/smartSlurm).**

slurmless lets you run SmartSlurm's **unmodified** tools — `runAsPipeline`,
`ssbatch`, `cleanUp.sh`, `checkRun`, `checkRunGrid`, `diagnose.sh` — on a machine
with **no Slurm scheduler** (a laptop, a CI runner, or an AI-agent sandbox), while
still executing real job payloads. That makes it possible to develop and
regression-test SmartSlurm — pipeline conversion, dependency wiring, failure
attribution, the status grid, the rerun menu, and the resource-estimation
learning loop — entirely off-cluster.

It works by interposing **mock Slurm commands** (`sbatch`, `sacct`, `srun`, …) on
`$PATH` ahead of the real binaries. SmartSlurm is **not modified**; the mocks
actually execute each job's payload, so the real `cleanUp.sh` trap and the real
`.success`/`.failed` flags are produced, exactly as on a cluster.

## Who it's for

- **SmartSlurm developers** — a fast, Slurm-free inner loop and regression net for
  work on `runAsPipeline`, resource estimation, logging, and the monitor UI.
- **Pipeline builders (bioinformaticians)** — a safe place to convert a hand-rolled
  SLURM workflow (or loose scripts) to `runAsPipeline` and prove it behaves before
  running it on real HPC, with the **Migration Assistant** (below) driving the
  conversion.
- **AI agents** — a self-contained environment to develop/test SmartSlurm and to
  run the Migration Assistant. Two skills ship under `skills/` and the repo installs
  as a Claude Code plugin; the same skills are usable by any agent via `AGENTS.md`
  (see [`skills/README.md`](skills/README.md)).

## Quick start

slurmless is vendored into SmartSlurm at `SmartSlurm/slurmless/`, so from a
SmartSlurm checkout:

```bash
cd SmartSlurm
./slurmless/test/run_tests.sh --list          # see every test block and what it checks
./slurmless/test/run_tests.sh clean_2step     # run the baseline block (must pass)
```

To run the standalone repo against a SmartSlurm checkout elsewhere:

```bash
export SMARTSLURM_ROOT=/path/to/SmartSlurm
/path/to/slurmless/test/run_tests.sh clean_2step
```

A single block is usually quick; the full suite is slow (the mock scheduler
sleeps through SmartSlurm's cleanup counters), so pass a block name in a
time-boxed environment.

## Migration Assistant

On top of testing SmartSlurm's own tools, slurmless now helps you **move a pipeline
onto SmartSlurm**. The **Migration Assistant (MA)** converts a pipeline to
`runAsPipeline` (rAP) and proves the port off-cluster in this sandbox before you
touch a cluster. It handles two starting points:

- **An existing SLURM workflow** (hand-rolled `sbatch` + `--dependency`) → MA
  recovers its job DAG, ports it to rAP, and **differentially tests** it (runs your
  original and the port side-by-side and confirms identical outputs under an
  identity contract you approve).
- **Loose scripts with no pipeline** → MA designs the DAG with you and writes the
  annotated rAP script.

MA converts **faithfully** (it never silently rewrites your logic; bugs are
preserved and *flagged* for your decision, light-touch) and self-tests the port in
a clean sandbox before handing it back. See
[`skills/migration-assistant/README.md`](skills/migration-assistant/README.md).

### What it adds to the harness

Two layers were added, both Slurm-free and jR-v2-independent:

- **User-workflow mode** — the mock `sbatch` can *register-don't-run* an
  **unmodified** hand-rolled submit script (so a real `jid=$(sbatch …)` +
  `--dependency=…` builds a true DAG), then a `wfdrain` executor runs that DAG
  topologically. It emulates the **full `--dependency` grammar** (after / afterany /
  afterok / afternotok / afterburstbuffer / aftercorr / singleton, `,`=AND / `?`=OR,
  `+time`, `--kill-on-invalid-dep`), job **arrays** (`--array`, element-wise
  `aftercorr`), `sbatch --wait`, and live `scancel`/`squeue` state — enabled behind
  `ss_wf_activate` so the default (SmartSlurm-testing) mode is unchanged. See
  [`test/README.md`](test/README.md#user-workflow-mode-migration-assistant).
- **Conversion tools** — `extract_dag`, `emit_rap`, `compare_outputs` under
  `skills/migration-assistant/tools/`, plus the `ma_*` test blocks and an
  agent-in-the-loop fuzz eval (`eval/migration/`).

Verify it end-to-end in one command (converts a tiny SLURM workflow and
differentially tests the port; exits 0 on success — works from a checkout or an
installed plugin):

```bash
bash skills/migration-assistant/examples/ma_dryrun.sh
```

### Install as a plugin

MA is packaged as a Claude Code plugin (named `slurmless`) that bundles both the
`migration-assistant` and `slurmless-harness` skills **and** the harness they run
on. Install it (Claude Code 2.1.275+):

```text
/plugin install slurmless --marketplace ntc-agnt/slurmless_smartSlurm_testing
```

or in two steps on any version:

```text
/plugin marketplace add ntc-agnt/slurmless_smartSlurm_testing
/plugin install slurmless@slurmless-marketplace
```

During install you can optionally set a **SmartSlurm checkout** path; if you leave
it blank the harness auto-resolves a vendored/sibling checkout or clones upstream
on first `bootstrap.sh`. Then point the agent at your pipeline and ask it to port
and test — the skill triggers automatically, or invoke it as `/migration-assistant`.

**Not just Claude.** The skills are plain-Markdown playbooks; any agent
(OpenAI / Codex, Gemini, Cursor, Aider, …) or a human can run the same workflow
with no plugin — see the *"Using the skills with any agent"* section of
[`AGENTS.md`](AGENTS.md). For the full user guide see
[`skills/migration-assistant/README.md`](skills/migration-assistant/README.md); for
how the plugin is built, kept agent-neutral, and will move into SmartSlurm, see
[`PACKAGING.md`](PACKAGING.md).

## Layout

| Path | What it is |
|------|------------|
| `test/run_tests.sh` | The pipeline + resource-model test suite (CI exit code). |
| `test/checkrun_ui_tests.sh` | The interactive-monitor (`checkRun`/`checkRunGrid`) UI suite, driven over a pseudo-terminal. |
| `test/lib/harness.sh` | Environment activation, the resource-model API, and assertion helpers. |
| `test/mockSlurm/` | The mock Slurm commands (`sbatch`, `sacct`, `srun`, `squeue`, …). |
| `test/fixtures/` | Tiny annotated pipelines exercising specific behaviors. |
| `test/models/` | Resource-curve models for the estimator learning-loop tests. |
| `test/version.txt` | The SmartSlurm baseline this harness was last validated against. |
| `test/mockSlurm/wfdrain` | The user-workflow DAG executor (register-don't-run mode). |
| `skills/` | The two shipped skills — `migration-assistant` and `slurmless-harness` (see [`skills/README.md`](skills/README.md)). |
| `.claude-plugin/` | Plugin + marketplace manifests that make this repo installable. |
| `eval/migration/` | Agent-in-the-loop fuzz eval for the conversion workflow. |
| `AGENTS.md` / `PACKAGING.md` | Agent-neutral entry guide / plugin-build + vendoring design. |

## Learn more

- **[skills/README.md](skills/README.md)** — the two skills, what each does, and how
  to use them (Claude plugin **or** any agent).
- **[CONTEXT.md](CONTEXT.md)** — setup, how it stays in sync with SmartSlurm, the
  harness API, failure injection, and gotchas.
- **[test/README.md](test/README.md)** — the full block catalog, the user-workflow
  mode, the resource-model layer, fidelity tiers, and how to add a fixture.
- **[test/mockSlurm/README.md](test/mockSlurm/README.md)** — the per-command mock
  contracts (default **and** user-workflow mode).
- **[PACKAGING.md](PACKAGING.md)** — how the plugin is built, the agent-neutral
  options and sync model, and the move into SmartSlurm.
