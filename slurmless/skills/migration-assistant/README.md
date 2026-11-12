# Migration Assistant — user guide

**What it does.** The Migration Assistant (MA) helps you move a compute pipeline onto
SmartSlurm's `runAsPipeline` (rAP) framework — and proves it works **off-cluster**,
in the Slurm-free *slurmless* sandbox, before you ever run it for real. It handles
two starting points:

- **An existing SLURM workflow** (hand-rolled `sbatch` + `--dependency` scripts) → it
  recovers your job DAG, ports it to rAP, and **differentially tests** it (runs your
  original and the port side-by-side and confirms identical outputs).
- **Loose scripts with no pipeline yet** → it designs the pipeline with you and builds
  the rAP script.

**What it will NOT do.** It converts *faithfully* — it never silently rewrites your
logic. If it spots an obvious, high-confidence error (e.g. mapping reads before
adapter trimming) it will ask, never "fix" things behind your back, and record your
decision. Your methods are assumed intentional; novel approaches are respected.

## Install

The fastest way to get the skill **and** the sandbox it needs is to install the
`slurmless` plugin (Claude Code 2.1.275+):

```text
/plugin install slurmless --marketplace ntc-agnt/slurmless_smartSlurm_testing
```

(or `/plugin marketplace add ntc-agnt/slurmless_smartSlurm_testing` then
`/plugin install slurmless@slurmless-marketplace`). The plugin bundles the harness,
so once installed you can be in **your own pipeline's directory** — you do not need
to clone this repo. Optionally set the **SmartSlurm checkout** path at install time;
leave it blank to auto-resolve or clone upstream on first run.

You can also use it without installing, by working inside a checkout of this repo
(or `SmartSlurm/slurmless/`) as your agent's working directory — the skill is
auto-discovered from `.claude/skills/`.

**Verify it works** — run the bundled smoke test, which converts a tiny SLURM
workflow to rAP and differentially tests it (exits 0 on MATCH). It self-locates the
harness, so the same command works from a repo checkout or an installed plugin:

```bash
bash skills/migration-assistant/examples/ma_dryrun.sh        # from a repo checkout
# or, installed as a plugin:
bash "${CLAUDE_PLUGIN_ROOT}/skills/migration-assistant/examples/ma_dryrun.sh"
```

## Prerequisites

1. The `slurmless` plugin installed (above), **or** a checkout of this repository
   (the slurmless harness + the skill), or SmartSlurm with `slurmless/` vendored in.
2. A **SmartSlurm checkout** the harness can point at (`SMARTSLURM_ROOT`; auto-resolved
   when vendored).
3. Your pipeline: the scripts (and, for a SLURM workflow, the submit script).
4. An **environment spec** so the real tools exist in the sandbox (a conda `env.yml`,
   a container, or a module list), plus either small real inputs to downsample or a
   note of each step's input format so tiny inputs can be synthesized.

## How to use it

1. **Open an AI coding agent (e.g. Claude Code)** with the `slurmless` plugin
   installed and your pipeline's directory as the working directory — or, without
   the plugin, with a checkout of this repo as the working directory.
2. **Prove the sandbox is set up:** have the agent run `bootstrap.sh` (it resolves
   the harness root; runs the `clean_2step` baseline). Green → ready.
3. **Ask the agent to convert your pipeline** — e.g. *"port this SLURM workflow to
   runAsPipeline and test it,"* pointing it at your scripts. The `migration-assistant`
   skill triggers automatically; you can also invoke it explicitly as
   `/migration-assistant`.
4. **Collaborate at the checkpoints.** The agent pauses for your sign-off at: the
   pipeline plan / DAG, the differential **identity contract** (which outputs count +
   how to normalize volatile fields), and the final review. It asks whenever a
   conversion decision is genuinely ambiguous.
5. **Receive a tested pipeline.** Before handing back, the agent self-tests the port
   in a clean sandbox (and, ideally, has an independent agent re-run it against your
   original). You get the rAP `pipeline.sh`, a `MIGRATION_NOTES.md` recording any
   flagged issues and your decisions, and the test results.

## What you get out of SmartSlurm/rAP

Resource learning (auto mem/time estimation), failure floors and auto-requeue,
checkRun/resume, and consolidated logging — without hand-wiring `--dependency`.

## Scope & limits

- The sandbox validates **orchestration and data flow** faithfully; it runs your real
  tools on tiny inputs where available, or stubs a step whose tool/input is missing
  (wiring is validated, not the science of that step).
- rAP's dependency model is `afterok`-style (per-step, with fan-in). Richer SLURM
  dependencies (`afternotok` cleanups, `?`-OR, `after`+time, `aftercorr`) don't map
  1:1; the agent surfaces these and proposes how to express them, rather than guessing.
- Scratch only — the harness never touches a real `~/.smartSlurm` or a real cluster.

See `SKILL.md` for the agent-facing workflow and `references/` for the rAP marker,
plan, and identity-contract formats. `eval/migration/` holds the agent-in-the-loop
fuzz tests that validate the conversion workflow itself.
