# skills/

Two task **skills** ship here. Each is a self-contained, plain-Markdown playbook
(`SKILL.md`) plus any `tools/`, `references/`, and `examples/` it needs. They are
usable three ways, from the same files:

- **Claude Code plugin** — this repo installs as the `slurmless` plugin; the skills
  auto-discover and are invocable as `/slurmless-harness` and `/migration-assistant`
  (see [`../PACKAGING.md`](../PACKAGING.md)).
- **In-repo** — working inside a checkout (or `SmartSlurm/slurmless/` once vendored),
  Claude Code auto-discovers them via the `.claude/skills/` dev symlinks.
- **Any other agent or a human** — follow the `SKILL.md` body as plain Markdown; the
  tools are plain bash. See the *"Using the skills with any agent"* section of
  [`../AGENTS.md`](../AGENTS.md).

> One source of truth: the canonical skill content lives here under `skills/`.
> `.claude/skills/<name>` are dev symlinks back to it (never shipped in the plugin).
> Do not keep a second copy — edit the files here.

## `slurmless-harness`

Run or develop with the slurmless Slurm-free SmartSlurm test harness — the mock-Slurm
suite, adding a fixture, reproducing an OOM/timeout/failure, resource estimation.
It's the task-oriented front door to the harness documented in
[`../AGENTS.md`](../AGENTS.md), [`../CONTEXT.md`](../CONTEXT.md), and
[`../test/README.md`](../test/README.md).

- Entry: [`slurmless-harness/SKILL.md`](slurmless-harness/SKILL.md)
- First move (from a repo checkout): `./bootstrap.sh` (green `clean_2step` ⇒ the
  sandbox is set up).

## `migration-assistant`

Convert a non-SmartSlurm pipeline to `runAsPipeline` (rAP) and validate it in the
sandbox before a real cluster — either an existing hand-rolled SLURM workflow
(`sbatch` + `--dependency`) or a pile of loose scripts. It recovers the DAG, ports
it, and **differentially tests** the port (original vs rAP, identical outputs under
an identity contract). It converts **faithfully** (never silently rewrites logic;
flags suspected bugs light-touch and asks) and self-tests before handback.

- User guide: [`migration-assistant/README.md`](migration-assistant/README.md)
- Agent workflow: [`migration-assistant/SKILL.md`](migration-assistant/SKILL.md)
- Tools (plain bash): `migration-assistant/tools/` — `extract_dag.sh`,
  `emit_rap.sh`, `compare_outputs.sh`
- Formats: `migration-assistant/references/` — rAP markers, the logical-plan format,
  the identity contract, the project-plan template
- Smoke test / worked example: `migration-assistant/examples/ma_dryrun.sh`
  (converts a tiny SLURM workflow and differentially tests it; exits 0 on MATCH;
  runs from a checkout or an installed plugin)

It builds on the harness's **user-workflow mode** (`ss_wf_activate` / `ss_wf_drain`);
see [`../test/README.md`](../test/README.md#user-workflow-mode-migration-assistant)
for that layer and [`../eval/migration/`](../eval/migration/) for the
agent-in-the-loop fuzz eval that validates the conversion workflow itself.

## Prerequisites (both skills)

A SmartSlurm checkout whose `bin/` holds the tools under test (auto-resolved when
vendored or as a sibling `../SmartSlurm`, else set `SMARTSLURM_ROOT`); `gnuplot` for
the curve-fit blocks; and a POSIX/bash environment where `/bin/sh` is bash (the
harness repoints it in a disposable sandbox). The harness runs in scratch
(`$SS_SCRATCH`) — never a real `~/.smartSlurm` or a real cluster.
