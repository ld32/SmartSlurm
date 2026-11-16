# CLAUDE.md

You are working in **slurmless**, a Slurm-free test harness for SmartSlurm.

**Read [`AGENTS.md`](AGENTS.md) first** — it is the model-neutral session-startup
guide (what this is, the golden-path commands, how to add a test, failure
injection, and the drift model). This file only restates the rules most likely
to trip you up.

## Critical rules

- **Do not modify SmartSlurm's `bin/` to make a test pass.** The harness runs the
  *unmodified* tools; a failing test that reveals a real tool bug is a finding to
  report, not a tool to quietly patch. Change a tool only when the task is
  explicitly to change that tool — then update the affected fixtures/asserts.
- **Prove the baseline first:** `./bootstrap.sh` (runs `clean_2step`). If it is
  red, the environment is not set up — stop and diagnose before other work.
- **Develop against one block** (`./test/run_tests.sh <block>` runs serial with live
  output); a bare `./test/run_tests.sh` runs the whole suite **in parallel** (one
  worker per core, capped at 8), so it's no longer slow — use `-j1` to force serial.
- **Scratch only.** Never aim the harness at a real `~/.smartSlurm` or a real
  cluster.
- **The harness is code too — self-test it; don't trust its own green.** A change to
  the test *runner* or its orchestration (the `-jN` dispatcher, `ss_parallel_run`,
  scratch/RAM selection, `ss_activate`) is **not** validated by "the suite still
  passes" — that signal is produced *by* the thing you changed. It needs its own
  block asserting its invariants **and** its own adversary pass, and you watch the
  **process tree**, not just the `RESULT:` line. (The parallel default once shipped a
  fork-recursion that "passed" exactly this way.) See AGENTS.md → *Hard rules* and the
  adversary skill's third bug class.
- **Source of truth: `SmartSlurm/slurmless/` (in-tree).** slurmless now lives in
  SmartSlurm as a plain subdirectory and is developed there directly, versioned with
  the engine — one commit, so the harness can't drift from the tools it tests. The
  former standalone repo (`slurmless_smartSlurm_testing`) is **archived**, and the old
  "develop standalone, then re-snapshot into SmartSlurm" workflow is **retired** — do
  not re-snapshot.

## Working with the SmartSlurm upstream (ld32 directives)

These rules come from ld32, who owns `master` on the upstream SmartSlurm repo.
They override anything above when they conflict.

- **Never touch `master` yourself.** Do **not** open a PR against `master` and do
  **not** push to it. We may create new branches on the main repo, but `master`
  is ld32's — they push to it directly. Land our work on branches and let ld32
  integrate.
- **Always pull before you push.** ld32 pushes to `master` on their own cadence,
  so any branch we build on can be stale. `git fetch` + rebase/merge the latest
  base **before** every push to avoid clobbering or diverging from their changes.
- **Priority is resource efficiency, not UI/UX.** Stop investing effort in
  interface/experience polish. Direct work toward improving SmartSlurm's resource
  efficiency (estimation accuracy, allocation tightness, wasted-reservation
  reduction).

### Current program state

- **jobRecords-v2 migration: wrapping up.** We are closing out the in-flight jR2
  correctness work and not starting new jR2 scope.
- **Feature development: parked.** Speculative features (e.g. the `#@pipeline`
  name tag, guided legacy curation) are built and testable on parallel branches
  but are **not** being pushed forward or merged for now. Resume only on explicit
  direction.

## Skills

Two skills live under `skills/` (the plugin-standard location; `.claude/skills/`
holds dev symlinks to them for in-repo auto-discovery — see `PACKAGING.md`):

- `skills/slurmless-harness/SKILL.md` — run/develop with this harness. Invoke in
  Claude Code with `/slurmless-harness`.
- `skills/migration-assistant/SKILL.md` — convert a pipeline to runAsPipeline and
  validate it in the sandbox. Invoke with `/migration-assistant`.

Both are usable as plain Markdown by any agent. This repo is installable as the
`slurmless` Claude Code plugin (manifest in `.claude-plugin/`); see `PACKAGING.md`.
