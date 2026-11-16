# Packaging: the `slurmless` Claude Code plugin

This repo is installable as a Claude Code **plugin** named `slurmless`. The plugin
bundles two skills and the harness runtime they depend on, so a user can install it
and port a pipeline from **their own** directory without cloning anything.

## What ships

```
.claude-plugin/
  plugin.json            # manifest (name "slurmless", version, userConfig)
  marketplace.json       # single-plugin marketplace, source "./" (installs straight from this repo)
skills/
  migration-assistant/   # SKILL.md + tools/ + references/   -> /migration-assistant
  slurmless-harness/      # SKILL.md                           -> /slurmless-harness
  smartslurm-adversary/   # SKILL.md (adversarial review gate) -> /smartslurm-adversary
agents/
  smartslurm-adversary.md # independent adversarial reviewer (invokes the skill)
test/                    # the mock-Slurm harness the skills run (sourced as $SL/test/lib/harness.sh)
bootstrap.sh             # one-command "is the sandbox green?" check ($SL/bootstrap.sh)
```

The repo **root is the plugin root**, because the Migration Assistant hard-depends
on the harness (`bootstrap.sh`, `test/lib/harness.sh`, `test/mockSlurm/*`): the
skill and the runtime are one installable unit. There is no root `bin/`, so the
plugin is installable on claude.ai / Cowork as well as the CLI.

## Install

The plugin ships inside SmartSlurm at `slurmless/`. A marketplace manifest in a repo
subdirectory is **not** reachable by the `owner/repo` shorthand (verified against the
plugin docs 2026-11-12), so add it from a local SmartSlurm checkout (2.1.275+):

```text
/plugin marketplace add ./SmartSlurm/slurmless
/plugin install slurmless@slurmless-marketplace
```

Validate the manifest before publishing:

```bash
claude plugin validate .            # add --strict in CI
```

## Two consumption modes, one source of truth

Skill content lives **once**, at `skills/<name>/` (the plugin-standard location).
For standalone development — using a checkout of this repo as the agent's working
directory — `.claude/skills/<name>` are **symlinks** back to `../../skills/<name>`,
so project-skill auto-discovery still finds them. The same holds for agents:
`agents/<name>.md` is canonical, and `.claude/agents/<name>.md` is a dev symlink back
to it so project-agent auto-discovery finds it in a standalone checkout. The symlinks are only for in-repo
dev; they are never part of the installed plugin (which uses `skills/` directly), so
the "no symlinks escaping the plugin root" rule does not apply to them.

Do not add a second real copy under `.claude/skills/` — that would drift. Edit the
canonical files under `skills/`.

## Path resolution (`$SL`) — why the skill does not assume a CWD

When installed as a plugin, the harness is **not** in the user's working directory
(their CWD is their own pipeline). `SKILL.md` therefore resolves the harness root
into `$SL` first:

- plugin install → `$SL = ${CLAUDE_PLUGIN_ROOT}` (substituted to the plugin's
  absolute path in the skill body).
- standalone repo → `$SL` = the repo root.

Every harness/tool call is `"$SL/bootstrap.sh"`, `"$SL/test/lib/harness.sh"`,
`"$SL/skills/migration-assistant/tools/<tool>.sh"`. (`${CLAUDE_PLUGIN_ROOT}` is
substituted in Markdown but is **not** exported into Bash-tool shells, so the skill
resolves it once into `$SL` rather than relying on the env var in commands.)

## The SmartSlurm dependency

The harness runs SmartSlurm's **unmodified** `bin/runAsPipeline` etc., so it needs a
SmartSlurm checkout. `bootstrap.sh` resolves one in this order:

1. `$SMARTSLURM_ROOT` if set;
2. a vendored parent (`.../SmartSlurm/slurmless/` → the enclosing SmartSlurm);
3. a sibling `../SmartSlurm`;
4. clone `https://github.com/ld32/smartSlurm` (depth 1) next to the harness.

The plugin also exposes an optional `userConfig.smartslurm_path` (a directory). When
the user sets it, the skill uses `${user_config.smartslurm_path}` for
`SMARTSLURM_ROOT`; otherwise the resolution chain above applies. No network and no
config → point `SMARTSLURM_ROOT` at a checkout by hand.

## Scratch / state

The harness runs in a scratch workspace (`SS_SCRATCH`), never in the repo or plugin
root. Keep accidental in-root debris out of git via the root `.gitignore`. For a
plugin that needs durable state, `${CLAUDE_PLUGIN_DATA}` persists across updates —
the harness does not currently use it (each run is scratch).

## Agent-neutral distribution

The Claude plugin is one surface, not the only one. The **substance is already
agent-neutral** — only the outer trigger/distribution layer is Claude-specific:

| Layer | What it is | Neutral? |
|---|---|---|
| Harness (`test/`, mock Slurm, `bootstrap.sh`) | pure bash | yes |
| Tools (`extract_dag`, `emit_rap`, `compare_outputs`) | pure bash | yes |
| `SKILL.md` **body** + `references/` | Markdown procedure | yes |
| `SKILL.md` front-matter, `.claude-plugin/`, `/plugin install`, `.claude/skills/` | Claude Code trigger/dist | no |

Options for non-Claude agents (OpenAI / Codex, Gemini, Cursor, Aider, …), cheapest
first:

- **A — `AGENTS.md` convention (in place today).** `AGENTS.md` has a *"Using the
  skills with any agent"* section; the `SKILL.md` files are plain-Markdown playbooks.
  Frameworks that read `AGENTS.md`/conventions files pick it up automatically; others
  load the relevant `SKILL.md` as system/developer instructions. The agent sets
  `$SL` = repo root and calls the bash tools directly. Zero new code.
- **B — MCP server (recommended next for real cross-vendor tool-calling).** Wrap the
  three tools + `bootstrap` as an MCP server and expose `SKILL.md` as an MCP *prompt*;
  MCP is spoken by Claude, OpenAI, Gemini, and most frameworks, so one server serves
  every client. It must call the *same* `tools/*.sh` — no logic is duplicated.
- **C — plain `migrate` CLI + pasteable playbook.** A wrapper CLI + the Markdown as a
  prompt pack, for frameworks with only a generic shell tool.

**Keeping surfaces in sync — don't branch; wrap.** The anti-pattern is divergent
plugin/neutral branches. Instead: one branch, one neutral core (`test/`,
`skills/*/tools`, `skills/*/references`, the `SKILL.md` bodies), and each surface is a
*thin adapter* beside it — `.claude-plugin/` (manifest only), a future `mcp/` server
(calls `tools/`), `AGENTS.md` (prose pointer). None re-implement logic, so there is
nothing to reconcile. The guard is **CI, not discipline**: the `ma_*` blocks and
`skills/migration-assistant/examples/ma_dryrun.sh` exercise the tools the same way
every adapter does, so drift breaks the build. If a surface ever needs a *generated*
artifact (e.g. an OpenAI function-schema JSON), generate it from the core in CI;
never hand-maintain it. A one-way distribution branch (plugin-only files) is fine if
produced from `main` by CI (`git archive`/subtree), never hand-edited.

## Looking ahead: vendoring into SmartSlurm

When slurmless moves to `SmartSlurm/slurmless/` and this standalone repo is retired
(per `CLAUDE.md`), the **runtime design is already forward-compatible** — the `$SL`
indirection and `bootstrap.sh`'s resolver chain were chosen for exactly this. No
`SKILL.md` change is needed. The open decisions are about *distribution*:

- **Plugin root = the vendored `SmartSlurm/slurmless/` subdir, not `SmartSlurm/`
  root.** Two reasons: (a) a plugin that contains a root `bin/` is reportedly **not
  installable on claude.ai / Cowork**, and SmartSlurm has `bin/`; (b) shipping all of
  SmartSlurm as "the migration plugin" conflates the tool with the skill. Keep
  `.claude-plugin/` inside `slurmless/`.
- **Distribute via a `git-subdir` marketplace source** (`url: ntc-agnt/SmartSlurm`,
  `path: slurmless`) so install carves out just the `slurmless/` subtree as the
  plugin. The standalone repo retires as the *dev* source; the *published* plugin is
  built from the subdir — same "generate, don't fork" rule as above.
- **SmartSlurm becomes "the parent," resolved at runtime.** Used in-place inside a
  SmartSlurm checkout, `bootstrap.sh`'s "walk up for `bin/runAsPipeline`" resolves the
  enclosing SmartSlurm (the vendored case is already tested). From a cache install
  (just the subdir, no parent present) it falls to `SMARTSLURM_ROOT` /
  `userConfig.smartslurm_path` / clone — already handled.
- **Top-level discovery for SmartSlurm devs.** Add a `SmartSlurm/.claude/skills/*`
  symlink into `slurmless/skills/*` (same pattern we use one level up) so an agent
  working at the SmartSlurm root still auto-discovers the skills.

Net new work at vendoring time: the marketplace `git-subdir` entry and the top-level
symlink — not a repackage.

> **Confirm before relying on these two:** the root-`bin/` claude.ai/Cowork install
> restriction and the `git-subdir` marketplace source come from a plugin-docs research
> pass, not from a run in this repo. Verify both against the current Claude Code plugin
> docs before committing to this distribution shape.
>
> **Verified 2026-11-12 (against the plugin docs):** the root-`bin/` claude.ai/Cowork
> restriction holds but does **not** bite here — the plugin root is `slurmless/`, which
> has no `bin/`. However, the **`git-subdir` *marketplace* source does not exist**:
> `git-subdir` is a plugin-*entry* source, not a marketplace locator, and a marketplace
> manifest in a subdirectory is **not** reachable by the `owner/repo` shorthand. So a
> remote one-line install needs a **root-level `.claude-plugin/marketplace.json` in
> SmartSlurm** (or an `extraKnownMarketplaces` settings entry). Until that exists,
> install from a local checkout (see **Install** above). Setting up the remote path is
> a distribution decision for the maintainers.
