[&larr; SmartSlurm README](../README.md) &middot; [Docs index](README.md)

# Installation and setup

Setup is two commands — clone, then run the setup script. Nothing is compiled or
copied into system directories, and **you never add an `export` or edit a file
yourself** (avoiding that manual `PATH` export is the whole point): the script
edits your `~/.bashrc` for you.

```bash
# 1. Download SmartSlurm
git clone https://github.com/ld32/SmartSlurm.git $HOME/SmartSlurm
cd $HOME/SmartSlurm

# 2. Run the setup script once — it adds the `smartslurm` command to your ~/.bashrc
bash smartSlurmSetup.sh

# 3. Open a new shell (or: source ~/.bashrc), then turn SmartSlurm on
smartslurm activate
```

`smartSlurmSetup.sh` is the whole install step. It appends a small, clearly
marked block to your `~/.bashrc` — you do **not** write this by hand:

```bash
# >>> smartSlurm initialize >>>
export SMARTSLURM_HOME="/home/you/SmartSlurm"
[ -f "$SMARTSLURM_HOME/bin/smartslurm.sh" ] && . "$SMARTSLURM_HOME/bin/smartslurm.sh"
# <<< smartSlurm initialize <<<
```

The script writes that block (and `smartslurm uninstall` removes it again). The
block only teaches your shell the `smartslurm` command; it does **not** put
SmartSlurm on your `PATH` until you run `smartslurm activate`. Re-running
`bash smartSlurmSetup.sh` is safe — it updates the block in place (backing up the
rc first) instead of adding a second copy. All the logic lives in
`bin/smartslurm.sh`, so a later `git pull` upgrades behavior without rewriting
your rc file.

## The `smartslurm` command

| Command | What it does |
|---------|--------------|
| `smartslurm activate` | Puts SmartSlurm on your `PATH` and **aliases `sbatch` → smart `ssbatch`**, so existing scripts get automatic resource estimation for free. It prints a notice that the alias is in effect (previously this was silent) and prefixes your prompt with `{smartSlurm} ` so you can always see when SmartSlurm is active. |
| `smartslurm deactivate` | Removes **only** SmartSlurm's own `PATH` entries and prompt prefix, restoring the system `sbatch`. Anything you added *after* activating — a `conda activate`, a `module load`, your own exports — is left untouched. |
| `smartslurm status` | Reports whether SmartSlurm is active and whether `sbatch` currently resolves to the smart alias or the system binary. |
| `smartslurm uninstall` | Removes the initialize block from your rc file (with a backup). New shells will no longer know the command. |

> [!NOTE]
> The prompt prefix uses curly braces — `{smartSlurm} ` — deliberately distinct
> from conda's `(env) `. The two coexist cleanly: activating or deactivating one
> never disturbs the other's prompt or `PATH`, in any order.

> [!TIP]
> Want SmartSlurm's tools without aliasing `sbatch` (for example, only inside one
> pipeline)? Skip `activate` and either call `ssbatch` / `runAsPipeline` by full
> path, or add just `$HOME/SmartSlurm/bin` to your `PATH` yourself — that exposes
> `ssbatch` while leaving the system `sbatch` in place. `smartslurm activate` is
> the convenient, fully reversible default.

## Turning SmartSlurm on and off

`smartslurm activate` turns SmartSlurm on; `smartslurm deactivate` restores the
system `sbatch` in the current shell. Deactivation removes only SmartSlurm's own
`PATH` entries and prompt prefix, so anything you added after activating is
preserved:

```bash
smartslurm deactivate
```

> [!NOTE]
> **`unExport` is deprecated.** The old `source unExport; unExport` still works
> as a thin shim that warns and delegates to `smartslurm deactivate`, but it will
> be removed in a future release — switch to `smartslurm deactivate`. (The legacy
> `unExport` stripped whatever directory held the current `sbatch`, which could
> remove the *system* `sbatch` from your `PATH` when SmartSlurm was not actually
> aliasing it; `smartslurm deactivate` does not have that failure mode.)

## Optional: Graphviz for workflow charts

Only needed for the workflow-chart (`w`/`g`) option in `checkRun`:

```bash
module load conda/miniforge3/24.11.3-0
mamba create -n smartSlurmEnv -c conda-forge graphviz
```

> [!NOTE]
> The `module load conda/...` line is specific to an HPC that provides conda as an environment module. On other systems, activate conda however your site does. If your site uses a different conda module name, adjust it here and in `checkRun`.

---


---

[&larr; SmartSlurm README](../README.md) &middot; [Docs index](README.md)
