#!/bin/bash
# smartSlurmSetup.sh -- one-time setup for the `smartslurm` shell command.
#
# This script does all the work: it adds a small "initialize" block to your
# shell rc file (~/.bashrc) so future shells know the `smartslurm` command
# (activate / deactivate / uninstall). You do NOT edit any file or add any
# `export` yourself -- avoiding that manual PATH export is the whole point.
# Nothing is compiled or copied into system directories.
#
# Install flow:
#     git clone https://github.com/ld32/SmartSlurm.git
#     cd SmartSlurm
#     bash smartSlurmSetup.sh
#
# Then open a new shell (or `source ~/.bashrc`) and run `smartslurm activate`.
# Activation stays explicit, exactly like `conda activate`.
#
# Optional: pass a different rc file as the first argument:
#     bash smartSlurmSetup.sh ~/.bash_profile
#
# This is the bash version. (zsh support is planned.)
set -u

# Resolve the SmartSlurm checkout from this script's own location (repo root).
home="$(cd "$(dirname "$0")" && pwd)"
hook="$home/bin/smartslurm.sh"

rc="${1:-$HOME/.bashrc}"

if [ ! -f "$hook" ]; then
    echo "smartSlurmSetup: hook not found: $hook" >&2
    echo "  (is this script still inside a SmartSlurm checkout?)" >&2
    exit 1
fi

begin='# >>> smartSlurm initialize >>>'
end='# <<< smartSlurm initialize <<<'

# The block is intentionally tiny: it only points at SMARTSLURM_HOME and sources
# the hook. All logic lives in bin/smartslurm.sh, so a `git pull` upgrades
# behavior without ever rewriting the rc file.
block="$begin
# Managed by smartSlurmSetup -- remove with 'smartslurm uninstall'.
# Enables: smartslurm activate | deactivate | uninstall | status
export SMARTSLURM_HOME=\"$home\"
[ -f \"\$SMARTSLURM_HOME/bin/smartslurm.sh\" ] && . \"\$SMARTSLURM_HOME/bin/smartslurm.sh\"
$end"

touch "$rc"

if grep -qF "$begin" "$rc"; then
    # Idempotent: replace the existing block in place (updates SMARTSLURM_HOME
    # if the checkout moved). Any stray duplicate blocks are collapsed into one.
    bak="$rc.smartslurm.bak.$(date +%s)"
    cp "$rc" "$bak"
    tmp="$(mktemp)"
    awk -v blk="$block" -v b="$begin" -v e="$end" '
        $0==b && !done { print blk; done=1; skip=1; next }  # replace first block
        $0==b          { skip=1; next }                      # swallow duplicates
        $0==e          { skip=0; next }
        skip!=1        { print }
    ' "$rc" > "$tmp"
    cat "$tmp" > "$rc"        # truncate-in-place: keep rc perms/inode
    rm -f "$tmp"
    echo "smartSlurmSetup: updated the smartSlurm block in $rc"
    echo "  backup: $bak"
else
    printf '\n%s\n' "$block" >> "$rc"
    echo "smartSlurmSetup: added the smartSlurm block to $rc"
fi

echo
echo "Next steps:"
echo "  1) start a new shell   (or run:  source \"$rc\")"
echo "  2) smartslurm activate"
