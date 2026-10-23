#!/bin/bash
# smartslurm.sh -- shell hook that defines the `smartslurm` command.
#
# This file is meant to be SOURCED (not executed), normally from the block that
# `smartSlurmSetup` writes into your ~/.bashrc:
#
#     # >>> smartSlurm initialize >>>
#     export SMARTSLURM_HOME="/path/to/SmartSlurm"
#     [ -f "$SMARTSLURM_HOME/bin/smartslurm.sh" ] && . "$SMARTSLURM_HOME/bin/smartslurm.sh"
#     # <<< smartSlurm initialize <<<
#
# It provides, analogous to `conda`:
#     smartslurm activate     - put SmartSlurm on PATH; alias sbatch -> ssbatch
#     smartslurm deactivate   - remove ONLY SmartSlurm's own PATH/prompt additions
#     smartslurm uninstall    - delete the initialize block from the rc file
#     smartslurm status       - show whether SmartSlurm is active
#     smartslurm help
#
# Defining `smartslurm` as a shell function (rather than a script on PATH) is
# what lets activate/deactivate change the *current* shell's PATH and prompt.

# --- locate SMARTSLURM_HOME -------------------------------------------------
# The rc block exports SMARTSLURM_HOME. If the hook is sourced directly without
# it (e.g. in tests), derive it from this file's own location (bin/..).
if [ -z "${SMARTSLURM_HOME:-}" ]; then
    if [ -n "${BASH_SOURCE:-}" ]; then
        SMARTSLURM_HOME="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
        export SMARTSLURM_HOME
    fi
fi

# Prompt token -- curly braces, deliberately distinct from conda's "(env) ".
# (Assigned plainly rather than via ${:=...}: a '}' in the default value would
# prematurely close that parameter expansion.)
if [ -z "${SMARTSLURM_PROMPT_TOKEN:-}" ]; then
    SMARTSLURM_PROMPT_TOKEN='{smartSlurm} '
fi

# --- helpers ----------------------------------------------------------------

# Rebuild PATH keeping every entry EXCEPT the exact dirs we added at activate.
# Surgical by design: anything added to PATH after activation (conda, module
# load, the user's own exports) is preserved untouched.
_ss_path_strip() {
    # args: colon-list of dirs to remove
    local remove="$1" oldifs="$IFS" out="" d e drop
    IFS=:
    for d in $PATH; do
        drop=""
        for e in $remove; do
            [ "$d" = "$e" ] && { drop=1; break; }
        done
        [ -n "$drop" ] && continue
        if [ -z "$out" ]; then out="$d"; else out="$out:$d"; fi
    done
    IFS="$oldifs"
    printf '%s' "$out"
}

_ss_sbatch_path() { command -v sbatch 2>/dev/null; }

# --- subcommands ------------------------------------------------------------

_ss_activate() {
    if [ -n "${_SMARTSLURM_ACTIVE:-}" ]; then
        echo "smartSlurm: already active" >&2
        return 0
    fi
    if [ -z "${SMARTSLURM_HOME:-}" ] || [ ! -d "$SMARTSLURM_HOME/bin" ]; then
        echo "smartSlurm: SMARTSLURM_HOME is not set to a valid SmartSlurm checkout" >&2
        return 1
    fi

    local bin="$SMARTSLURM_HOME/bin" sb="$SMARTSLURM_HOME/sbatchBin"
    local entries="$bin"
    [ -d "$sb" ] && entries="$bin:$sb"

    export PATH="$entries:$PATH"
    export _SMARTSLURM_PATH_ENTRIES="$entries"

    # Prompt: prepend our token once (idempotent).
    export _SMARTSLURM_PROMPT="$SMARTSLURM_PROMPT_TOKEN"
    case "$PS1" in
        *"$_SMARTSLURM_PROMPT"*) ;;                       # already present
        *) PS1="$_SMARTSLURM_PROMPT$PS1" ;;
    esac

    export _SMARTSLURM_ACTIVE=1

    # Announce -- especially the sbatch aliasing, which was previously silent.
    echo "smartSlurm activated."
    local sbp; sbp="$(_ss_sbatch_path)"
    case "$sbp" in
        "$sb/"*)
            echo "  sbatch -> $sbp"
            echo "  (plain 'sbatch' now runs smart ssbatch; 'smartslurm deactivate' restores the system sbatch)"
            ;;
        *)
            echo "  ssbatch and tools are on PATH."
            [ -d "$sb" ] || echo "  (no sbatchBin/ in this checkout; plain 'sbatch' is unchanged)"
            ;;
    esac
}

_ss_deactivate() {
    if [ -z "${_SMARTSLURM_ACTIVE:-}" ]; then
        echo "smartSlurm: not active" >&2
        return 0
    fi

    # Surgical PATH removal: drop only our recorded dirs.
    export PATH="$(_ss_path_strip "${_SMARTSLURM_PATH_ENTRIES:-}")"

    # Prompt: remove only our exact token; leave any other prefix (e.g. conda's
    # "(env) ") in place.
    if [ -n "${_SMARTSLURM_PROMPT:-}" ]; then
        PS1="${PS1/$_SMARTSLURM_PROMPT/}"
    fi

    unset _SMARTSLURM_ACTIVE _SMARTSLURM_PATH_ENTRIES _SMARTSLURM_PROMPT
    echo "smartSlurm deactivated. System sbatch restored."
}

_ss_uninstall() {
    # Clean the current shell first.
    [ -n "${_SMARTSLURM_ACTIVE:-}" ] && _ss_deactivate

    local rc="${1:-$HOME/.bashrc}"
    local begin='# >>> smartSlurm initialize >>>'
    local end='# <<< smartSlurm initialize <<<'

    if [ ! -f "$rc" ]; then
        echo "smartSlurm: no rc file at $rc" >&2
        return 1
    fi
    if ! grep -qF "$begin" "$rc"; then
        echo "smartSlurm: no initialize block found in $rc (nothing to remove)" >&2
        return 0
    fi

    local bak="$rc.smartslurm.bak.$(date +%s)"
    cp "$rc" "$bak"

    local tmp; tmp="$(mktemp)"
    awk -v b="$begin" -v e="$end" '
        $0==b { skip=1; next }
        $0==e { skip=0; next }
        skip!=1 { print }
    ' "$rc" > "$tmp"
    cat "$tmp" > "$rc"        # truncate-in-place: keep rc perms/inode
    rm -f "$tmp"

    echo "smartSlurm: removed initialize block from $rc (backup: $bak)"
    echo "  New shells will no longer load smartSlurm."
}

_ss_status() {
    if [ -n "${_SMARTSLURM_ACTIVE:-}" ]; then
        echo "smartSlurm: ACTIVE"
        echo "  home:    ${SMARTSLURM_HOME:-<unset>}"
        echo "  on PATH: ${_SMARTSLURM_PATH_ENTRIES:-<none>}"
        local sbp; sbp="$(_ss_sbatch_path)"
        echo "  sbatch:  ${sbp:-<none>}"
        case "$sbp" in
            "$SMARTSLURM_HOME/sbatchBin/"*) echo "           (aliased to smart ssbatch)" ;;
            *)                              echo "           (system sbatch)" ;;
        esac
    else
        echo "smartSlurm: inactive"
        echo "  home:    ${SMARTSLURM_HOME:-<unset>}"
    fi
}

_ss_help() {
    cat <<'EOF'
smartslurm -- manage SmartSlurm activation in your shell

Usage:
  smartslurm activate      Add SmartSlurm to PATH and alias sbatch -> ssbatch.
                           Prefixes the prompt with "{smartSlurm} ".
  smartslurm deactivate    Remove only SmartSlurm's own PATH and prompt
                           additions (anything you added afterward is kept).
  smartslurm uninstall     Remove the smartSlurm block from your shell rc file.
  smartslurm status        Show whether SmartSlurm is currently active.
  smartslurm help          Show this message.
EOF
}

# --- dispatcher -------------------------------------------------------------
smartslurm() {
    local cmd="${1:-help}"
    [ "$#" -gt 0 ] && shift
    case "$cmd" in
        activate)   _ss_activate "$@" ;;
        deactivate) _ss_deactivate "$@" ;;
        uninstall)  _ss_uninstall "$@" ;;
        status)     _ss_status "$@" ;;
        help|-h|--help) _ss_help ;;
        *) echo "smartslurm: unknown command '$cmd'" >&2; _ss_help; return 2 ;;
    esac
}
