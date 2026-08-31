#!/bin/sh

# SmartSlurm AI helper library.
#
# Source this file to get optional AI support. Every function here is safe to
# call when AI is unavailable: it returns non-zero and prints nothing to stdout,
# so callers can keep working exactly as before.
#
#   . $(dirname $0)/aiHelper.sh
#   aiEnabled && echo "$prompt" | aiAsk > answer.txt
#
# The AI token never appears here: it lives inside the setuid ai_call helper.

# Path to the setuid AI helper binary. Config or environment can override it.
if [ -z "$smartSlurmAiCall" ]; then
    for _aiTry in "$HOME/.smartSlurm/bin/ai_call" \
                  "$(dirname "$0")/ai_call" \
                  "/home/ld32/.smartSlurm/bin/ai_call"; do
        [ -x "$_aiTry" ] && { smartSlurmAiCall="$_aiTry"; break; }
    done
    unset _aiTry
fi
_AI_CALL="$smartSlurmAiCall"

# Seconds to wait for an answer before giving up, so AI never hangs a job.
[ -z "$smartSlurmAiTimeout" ] && smartSlurmAiTimeout=90

# Largest prompt ai_call accepts is 64K; stay well under it.
[ -z "$smartSlurmAiMaxPrompt" ] && smartSlurmAiMaxPrompt=48000

# aiEnabled: true only when AI is turned on and the helper is really runnable.
aiEnabled() {
    [ "$smartSlurmAiEnabled" = "no" ] && return 1
    [ -n "$_AI_CALL" ] && [ -x "$_AI_CALL" ]
}

# aiWhyNot: one line explaining why AI is unavailable, for logs and menus.
aiWhyNot() {
    if [ "$smartSlurmAiEnabled" = "no" ]; then
        echo "AI is disabled by config (smartSlurmAiEnabled=no)."
    elif [ -z "$_AI_CALL" ] || [ ! -x "$_AI_CALL" ]; then
        echo "AI helper ai_call not found or not executable. Ask your SmartSlurm admin to install it and add '$USER' to the owner's ~/.smartSlurm/allowed_users.txt"
    else
        echo "AI call failed. If you see 'Access denied', ask your SmartSlurm admin to add '$USER' to the owner's ~/.smartSlurm/allowed_users.txt"
    fi
}

# aiAsk: read a prompt on stdin, print the answer on stdout.
# Returns non-zero and prints nothing on stdout if AI is unavailable or fails.
# stderr from the helper is kept out of the answer so it cannot pollute emails.
aiAsk() {
    aiEnabled || return 1

    _aiIn=`mktemp 2>/dev/null` || return 1
    _aiOut=`mktemp 2>/dev/null` || { rm -f "$_aiIn"; return 1; }
    _aiErr=`mktemp 2>/dev/null` || { rm -f "$_aiIn" "$_aiOut"; return 1; }

    head -c "$smartSlurmAiMaxPrompt" > "$_aiIn"

    if [ ! -s "$_aiIn" ]; then
        rm -f "$_aiIn" "$_aiOut" "$_aiErr"
        return 1
    fi

    if command -v timeout >/dev/null 2>&1; then
        timeout "$smartSlurmAiTimeout" "$_AI_CALL" < "$_aiIn" > "$_aiOut" 2> "$_aiErr"
    else
        "$_AI_CALL" < "$_aiIn" > "$_aiOut" 2> "$_aiErr"
    fi
    _aiRc=$?

    if [ "$_aiRc" -eq 0 ] && [ -s "$_aiOut" ]; then
        cat "$_aiOut"
    else
        # Log the reason where it is useful, but never on stdout.
        echo "AI call skipped (exit $_aiRc): `head -n 3 \"$_aiErr\" | tr '\n' ' '`" >&2
        _aiRc=1
    fi

    rm -f "$_aiIn" "$_aiOut" "$_aiErr"
    return $_aiRc
}
