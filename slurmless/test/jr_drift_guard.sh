#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# jr_drift_guard.sh -- fail if any SmartSlurm tool hand-rolls the jobRecord
# format instead of going through the jobRecordSchema.sh API. This is what makes
# the 18-col-header-vs-21-col-data drift (and the raw-append hybrid bug)
# structurally impossible to reintroduce: the tree does not go green if a tool
# writes a header literal, redirects to a store, or reads it by raw column number.
#
#   SMARTSLURM_ROOT=/path/to/smartslurm test/jr_drift_guard.sh
#
# Exit 0 = clean; non-zero = at least one violation (printed with file:line).
# The store-management layer itself is exempt: jobRecordSchema.sh (the API) and
# jrRepair.sh (the sanctioned manual low-level repair tool).
# ---------------------------------------------------------------------------
set -u
ROOT="${SMARTSLURM_ROOT:-}"
[ -n "$ROOT" ] && [ -d "$ROOT/bin" ] || { echo "jr_drift_guard: set SMARTSLURM_ROOT to a SmartSlurm checkout" >&2; exit 2; }
BIN="$ROOT/bin"
EXEMPT='jobRecordSchema.sh jrRepair.sh'

viol=0
# scan <awk-regex> <message> [file-glob] -- print non-comment offending lines.
scan() {
    local re=$1 msg=$2 only=${3:-}
    local f base out
    for f in "$BIN"/*; do
        [ -f "$f" ] || continue
        base=$(basename "$f")
        case " $EXEMPT " in *" $base "*) continue ;; esac
        [ -n "$only" ] && [ "$base" != "$only" ] && continue
        out=$(awk -v re="$re" '!/^[[:space:]]*#/ && $0 ~ re { printf "    %s:%d: %s\n", FILENAME, FNR, $0 }' "$f")
        if [ -n "$out" ]; then
            [ "$viol" -eq 0 ] && echo "jr_drift_guard: VIOLATIONS"
            echo "  [$msg]"; echo "$out"; viol=1
        fi
    done
}

# 1) No tool but the API may REDIRECT to a jobRecord store (append, truncate, or a
#    .new temp). Writers must use jrAppendLine/jrAppendRow; rewrites use jrRewrite.
scan '>[[:space:]]*"?\$(smartSlurmJobRecordDir|\{smartSlurmJobRecordDir\})/jobRecord' \
     'raw redirect to a jobRecord store -- use jrAppendLine/jrAppendRow/jrRewrite'

# 2) No tool but the API may write a HEADER literal. Use JR_HEADER_V2/jrEnsureHeader.
scan '1jobID,2inputSize|pipeline,step,ref,sample,jobID' \
     'jobRecord header written as a literal -- use jrEnsureHeader (JR_HEADER_V2)'

# 3) No positional read of the store date/epoch columns ($20/$21) in cleanUp -- the
#    exact drift that let the prune mis-read a hybrid row. Readers address by name
#    via jrEmitV2 + JR_AWK_COLMAP (col["epoch"], col["status"], ...).
scan '\$2[01]([^0-9]|$)' \
     'positional $20/$21 (date/epoch) read -- address by name via col["..."]' cleanUp.sh

if [ "$viol" -eq 0 ]; then
    echo "jr_drift_guard: clean -- all store I/O goes through the jobRecordSchema.sh API."
    exit 0
fi
exit 1
