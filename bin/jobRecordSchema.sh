#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# jobRecordSchema.sh -- single source of truth for the SmartSlurm jobRecord
# schema. Source it (it defines variables and functions, no side effects):
#
#     . "$(dirname "$0")/jobRecordSchema.sh"
#
# Writers (ssbatch/ssbatchSimple header, cleanUp row), readers
# (estimateResource, calculateMemTime, cleanUp prune) and the repair tool
# (jrRepair.sh) all take the column list and the name->index map from here, so
# the layout is defined in exactly one place.
#
# Versioning is by HEADER SIGNATURE, not a marker line: a v2 file's header
# begins with "pipeline,"; a legacy v1 file's header began with "1jobID,".
# ---------------------------------------------------------------------------

# --- v2 schema (current) ---------------------------------------------------
# This order IS the on-disk column order. Keep names descriptive; the header
# line is literally this list, comma-joined.
JR_COLUMNS_V2="pipeline step ref sample jobID inputSize memUsed memAllocated timeUsed timeAllocated cores status user node epoch isoDate flag"
JR_HEADER_V2="$(printf '%s' "$JR_COLUMNS_V2" | tr ' ' ',')"

# Required subset (every consumer depends on these; the rest are
# readability/future and no reader may assume their presence/position).
JR_REQUIRED_V2="pipeline step ref jobID inputSize memUsed timeUsed memAllocated timeAllocated status epoch"

# --- v1 (legacy) layout, for detection + migration only --------------------
# The REAL 21-field order old cleanUp.sh wrote (the old header line was stale
# and only named 18 of them). jrRepair reads old rows positionally via this.
JR_COLUMNS_V1="jobID inputSize memDefault timeDefault memAllocated timeAllocated memUsed timeUsed status user saccMem program ref flag cores extraMem extraTime ratioM ratioT date epoch"

# awk prelude: build col["name"]=index from the header row, then skip it.
# Usage:  awk -F, "$JR_AWK_COLMAP"'{ print $(col["memUsed"]) }' jobRecord.txt
JR_AWK_COLMAP='NR==1{for(i=1;i<=NF;i++)col[$i]=i; next}'

# jrSchemaVersion <file> -> prints 2, 1, or 0 (unknown/missing).
jrSchemaVersion() {
    local f=$1 h
    [ -s "$f" ] || { echo 0; return; }
    IFS= read -r h < "$f"
    case "$h" in
        pipeline,*)          echo 2 ;;
        1jobID,*|jobID,*)    echo 1 ;;
        *)                   echo 0 ;;
    esac
}

# jrEnsureHeader <file> -- create the file with the v2 header if it is absent.
jrEnsureHeader() { [ -f "$1" ] || printf '%s\n' "$JR_HEADER_V2" > "$1"; }

# jrSampleFromFlag <flag> -> the per-sample identifier, or "" (experiment-level).
# Flag grammar (rAP): step.dep[.dep...].name[.sample[.fanout...]] where step and
# the deps are numeric. So: skip field 1 (step) and any numeric deps, the next
# field is the step name, and the one after that is the sample. (jrRepair.sh does
# the same in awk for stream conversion; keep the two in sync.)
jrSampleFromFlag() {
    local flag=$1 i n s
    local IFS=.
    set -- $flag
    n=$#
    [ "$n" -ge 2 ] || { printf ''; return; }
    i=2
    while [ "$i" -le "$n" ]; do
        eval "s=\${$i}"
        case "$s" in (*[!0-9]*|'') break;; esac   # stop at first non-numeric (the name)
        i=$((i+1))
    done
    i=$((i+1))                                     # past the name -> sample
    if [ "$i" -le "$n" ]; then eval "printf '%s' \"\${$i}\""; else printf ''; fi
}

# jrcol <name> <headerLine> -> prints the 1-based column index of <name> in a
# comma-separated header line, or empty if absent. For bash-side callers.
jrcol() {
    local name=$1 header=$2 i=1 c
    local IFS=,
    for c in $header; do
        [ "$c" = "$name" ] && { printf '%s' "$i"; return 0; }
        i=$((i+1))
    done
    return 1
}
