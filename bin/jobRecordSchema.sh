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

# jrSchemaVersion <file> -> prints 2, 1, or 0 (unknown/missing). Handles files
# with a header AND legacy headerless v1 files (the old 7-day prune dropped the
# stale header, so real v1 stores are often headerless): infer from the first
# row's shape -- numeric first field + 21 columns = v1; non-numeric first field +
# 17 columns = v2.
jrSchemaVersion() {
    local f=$1 h f1 nf
    [ -s "$f" ] || { echo 0; return; }
    IFS= read -r h < "$f"
    case "$h" in
        pipeline,*)   echo 2; return ;;
        1jobID,*)     echo 1; return ;;
    esac
    nf=$(printf '%s\n' "$h" | awk -F, '{print NF; exit}')
    f1=${h%%,*}
    case "$f1" in
        ''|*[!0-9]*)  [ "${nf:-0}" -eq 17 ] && echo 2 || echo 0 ;;   # non-numeric id -> v2-shaped
        *)            [ "${nf:-0}" -eq 21 ] && echo 1 || echo 0 ;;   # numeric jobID -> v1-shaped
    esac
}

# jrEnsureHeader <file> -- create the file with the v2 header if it is absent.
jrEnsureHeader() { [ -f "$1" ] || printf '%s\n' "$JR_HEADER_V2" > "$1"; }

# jrRow name=value ... -> one CSV row with the fields placed in JR_COLUMNS_V2 order.
# Unspecified columns are emitted empty. This is the single writer: callers name
# fields, so the on-disk column order lives only in JR_COLUMNS_V2 and a row can
# never silently drift from the header. Values must not contain commas/newlines
# (CSV). Portable (no associative arrays): 17 columns x a few args is trivial.
jrRow() {
    local c kv val out=""
    for c in $JR_COLUMNS_V2; do
        val=""
        for kv in "$@"; do
            case "$kv" in "$c="*) val="${kv#*=}"; break ;; esac
        done
        out="$out,$val"
    done
    printf '%s\n' "${out#,}"
}

# jrEmitV2 <file> -- write the store to stdout in v2 form (header + rows),
# converting a legacy v1 file on the fly. This is the compat shim for READ-ONLY
# readers: they pipe through it and parse v2 by name, working on either version
# without touching the file. It is also the single definition of the v1->v2 field
# mapping (jrRepair.sh reuses it). A v2 file passes through; a v1 file is mapped
# (malformed, non-21-field rows are dropped here -- they are unfittable; jrRepair
# diverts them to a .skipped file instead). Unknown header -> nothing.
jrEmitV2() {
    local f=$1
    [ -s "$f" ] || return 0
    case "$(jrSchemaVersion "$f")" in
        2) cat "$f" ;;
        1) awk -F, -v OFS=, -v hdr="$JR_HEADER_V2" '
             function sampleOf(fl,  n,a,i){ n=split(fl,a,"."); if(n<2)return ""; i=2;
                 while(i<=n && a[i]~/^[0-9]+$/) i++; i++; return (i<=n)?a[i]:"" }
             BEGIN{ print hdr }
             # only real 21-field v1 data rows map; the 18-field stale header (if
             # present) and any malformed row are naturally excluded. Works for
             # header-present and legacy headerless v1 files alike.
             NF==21 { print "legacy",$12,$13,sampleOf($14),$1,$2,$7,$5,$8,$6,$15,$9,$10,"",$21,$20,$14 }
           ' "$f" ;;
    esac
}

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

# --- optional ref size-class pooling ---------------------------------------
# Some refs are close enough in size that a tool's resource use is comparable;
# pooling their history gives better fits. JR_REF_CLASS maps "ref=class" pairs
# (space-separated). Default pools the common mammalian genomes (mouse + human,
# all ~3 Gbp). Override by setting JR_REF_CLASS in the environment/config, or ""
# to disable (each ref is its own group). "none" (ref-less tools) never pools.
JR_REF_CLASS="${JR_REF_CLASS-mm9=mammal mm10=mammal mm39=mammal hg19=mammal hg38=mammal GRCm38=mammal GRCm39=mammal GRCh37=mammal GRCh38=mammal}"

# jrRefClass <ref> -> the group name for <ref> (the mapped class, or <ref> itself).
jrRefClass() {
    local r=$1 pair
    [ "$r" = none ] && { printf 'none'; return; }
    for pair in $JR_REF_CLASS; do
        case "$pair" in "$r"=*) printf '%s' "${pair#*=}"; return ;; esac
    done
    printf '%s' "$r"
}

# jrRefMatchRe <ref> -> an ERE matching every ref in the same group (so an extract
# can pull the whole class). Unmapped ref -> exact match on itself.
jrRefMatchRe() {
    local r=$1 cls pair members=""
    cls=$(jrRefClass "$r")
    [ "$cls" = "$r" ] && { printf '^%s$' "$r"; return; }
    for pair in $JR_REF_CLASS; do
        case "$pair" in *="$cls") members="$members|${pair%=*}" ;; esac
    done
    printf '^(%s)$' "${members#|}"
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
