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

# jrSchemaVersion <file> -> prints 2, 1, or 0 (unknown/missing). Classifies by
# scanning EVERY data row, not just the first line, so ORDERING and HEADERS never
# fool it. This matters because a store can go hybrid in either direction: the old
# writer appended v2 rows under a v1 file, and -- during a mixed-version cluster
# rollout -- an old node can append v1 rows into a store a new node already
# migrated to v2. A line-1-only check called the latter "pure v2", the migration
# was skipped, and the v2 prune then read those v1 rows through the wrong columns
# and deleted COMPLETED history. So:
#   * ANY v1-shaped data row present (numeric jobID, 20 or 21 cols -- the two real
#     widths: 21 w/ trailing epoch, 20 w/ date and no epoch) -> 1 (needs migration),
#     even when v2 rows or a v2 header are also present (a mixed store IS v1-dirty).
#   * otherwise, any v2-shaped data row OR a lone v2 header -> 2 (clean v2).
#   * otherwise (only a v1 header, or nothing recognizable) -> 1 if a v1 header is
#     present, else 0. Blank lines and junk rows do not vote.
jrSchemaVersion() {
    local f=$1
    [ -s "$f" ] || { echo 0; return; }
    # Canonical v2 (a no-op for migration) means the v2 header LEADS the file and
    # no v1 row appears anywhere. Everything else recognizable is "needs
    # normalization" (1): v1/mixed stores, AND a headerless or mis-ordered v2 file
    # -- because the v2 prune treats line 1 as the header (NR==1), so a leading v2
    # DATA row would be eaten. Migration (jrEmitV2) rewrites any of these to the
    # canonical header-led form; a file with nothing recognizable is left alone (0).
    awk -F, '
        /^[[:space:]]*$/            { next }                       # blank line (no vote, no "first")
        !seen { seen=1; firstIsV2H = ($1=="pipeline") }           # is the FIRST real line the v2 header?
        $1=="pipeline"             { h2=1; next }                  # v2 header
        $1=="1jobID"               { h1=1; next }                  # v1 header (18- or 21-col)
        $1 ~ /^[0-9]+$/ && NF>=20  { v1=1; next }                  # v1 data row (20/21 col; >21 if a comma crept into the old date)
        $1 !~ /^[0-9]+$/ && NF==17 { v2=1; next }                  # v2 data row
                                   { }                             # junk: does not vote
        END {
            if (v1)                 { print 1; exit }   # any v1 row -> migrate (even if mixed)
            if (firstIsV2H)         { print 2; exit }   # header-led, no v1 -> canonical v2 (junk tolerated)
            if (v2 || h1 || h2)     { print 1; exit }   # headerless/mis-ordered v2, or a bare v1/v2 header -> normalize
            print 0                                     # nothing recognizable -> leave alone
        }' "$f"
}

# jrCountUnrecognized <file> -> count of non-blank, non-header rows that are
# NEITHER a v1 record (numeric jobID, 20 or 21 cols) NOR a v2 record (non-numeric
# first field, 17 cols). These are the rows jrEmitV2 drops on conversion -- the
# ones that must never vanish silently. Used for the migration warning and lint.
jrCountUnrecognized() {
    [ -s "$1" ] || { echo 0; return; }
    awk -F, '
        /^[[:space:]]*$/            { next }
        $1=="pipeline" || $1=="1jobID" { next }
        $1 ~ /^[0-9]+$/ && (NF==20 || NF==21) { next }
        $1 !~ /^[0-9]+$/ && NF==17 { next }
        { c++ } END { print c+0 }' "$1" 2>/dev/null
}

# jrLint <file> -- print a human health report of a store to stdout (schema
# version + a row-shape census + the line numbers of unrecognized rows). Pure
# read; changes nothing. The audit tool behind `jrRepair.sh --check` so a stale
# or malformed store on a real cluster is a visible finding, never a silent one.
# Returns 0 if the store is clean (v1 or v2, no unrecognized rows), 1 otherwise.
jrLint() {
    local f=$1 ver
    if [ ! -f "$f" ]; then echo "  $f: (absent)"; return 0; fi
    if [ ! -s "$f" ]; then echo "  $f: (empty)"; return 0; fi
    ver=$(jrSchemaVersion "$f")
    awk -F, -v file="$f" -v ver="$ver" '
        /^[[:space:]]*$/            { blank++; next }
        $1=="pipeline"             { h2++; next }
        $1=="1jobID"               { h1++; next }
        $1 ~ /^[0-9]+$/ && NF>=21  { v1_21++; next }
        $1 ~ /^[0-9]+$/ && NF==20  { v1_20++; next }
        $1 !~ /^[0-9]+$/ && NF==17 { v2++; next }
        { bad++; if (nbad<5) badln[++nbad]=NR" (NF="NF"): "substr($0,1,60) }
        END {
            tag = (ver==2?"v2 (current)":(ver==1?"v1/mixed -> will auto-migrate":"UNRECOGNIZED"))
            printf "  %s\n", file
            printf "    schema version : %s  [%s]\n", ver, tag
            printf "    v2 rows        : %d\n", v2+0
            if (v1_21||v1_20) printf "    v1 rows        : %d  (21-col:%d  20-col/no-epoch:%d)\n", v1_21+v1_20, v1_21+0, v1_20+0
            if (h2||h1)       printf "    header lines   : v2:%d v1:%d\n", h2+0, h1+0
            if (blank)        printf "    blank lines    : %d\n", blank
            if (bad) {
                printf "    UNRECOGNIZED   : %d row(s) -- will NOT migrate (kept in backup); first few:\n", bad
                for (i=1;i<=nbad;i++) printf "        line %s\n", badln[i]
            }
            exit (bad>0 ? 1 : 0)
        }' "$f"
}

# _jrWarnUnknown <file> -- warn (once per file state) that a store matches no known
# schema. Deduped via a marker so a per-job hot path (ssbatch/cleanUp) does not
# flood: we warn again only after the file changes (marker older than the file).
_jrWarnUnknown() {
    local f=$1 m="$1.schemaWarning"
    [ -f "$m" ] && [ ! "$f" -nt "$m" ] && return 0
    echo "SmartSlurm: WARNING: jobRecord store matches no known schema (not v1 or v2); left untouched: $f" >&2
    echo "SmartSlurm:          inspect it with: jrRepair.sh --check \"$(dirname -- "$f")\"" >&2
    : > "$m" 2>/dev/null
}

# ===========================================================================
# ONE lock, ONE I/O API. Every mutation of a store (migrate, append, prune/rewrite)
# goes through these so no two writers can interleave and no tool can hand-roll the
# format. The lock is a directory (mkdir is atomic even on NFS, unlike flock), held
# only for the critical section. A crashed holder cannot wedge the cluster: waiters
# reclaim a lock older than JR_LOCK_STALE seconds, and give up waiting after a bound.
# ---------------------------------------------------------------------------
JR_LOCK_STALE="${JR_LOCK_STALE:-30}"   # seconds before a held lock is treated as stale
JR_LOCK_TRIES="${JR_LOCK_TRIES:-200}"  # poll attempts (x0.05s) before giving up waiting

# _jrLock <file> -> 0 if the lock was acquired (caller MUST _jrUnlock), 1 on timeout.
_jrLock() {
    local lock="$1.lock" i=0 age now
    while :; do
        mkdir "$lock" 2>/dev/null && return 0
        # reclaim a stale lock (holder likely died mid-write)
        now=$(date +%s); age=$(( now - $(stat -c %Y "$lock" 2>/dev/null || echo "$now") ))
        [ "$age" -ge "$JR_LOCK_STALE" ] && rmdir "$lock" 2>/dev/null && continue
        i=$((i+1)); [ "$i" -ge "$JR_LOCK_TRIES" ] && return 1
        sleep 0.05
    done
}
_jrUnlock() { rmdir "$1.lock" 2>/dev/null; }

# _jrMigrateUnlocked <file> -- the migration body WITHOUT taking the lock (the
# caller already holds it). Converts a v1/mixed/headerless-v2 store to canonical
# v2 in place: back up, convert via jrEmitV2, atomic rename. NOT SILENT about
# dropped (unrecognized) rows or an unrecognized file.
_jrMigrateUnlocked() {
    local f=$1 ver bak tmp bad
    ver=$(jrSchemaVersion "$f")
    [ "$ver" = 2 ] && return 0
    if [ "$ver" = 0 ]; then _jrWarnUnknown "$f"; return 0; fi
    bak="$f.back.$(date +%Y-%m-%d_%H-%M-%S).txt"
    cp -p "$f" "$bak" 2>/dev/null
    bad=$(jrCountUnrecognized "$f")
    tmp="$f.v2.$$"
    if jrEmitV2 "$f" > "$tmp" && [ -s "$tmp" ]; then
        mv "$tmp" "$f"
        if [ "${bad:-0}" -gt 0 ]; then
            echo "SmartSlurm: WARNING: migrated $f to v2, but $bad row(s) matched no known schema and were NOT carried into the live store (kept in backup: $bak). Inspect with: jrRepair.sh --check \"$(dirname -- "$f")\"" >&2
        else
            echo "SmartSlurm: migrated legacy jobRecord to v2 schema: $f (backup: $bak)" >&2
        fi
    else
        rm -f "$tmp"
    fi
}

# jrMigrateStore <file> -- public, locked, idempotent migration of a legacy (v1)/
# mixed store to pure v2, so the v2 writer/prune never operate on a mixed file (the
# cause of the data-loss bug: the v1-positional prune deleted appended v2 rows).
# The backup (<file>.back.<ts>.txt) is the rollback, same as jrRepair.
jrMigrateStore() {
    local f=$1 ver
    [ -s "$f" ] || return 0
    ver=$(jrSchemaVersion "$f")
    [ "$ver" = 2 ] && return 0                 # already clean v2 -> no-op (no lock needed)
    if [ "$ver" = 0 ]; then _jrWarnUnknown "$f"; return 0; fi
    if _jrLock "$f"; then
        _jrMigrateUnlocked "$f"                # re-checks version under the lock
        _jrUnlock "$f"
    fi
    # lock timeout -> another job is migrating; it will become v2 shortly. A reader
    # (jrEmitV2) is hybrid-safe meanwhile, so skipping here loses nothing.
}

# jrEnsureHeader <file> -- make sure the store exists AND is v2 before writing:
# create it with the v2 header if absent, else migrate a v1/mixed store in place.
jrEnsureHeader() {
    if [ ! -f "$1" ]; then printf '%s\n' "$JR_HEADER_V2" > "$1"; return 0; fi
    jrMigrateStore "$1"
}

# jrAppendLine <file> <preformatted-v2-row> -- the ONE sanctioned way to append a
# row. Under the lock it: creates the v2 header if the file is new, migrates a
# v1/mixed store to v2 first (so we never append a v2 row onto a v1 file and make a
# hybrid), then appends atomically. If the lock cannot be taken in time it appends
# anyway and warns -- a possibly-interleaved row is far better than a LOST record.
jrAppendLine() {
    local f=$1 row=$2
    if _jrLock "$f"; then
        [ -f "$f" ] || printf '%s\n' "$JR_HEADER_V2" > "$f"
        _jrMigrateUnlocked "$f"
        printf '%s\n' "$row" >> "$f"
        _jrUnlock "$f"
    else
        [ -f "$f" ] || printf '%s\n' "$JR_HEADER_V2" > "$f"
        printf '%s\n' "$row" >> "$f"
        echo "SmartSlurm: WARNING: jobRecord append proceeded WITHOUT the lock (contention/stale): $f" >&2
    fi
}

# jrAppendRow <file> name=value ... -- build a sanitized v2 row with jrRow and
# append it via jrAppendLine. The convenience front door for writers.
jrAppendRow() {
    local f=$1; shift
    jrAppendLine "$f" "$(jrRow "$@")"
}

# jrRewrite <file> <awk-row-program> [awk -v args ...] -- the ONE sanctioned way to
# prune/rewrite a store. Under the lock it migrates to v2, then runs the row program
# over every DATA row with col["name"] available (OFS=","), keeping rows the program
# `print`s; the header is always preserved. Atomic (temp + rename). The two prunes
# (7-day recency, fail-prune) are one call each -- no positional columns, no second
# colmap, no v1 fallback. Skips a store that is not (or cannot become) v2.
jrRewrite() {
    local f=$1 prog=$2; shift 2
    [ -s "$f" ] || return 0
    _jrLock "$f" || return 0
    _jrMigrateUnlocked "$f"
    if [ "$(jrSchemaVersion "$f")" = 2 ]; then
        local tmp="$f.rw.$$"
        if awk -F, -v OFS=, "$@" '
                NR==1 { for(i=1;i<=NF;i++) col[$i]=i; print; next }
                { '"$prog"' }
            ' "$f" > "$tmp" && [ -s "$tmp" ]; then
            mv "$tmp" "$f"
        else
            rm -f "$tmp"
        fi
    fi
    _jrUnlock "$f"
}

# jrRow name=value ... -> one CSV row with the fields placed in JR_COLUMNS_V2 order.
# Unspecified columns are emitted empty. This is the single writer: callers name
# fields, so the on-disk column order lives only in JR_COLUMNS_V2 and a row can
# never silently drift from the header. Portable (no associative arrays).
#
# SANITIZED: every value has commas, spaces, tabs and CRs/newlines replaced with
# '_', so a row is ALWAYS exactly 17 clean comma-separated fields. This is what
# keeps field parsing accurate -- the legacy store's human `date` column (from the
# `date` command) carried spaces, and in some locales a comma, which is how a row
# could gain an extra field and shift the rest. A v2 row can never do that. Keys
# (pipeline/step/ref/sample) are tokens and should not contain spaces anyway; if
# one does, it is normalized identically on every write so it stays self-consistent.
jrRow() {
    local c kv val out=""
    for c in $JR_COLUMNS_V2; do
        val=""
        for kv in "$@"; do
            case "$kv" in "$c="*) val="${kv#*=}"; break ;; esac
        done
        val=${val//,/_}; val=${val// /_}; val=${val//$'\t'/_}
        val=${val//$'\n'/_}; val=${val//$'\r'/_}
        out="$out,$val"
    done
    printf '%s\n' "${out#,}"
}

# jrEmitV2 <file> -- write the store to stdout in v2 form (header + rows),
# converting legacy v1 rows on the fly. This is the compat shim for READ-ONLY
# readers and the single definition of the v1->v2 field mapping (jrRepair.sh and
# the auto-migration reuse it).
#
# It classifies PER ROW, not per file, so it is correct for every real variant:
# a pure v1 file (with a 21- or 18-column header, or headerless), a pure v2 file,
# and -- crucially -- a MIXED file (v1 rows with v2 rows appended, the state the
# old writer produced on an un-migrated store). Rules, in order:
#   * a v1 header ("1jobID,...") or v2 header ("pipeline,...") line -> skipped
#   * an already-v2 data row (17 fields) -> passed through unchanged
#   * a v1 data row (numeric jobID in field 1, 20+ fields) -> mapped to v2
#     (inputSize is carried verbatim, even non-numeric like "missingInputFile")
#   * anything else (junk / unexpected shape) -> skipped (jrRepair diverts these
#     to a .skipped file so they are visible, not silently lost)
#
# ROBUST field extraction: the legacy human `date` column (from `date`) carries
# spaces and, in some locales, a COMMA -- which would shift the trailing columns
# under -F,. So epoch is taken from the LAST field only if it looks like a unix
# timestamp (>=9 digits), never from a fixed position; the mapped fields we need
# (jobID, inputSize, mem/time, status, user, program, ref, cores, flag) all live
# in columns 1..15, BEFORE the date, so they are position-stable regardless. The
# carried text fields are sanitized (commas/spaces -> '_') so the emitted v2 row
# is always exactly 17 clean fields. isoDate is carried (sanitized) only when the
# row has a clean width (<=21); on a comma-shifted row it is left blank -- epoch,
# the authoritative timestamp, is preserved either way.
jrEmitV2() {
    local f=$1
    [ -s "$f" ] || return 0
    awk -F, -v OFS=, -v hdr="$JR_HEADER_V2" '
      function sampleOf(fl,  n,a,i){ n=split(fl,a,"."); if(n<2)return ""; i=2;
          while(i<=n && a[i]~/^[0-9]+$/) i++; i++; return (i<=n)?a[i]:"" }
      function san(x){ gsub(/[, \t\r]/,"_",x); return x }
      BEGIN{ print hdr }
      $1=="pipeline" || $1=="1jobID" { next }                 # header line (either schema)
      $1 !~ /^[0-9]+$/ && NF==17 { print; next }              # already-v2 data row -> pass through
      $1 ~ /^[0-9]+$/ && NF>=20 {                             # v1 data row -> map
          ep = ($NF ~ /^[0-9]{9,}$/) ? $NF : ""               # epoch = numeric tail, else none (20-col)
          iso = (NF<=21) ? san($20) : ""                       # human date only on a clean-width row
          print "legacy",san($12),san($13),san(sampleOf($14)),$1,san($2),$7,$5,$8,$6,$15,san($9),san($10),"",ep,iso,san($14); next }
      # else: junk / unexpected -> skip
    ' "$f"
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
