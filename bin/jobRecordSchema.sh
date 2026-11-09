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
        !seen { seen=1; firstIsV2H = ($1=="pipeline" && $2=="step") }   # is the FIRST real line the v2 header?
        $1=="pipeline" && $2=="step" { h2=1; next }               # v2 header (NOT a row whose pipeline is named "pipeline")
        $1=="1jobID"               { h1=1; next }                  # v1 header (18- or 21-col)
        $1 ~ /^[0-9]+$/ && NF>=20  { v1=1; next }                  # v1 data row (20/21 col; >21 if a comma crept into the old date)
        NF==17 && (h2 || $1 !~ /^[0-9]+$/) { v2=1; next }          # v2 data row (17 col). An all-digit $1 is ambiguous (v2 pipeline name vs corrupt v1 jobID), so it votes v2 ONLY under a v2 header already seen above; headerless, it stays unrecognized.
                                   { }                             # junk: does not vote
        END {
            if (v1)                 { print 1; exit }   # any v1 row -> migrate (even if mixed)
            if (firstIsV2H)         { print 2; exit }   # header-led, no v1 -> canonical v2 (junk tolerated)
            if (v2 || h1 || h2)     { print 1; exit }   # headerless/mis-ordered v2, or a bare v1/v2 header -> normalize
            print 0                                     # nothing recognizable -> leave alone
        }' "$f"
}

# A row is RECOGNIZED iff it is a v1 record (numeric jobID, 20+ cols -- 20/21, or
# more if a comma crept into the old date) or a v2 record (17 cols). NF separates
# the clear cases (v1>=20, v2==17). The one ambiguity is a 17-col row with an
# all-digit field 1: that is EITHER a v2 row whose pipeline name is all-digit OR a
# corrupt/truncated v1 row (jobID-first). We disambiguate by HEADER CONTEXT -- such
# a row counts as v2 only when a v2 header has already appeared above it in the store
# (the shape every sanctioned v2 writer produces); headerless, it stays UNRECOGNIZED
# rather than being guessed into v2. A non-digit-$1 17-col row is always v2.
# Everything else (blank lines aside) is UNRECOGNIZED: jrEmitV2 drops it,
# so it must never vanish silently. These three helpers share ONE definition of
# "recognized", matching jrEmitV2 exactly, so the migration integrity check and the
# warnings can't disagree with what conversion actually did.

# jrUnrecognizedLines <file> -> emit the unrecognized rows verbatim (for a .skipped
# / .unrecognized divert file and for jrRepair). Blank lines and headers excluded.
jrUnrecognizedLines() {
    [ -s "$1" ] || return 0
    awk -F, '
        /^[[:space:]]*$/            { next }
        $1=="pipeline" && $2=="step" { v2h=1; next }
        $1=="1jobID"                { next }
        $1 ~ /^[0-9]+$/ && NF>=20   { next }
        NF==17 && (v2h || $1 !~ /^[0-9]+$/) { next }
        { print }' "$1" 2>/dev/null
}
# jrCountUnrecognized <file> -> how many such rows.
jrCountUnrecognized() { jrUnrecognizedLines "$1" | grep -c . ; }
# jrCountRecognized <file> -> how many v1+v2 data rows (the ones jrEmitV2 keeps).
jrCountRecognized() {
    [ -s "$1" ] || { echo 0; return; }
    awk -F, '
        $1=="pipeline" && $2=="step" { v2h=1; next }
        $1=="1jobID"                { next }
        ($1 ~ /^[0-9]+$/ && NF>=20) || (NF==17 && (v2h || $1 !~ /^[0-9]+$/)) { c++ }
        END { print c+0 }' "$1" 2>/dev/null
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
        $1=="pipeline" && $2=="step" { h2++; next }
        $1=="1jobID"               { h1++; next }
        $1 ~ /^[0-9]+$/ && NF>=21  { v1_21++; next }
        $1 ~ /^[0-9]+$/ && NF==20  { v1_20++; next }
        NF==17 && (h2 || $1 !~ /^[0-9]+$/) { v2++; next }
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
            # A valid store has exactly ONE header. More than one header-shaped line is
            # either a concatenated duplicate header OR a data row whose pipeline+step
            # literally equal the header ("pipeline"+"step", or "1jobID") -- which is
            # byte-identical to the header and cannot be told apart from it, so it is
            # read as a header and its data silently excluded. Make that loud.
            if (h2>1 || h1>1) {
                printf "    DUPLICATE HEADER: %d extra header-shaped line(s) -- a store has ONE header. An extra one is a stray duplicate header, OR a data row whose pipeline/step are literally the header words (rename that pipeline: it is indistinguishable from the header and its row is not counted as data).\n", (h2>1?h2-1:0)+(h1>1?h1-1:0)
            }
            exit ((bad>0 || h2>1 || h1>1) ? 1 : 0)
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
# v2 in place: back up, convert via jrEmitV2, VERIFY the result, atomic rename.
#
# Post-conversion integrity gate (never commit a bad conversion): the temp output
# must (a) be canonical v2, and (b) preserve EVERY recognizable (v1+v2) input row
# -- output data-row count == recognized input-row count. Unrecognized junk is
# allowed to drop (that is expected, and it is saved separately); but if the output
# is not clean v2, or a RECOGNIZABLE row went missing, the migration ABORTS: the
# live store is left untouched, the redundant backup is removed, and it errors
# loudly. NOT SILENT about dropped rows (saved to a .unrecognized file + backup) or
# an unrecognized file.
_jrMigrateUnlocked() {
    local f=$1 ver bak tmp bad recog outrows unrec
    ver=$(jrSchemaVersion "$f")
    [ "$ver" = 2 ] && return 0
    if [ "$ver" = 0 ]; then _jrWarnUnknown "$f"; return 0; fi
    bak="$f.back.$(date +%Y-%m-%d_%H-%M-%S).txt"
    cp -p "$f" "$bak" 2>/dev/null
    bad=$(jrCountUnrecognized "$f")
    recog=$(jrCountRecognized "$f")
    tmp="$f.v2.$$"
    if ! jrEmitV2 "$f" > "$tmp" || [ ! -s "$tmp" ]; then
        rm -f "$tmp" "$bak"
        echo "SmartSlurm: ERROR: jobRecord migration produced no output; store left UNTOUCHED: $f" >&2
        return 1
    fi
    # integrity gate ------------------------------------------------------------
    outrows=$(( $(grep -c . "$tmp") - 1 ))              # data rows (minus the header)
    if [ "$(jrSchemaVersion "$tmp")" != 2 ] || [ "${outrows:-0}" -ne "${recog:-0}" ]; then
        rm -f "$tmp" "$bak"
        echo "SmartSlurm: ERROR: jobRecord migration failed its integrity check (output canonical-v2? / kept ${outrows:-0} of ${recog:-0} recognizable rows); store left UNTOUCHED: $f -- inspect with: jrRepair.sh --check \"$(dirname -- "$f")\"" >&2
        return 1
    fi
    # commit --------------------------------------------------------------------
    mv "$tmp" "$f"
    if [ "${bad:-0}" -gt 0 ]; then
        unrec="$f.unrecognized.$(date +%Y-%m-%d_%H-%M-%S).txt"
        jrUnrecognizedLines "$bak" > "$unrec" 2>/dev/null
        echo "SmartSlurm: WARNING: migrated $f to v2; $bad row(s) matched no known schema and were NOT carried over (saved to $unrec; full original in $bak). Inspect with: jrRepair.sh --check \"$(dirname -- "$f")\"" >&2
    else
        echo "SmartSlurm: migrated legacy jobRecord to v2 schema: $f (backup: $bak)" >&2
    fi
    # Case 2 (header collision): a data row whose pipeline/step are literally the header
    # words ("pipeline"/"step", or "1jobID") is byte-identical to the header, so it was
    # read as a (duplicate) header and its data excluded -- not counted in recog or bad,
    # so the gate above can't see it. Surface it loudly and preserve it: divert every
    # header-shaped line BEYOND THE FIRST (the real header) to the .unrecognized file.
    local dup; dup=$(awk -F, '(($1=="pipeline"&&$2=="step")||$1=="1jobID"){ if(seen++) print }' "$bak" 2>/dev/null)
    if [ -n "$dup" ]; then
        [ -n "${unrec:-}" ] || unrec="$f.unrecognized.$(date +%Y-%m-%d_%H-%M-%S).txt"
        printf '%s\n' "$dup" >> "$unrec"
        echo "SmartSlurm: WARNING: $f had $(printf '%s\n' "$dup" | grep -c .) line(s) identical to the jobRecord header beyond the first -- a stray duplicate header, or a data row whose pipeline/step are literally the header words (rename that pipeline). NOT migrated as data; saved to $unrec (full original in $bak)." >&2
    fi
    # Curation pointer: v1 rows have no pipeline column, so jrEmitV2 tags every migrated
    # legacy row pipeline=legacy. Those rows still feed the standalone/legacy estimate
    # bucket, but all pipelines' history is pooled under one "legacy" id, which blunts
    # PER-PIPELINE accuracy. Point the user at the (optional, one-time) curation step so
    # the v2 experience is complete, not a dead end. Fires only when this migration
    # actually produced legacy rows, so a store with none is left silent. The explicit
    # command is named only when jrCurate.sh is actually installed (it ships separately);
    # otherwise the note stays generic so we never point at a tool that is not present.
    if awk -F, 'NR>1 && $1=="legacy"{f=1; exit} END{exit !f}' "$f" 2>/dev/null; then
        if command -v jrCurate.sh >/dev/null 2>&1; then
            echo "SmartSlurm: NOTE: for optimal per-pipeline performance after migrating records, migrated history is tagged pipeline=legacy (v1 stored no pipeline id) and pooled. It still informs estimates, but assigning real pipeline names sharpens per-pipeline accuracy. Optional one-time curation: 'jrCurate.sh export \"$(dirname -- "$f")\"', edit the manifest, then 'jrCurate.sh import <manifest>'." >&2
        else
            echo "SmartSlurm: NOTE: for optimal per-pipeline performance after migrating records, migrated history is tagged pipeline=legacy (v1 stored no pipeline id) and pooled under one id. It still informs estimates, but assigning real per-pipeline names sharpens accuracy -- ask your SmartSlurm maintainer about curating legacy records." >&2
        fi
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
# RETURN CONTRACT: 0 only when the rewrite was APPLIED (atomic rename done) or the
# store was empty (honest no-op); NON-ZERO when it did NOT write -- lock timeout, awk
# error, empty output, or a store that is not (and cannot become) v2. Prune-only
# callers (cleanUp) ignore the status and safely skip on contention; a caller that
# reports success to a user (e.g. a curation import) MUST check it, because a false
# "applied" on a silent skip is a data-integrity trap (adversary Finding 1).
jrRewrite() {
    local f=$1 prog=$2; shift 2
    [ -s "$f" ] || return 0                        # nothing to rewrite -> honest no-op
    _jrLock "$f" || return 1                        # could NOT take the lock -> did NOT write
    local rc=1                                      # "did not apply" until a successful rename
    _jrMigrateUnlocked "$f"
    if [ "$(jrSchemaVersion "$f")" = 2 ]; then
        local tmp="$f.rw.$$"
        if awk -F, -v OFS=, "$@" '
                NR==1 { for(i=1;i<=NF;i++) col[$i]=i; print; next }
                { '"$prog"' }
            ' "$f" > "$tmp" && [ -s "$tmp" ]; then
            mv "$tmp" "$f" && rc=0                  # applied iff the atomic rename succeeded
        else
            rm -f "$tmp"
        fi
    fi
    _jrUnlock "$f"
    return $rc
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
      $1=="pipeline" && $2=="step" { v2h=1; next }           # v2 header line (NOT a data row whose pipeline is named "pipeline")
      $1=="1jobID" { next }                                  # v1 header line
      NF==17 && (v2h || $1 !~ /^[0-9]+$/) { print; next }    # already-v2 data row -> pass through. An all-digit $1 is kept ONLY under a v2 header seen above (else it is an ambiguous corrupt v1 row, skipped below).
      $1 ~ /^[0-9]+$/ && NF>=20 {                             # v1 data row -> map
          # epoch = numeric tail if it looks like a unix timestamp, else none (20-col).
          # Use length(), NOT an {n,} interval regex: mawk (Debian/Ubuntu default awk)
          # ignores interval regexes, which would silently blank EVERY migrated epoch.
          ep = ($NF ~ /^[0-9]+$/ && length($NF) >= 9) ? $NF : ""
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
