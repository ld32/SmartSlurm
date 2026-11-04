#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# jrRepair.sh -- one-time, MANUAL converter of a legacy (v1) jobRecord store to
# the v2 schema. Run once per store after upgrading SmartSlurm.
#
#   jrRepair.sh [jobRecordDir]     # default: $smartSlurmJobRecordDir
#   jrRepair.sh --dry-run [dir]    # report what would change; write nothing
#   jrRepair.sh --check [dir]      # health report only (schema + row census +
#                                  # unrecognized rows); writes nothing, exits
#                                  # non-zero if any store has unrecognized rows
#
# It converts jobRecord.txt and jobRecord.oom.oot.over.txt in place, after
# backing each up to <name>.back.<timestamp>.txt (never deleted -- that backup
# is your rollback). Idempotent: a file already in v2 is left untouched.
#
# v1 -> v2 mapping (by v1 position; see jobRecordSchema.sh):
#   pipeline = "legacy"   (old rows carry no pipeline id)
#   step=program(12) ref(13) sample=from flag(14) jobID(1) inputSize(2)
#   memUsed(7) memAllocated(5) timeUsed(8) timeAllocated(6) cores(15)
#   status(9) user(10) node=""  epoch(21)  isoDate=old date(20)  flag(14)
# Dropped: memDefault timeDefault saccMem extraMem extraTime ratioM ratioT.
# Rows that are not clean 21-field v1 rows are NOT converted: they are copied to
# <name>.jrRepair.skipped.<ts> and counted, so nothing is silently lost.
# ---------------------------------------------------------------------------
set -uo pipefail
. "$(dirname "$0")/jobRecordSchema.sh"

dryrun=""
check=""
dir=""
for a in "$@"; do
    case "$a" in
        --dry-run) dryrun=1 ;;
        --check)   check=1 ;;
        -h|--help) sed -n '2,17p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        -*) echo "jrRepair: unknown option: $a" >&2; exit 2 ;;
        *)  dir=$a ;;
    esac
done
[ -n "$dir" ] || dir="${smartSlurmJobRecordDir:-}"
if [ -z "$dir" ]; then
    # fall back to the configured default if we can read it
    cfg="$(dirname "$0")/../config/config.txt"
    [ -f "$cfg" ] && dir="$(. "$cfg" 2>/dev/null; echo "${smartSlurmJobRecordDir:-}")"
fi
[ -n "$dir" ] || { echo "jrRepair: no jobRecord dir (pass one, or set smartSlurmJobRecordDir)" >&2; exit 2; }
[ -d "$dir" ] || { echo "jrRepair: not a directory: $dir" >&2; exit 2; }

# --check: read-only health report (schema version + row census + unrecognized
# rows) for every store in the dir. Writes nothing. Exits non-zero if any store
# has unrecognized rows, so it can gate CI / a cluster audit.
if [ -n "$check" ]; then
    echo "jrRepair: health check of $dir"
    rc=0
    for s in jobRecord.txt jobRecord.oom.oot.over.txt; do
        [ -e "$dir/$s" ] || continue
        jrLint "$dir/$s" || rc=1
    done
    [ "$rc" = 0 ] && echo "jrRepair: all stores clean." || echo "jrRepair: some stores have unrecognized rows (see above)."
    exit $rc
fi

# convertOne <file> : backs up and rewrites one v1 file to v2 (honors $dryrun).
convertOne() {
    local f=$1 ver ts bak out skip conv total
    [ -f "$f" ] || return 0
    ver=$(jrSchemaVersion "$f")
    case "$ver" in
        2) echo "  $f: already v2 -- skipped"; return 0 ;;
        0) echo "  $f: unrecognized header -- NOT touched (inspect manually)"; return 0 ;;
    esac
    ts=$(date +%Y-%m-%d_%H-%M-%S)
    out="$f.jrRepair.new.$$"
    skip="$f.jrRepair.skipped.$ts"

    # Divert only rows jrEmitV2 actually drops (genuine junk): NOT a header, NOT an
    # already-v2 row (17 fields), NOT a v1 data row (numeric jobID, 20+ fields --
    # both the 20-col no-epoch and 21-col epoch variants convert).
    # (The old test diverted every non-21-field row, which wrongly skipped v2 rows.)
    awk -F, 'NF>0 && $1!="pipeline" && $1!="1jobID" \
             && !($1 !~ /^[0-9]+$/ && NF==17) \
             && !($1 ~ /^[0-9]+$/ && NF>=20)' "$f" > "$skip"
    skip_n=$(wc -l < "$skip"); skip_n=${skip_n//[^0-9]/}
    [ "${skip_n:-0}" -gt 0 ] || rm -f "$skip"

    # Convert via the shared shim (single source of the v1->v2 mapping).
    jrEmitV2 "$f" > "$out"
    conv=$(( $(wc -l < "$out") - 1 ))   # minus the header line
    total=$((conv + ${skip_n:-0}))

    if [ -n "$dryrun" ]; then
        echo "  $f: v1 -> v2 (dry-run): $conv rows convert, ${skip_n:-0} would divert, (of $total)"
        rm -f "$out" "$skip"          # dry-run writes nothing
        return 0
    fi

    bak="$f.back.$ts.txt"
    cp -p "$f" "$bak" || { echo "  $f: backup failed, aborting this file" >&2; rm -f "$out"; return 1; }
    mv "$out" "$f"
    echo "  $f: converted $conv rows to v2 (backup: $bak)"
    [ "$skip_n" -gt 0 ] && echo "    note: $skip_n non-standard row(s) diverted to $skip (not lost)"
    [ -s "$skip" ] || rm -f "$skip"
    return 0
}

echo "jrRepair: ${dryrun:+(dry-run) }scanning $dir"
rc=0
convertOne "$dir/jobRecord.txt"              || rc=1
convertOne "$dir/jobRecord.oom.oot.over.txt" || rc=1
echo "jrRepair: done."
exit $rc
