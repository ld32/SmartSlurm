#!/bin/bash
# jobRecordsOutliers.sh -- move outlier records out of jobRecord.txt (never deletes).
#
# Usage:
#   jobRecordsOutliers.sh [--dry-run] <pipeline> <step> <ref> [jobRecord.txt]
#   jobRecordsOutliers.sh [--dry-run] --all [jobRecord.txt]     # every group
#
# A COMPLETED record is an outlier when its memory or time use is far from the
# rest of its group (same pipeline + step + ref), BOTH too big and too small:
#   - compared per unit of input (memUsed/inputSize, timeUsed/inputSize) so big
#     inputs using more is normal; records without an input size are compared on
#     raw memUsed/timeUsed, separately;
#   - outlier if   |x - median| > K * 1.4826 * MAD     (robust z-score > K)
#            and   x > F * median  or  x < median / F   (clearly different)
#            and   the absolute difference is > JOB_RECORDS_OUTLIER_MIN_MEM M / _MIN_TIME min
#     (MAD = median absolute deviation; the absolute floors stop 1-vs-4 minute or
#      20-vs-80 MB jobs from being flagged as noise);
#   - only groups with >= JOB_RECORDS_OUTLIER_MIN_N records are checked.
# Tunables (env): JOB_RECORDS_OUTLIER_K=5 JOB_RECORDS_OUTLIER_F=3 JOB_RECORDS_OUTLIER_MIN_N=5
#                 JOB_RECORDS_OUTLIER_MIN_MEM=200 JOB_RECORDS_OUTLIER_MIN_TIME=5  JOB_RECORDS_OUTLIER=off disables.
#
# Moved rows are appended to jobRecords.outliers.txt with the reason added to the
# notes column; every move is also logged to outliers.log next to it. To undo,
# copy the line back into jobRecord.txt (and delete its "outlier:" note).
# Prints one human-readable line per outlier on stdout (nothing if none).

. "$(dirname "$0")/jobRecordSchema.sh" 2>/dev/null || { echo "jobRecordsOutliers: cannot load jobRecordSchema.sh" >&2; exit 0; }

[ "${JOB_RECORDS_OUTLIER:-on}" = off ] && exit 0

dry=""; [ "$1" = --dry-run ] && { dry=1; shift; }
if [ "$1" = --all ]; then
    all=1; jobRecords=${2:-$smartSlurmJobRecordDir/jobRecord.txt}
else
    all=""; p=$1; a=$2; b=$3; jobRecords=${4:-$smartSlurmJobRecordDir/jobRecord.txt}
    [ -z "$b" ] && { sed -n '2,12p' "$0" | sed 's/^# \{0,1\}//'; exit 1; }
fi
[ -s "$jobRecords" ] || exit 0
outF="${jobRecords%.txt}.outliers.txt"; logF="$(dirname "$jobRecords")/outliers.log"

# Pass 1 (read-only): find outliers. Output: jobID<TAB>reason<TAB>human message
found=$(jobRecordsEmitV2 "$jobRecords" 2>/dev/null | awk -F, -v all="$all" -v p="$p" -v a="$a" -v b="$b" \
    -v K="${JOB_RECORDS_OUTLIER_K:-5}" -v F="${JOB_RECORDS_OUTLIER_F:-3}" -v N="${JOB_RECORDS_OUTLIER_MIN_N:-5}" \
    -v MM="${JOB_RECORDS_OUTLIER_MIN_MEM:-200}" -v MT="${JOB_RECORDS_OUTLIER_MIN_TIME:-5}" "$JOB_RECORDS_AWK_COLMAP"'
    function med(arr, n,   i, j, t, s) {   # median of arr[1..n] (sorts a copy)
        for (i = 1; i <= n; i++) s[i] = arr[i]
        for (i = 2; i <= n; i++) { t = s[i]; for (j = i - 1; j >= 1 && s[j] > t; j--) s[j+1] = s[j]; s[j+1] = t }
        return (n % 2) ? s[(n+1)/2] : (s[n/2] + s[n/2+1]) / 2
    }
    $(col["status"]) == "COMPLETED" && (all || ($(col["pipeline"]) == p && $(col["step"]) == a && $(col["ref"]) == b)) {
        sz = $(col["inputSize"]); hasIn = (sz ~ /^[0-9]+$/ && sz + 0 > 0)
        g = $(col["pipeline"]) SUBSEP $(col["step"]) SUBSEP $(col["ref"]) SUBSEP hasIn
        k = ++n[g]
        id[g, k] = $(col["jobID"]); size[g, k] = hasIn ? sz : 1
        v["mem", g, k] = $(col["memUsed"]) + 0; v["time", g, k] = $(col["timeUsed"]) + 0
        label[g] = $(col["pipeline"]) "/" $(col["step"]) "/" $(col["ref"]) (hasIn ? "" : " (no input size)")
    }
    END {
        for (g in n) {
            if (n[g] < N) continue
            split(g, gp, SUBSEP); hasIn = gp[4]
            for (mi = 1; mi <= 2; mi++) {
                m = (mi == 1) ? "mem" : "time"; floor = (m == "mem") ? MM : MT; unit = (m == "mem") ? "M" : "min"
                for (k = 1; k <= n[g]; k++) x[k] = v[m, g, k] / size[g, k]
                md = med(x, n[g])
                for (k = 1; k <= n[g]; k++) d[k] = (x[k] > md) ? x[k] - md : md - x[k]
                mad = 1.4826 * med(d, n[g])
                for (k = 1; k <= n[g]; k++) {
                    if (md <= 0) continue
                    far   = (mad > 0) ? (d[k] > K * mad) : (d[k] > 0)
                    ratio = x[k] / md
                    diff  = (x[k] > md ? x[k] - md : md - x[k]) * size[g, k]   # back in M / min
                    if (far && (ratio > F || ratio < 1 / F) && diff > floor) {
                        want = md * size[g, k]
                        dir = (ratio > 1) ? "too-big" : "too-small"
                        r = sprintf("outlier:%s:%s:%.2fx", m, dir, ratio)
                        msg = sprintf("%s job %s: %s used %d%s, expected ~%d%s (%.2fx the group median, %d records)", \
                            label[g], id[g, k], m, v[m, g, k], unit, want, unit, ratio, n[g])
                        if (!((g, id[g, k]) in seen)) { seen[g, id[g, k]] = 1; print id[g, k] "\t" r "\t" msg }
                    }
                }
            }
        }
    }')

[ -z "$found" ] && exit 0
printf '%s\n' "$found" | cut -f3 | sed 's/^/Outlier removed from jobRecords: /'
[ -n "$dry" ] && { echo "(dry run: nothing moved)"; exit 0; }

# Pass 2 (locked, atomic): move the flagged rows, reason appended to notes.
ids=$(printf '%s\n' "$found" | cut -f1,2 | tr '\t\n' '= ' )
[ -f "$outF" ] || printf '%s\n' "$JOB_RECORDS_HEADER_V2" > "$outF"
jobRecordsRewrite "$jobRecords" '
    if (!init) { init = 1; nk = split(ids, kv, " "); for (i = 1; i <= nk; i++) { split(kv[i], q, "="); why[q[1]] = q[2] } }
    j = $(col["jobID"])
    if ($(col["status"]) == "COMPLETED" && (j in why) && (all || ($(col["pipeline"]) == p && $(col["step"]) == a && $(col["ref"]) == b))) {
        $(col["notes"]) = ($(col["notes"]) == "" ? "" : $(col["notes"]) ";") why[j]
        print >> of
        next
    }
    print' -v ids="$ids" -v of="$outF" -v all="$all" -v p="$p" -v a="$a" -v b="$b"

printf '%s\n' "$found" | cut -f3 | sed "s|^|$(date +%Y-%m-%dT%H:%M:%S) moved to $outF: |" >> "$logF"
exit 0
