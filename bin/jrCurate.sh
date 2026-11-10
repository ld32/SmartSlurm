#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# jrCurate.sh -- guided manual curation of migrated legacy jobRecords.
#
# After a v1 store is migrated, every old row carries pipeline=legacy: v1 had no
# pipeline column, so a standalone run and a pipeline step cannot be told apart
# automatically (the tool a step ran is never recorded). This tool lets a human
# assign real pipeline names to that history so it becomes first-class per-pipeline
# training data -- with as much machine guidance as the data allows.
#
#   jrCurate.sh export [jobRecordDir]        # read-only. Clusters the legacy rows
#                                            # into candidate pipelines and writes:
#                                            #   jobRecord.curate.csv       (the full
#                                            #       clustered rows; machine record)
#                                            #   jobRecord.curate.manifest  (a tiny
#                                            #       TEMPLATE: one line per candidate
#                                            #       -- this is what you edit)
#   jrCurate.sh import <file> [jobRecordDir] # apply your names. <file> is either the
#                                            # edited manifest (recommended) or the
#                                            # edited CSV (power users). By jobID,
#                                            # through the locked API, with a backup
#                                            # and an integrity gate.
#
# GUIDED CLUSTERING (export): rows whose flag is rAP-structured
# (step#.dep#.name[.sample...]) are pipeline steps; rows with a free-form flag are
# treated as standalone and left as "unknown" (they keep feeding the standalone
# estimation bucket). Pipeline steps are clustered into candidate pipelines by
# SAMPLE-CHAIN SHAPE: the set of step names that share a sample is one pipeline's
# shape, so two pipelines that merely share a step name (e.g. "mapping") stay
# separate as long as their overall shapes differ. Each distinct shape gets a
# placeholder id (legacyPipeN). The clustering is a guide, not an oracle.
#
# THE MANIFEST is the small edit surface -- one "|"-separated row per candidate, with
# columns:   currentName | newName | rows | steps | refs
#     legacyPipe1   | -       |     3 | align,count,trim | hg38,none
# Edit the newName column only: replace "-" with a real pipeline name to keep those
# rows as that pipeline; leave "-" to keep them as legacy. The other columns are
# read-only hints. Then: jrCurate.sh import <manifest>.
# ---------------------------------------------------------------------------
set -uo pipefail
. "$(dirname "$0")/jobRecordSchema.sh"

PROG=jrCurate
die()  { echo "$PROG: $*" >&2; exit 1; }
warn() { echo "$PROG: $*" >&2; }
usage() { sed -n '2,41p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }

sub="${1:-}"; [ -n "$sub" ] || usage 2
case "$sub" in -h|--help|help) usage 0 ;; esac
shift

_resolve_dir() {
    local d="${1:-}"
    [ -n "$d" ] || d="${smartSlurmJobRecordDir:-}"
    if [ -z "$d" ]; then
        local cfg; cfg="$(dirname "$0")/../config/config.txt"
        [ -f "$cfg" ] && d="$(. "$cfg" 2>/dev/null; echo "${smartSlurmJobRecordDir:-}")"
    fi
    [ -n "$d" ] || die "no jobRecord dir (pass one, or set smartSlurmJobRecordDir)"
    [ -d "$d" ] || die "not a directory: $d"
    printf '%s' "$d"
}

# A chosen pipeline name must be flag/CSV-safe, contain a non-digit (the schema
# drops an all-digit first column), and not be a reserved marker. Echoes a reason if
# bad, nothing if OK.
_badname() {
    case "$1" in
        *'|'*)            echo "contains a pipe (the manifest column separator -- would be truncated)" ;;
        *[,\ ]*|*"	"*) echo "contains a space/comma/tab" ;;
        legacy|unknown)   echo "is a reserved marker" ;;
        legacyPipe[0-9]*) echo "collides with jrCurate's placeholder naming (legacyPipeN); choose a real pipeline name" ;;
        1jobID|pipeline)  echo "collides with a jobRecord header line (the schema would silently drop its rows)" ;;
        *[!0-9]*)         : ;;                       # has a non-digit -> ok
        '')               echo "is empty" ;;
        *)                echo "is all-digits (the schema would drop its rows)" ;;
    esac
}

# ===========================================================================
# export  -- read-only; writes the clustered CSV + the manifest template.
# v2 column order is canonical (pipeline=1 ... jobID=5 ... flag=17), guaranteed by
# jrEmitV2, so the awks address columns by that fixed position.
# ===========================================================================
cmd_export() {
    local dir; dir="$(_resolve_dir "${1:-}")"
    local jr="$dir/jobRecord.txt"
    [ -s "$jr" ] || die "no jobRecord.txt in $dir"
    local out="$dir/jobRecord.curate.csv" man="$dir/jobRecord.curate.manifest"
    # An export-id stamped into BOTH the CSV and the manifest binds them together, so
    # import can refuse a manifest that was edited against a DIFFERENT (e.g. newer)
    # export -- otherwise stale placeholder numbering could rename the wrong rows to a
    # perfectly valid name (silent misassignment).
    local exid; exid="$(date +%Y%m%d%H%M%S)-$$-${RANDOM}"

    local legacy_n
    legacy_n=$(jrEmitV2 "$jr" 2>/dev/null | awk -F, 'NR>1 && $1=="legacy"{c++} END{print c+0}')
    [ "${legacy_n:-0}" -gt 0 ] || { echo "$PROG: no legacy rows to curate in $jr (nothing migrated, or already curated)."; return 0; }

    # Pass 1: cluster + write the CSV (pipeline col = placeholder) + a stderr summary.
    # The CSV leads with the export-id comment (import skips '#' lines on read).
    { printf '#jrCurate-export-id:%s\n' "$exid"
    jrEmitV2 "$jr" 2>/dev/null | awk -F, -v OFS=, -v hdr="$JR_HEADER_V2" '
      function parseFlag(fl,  n,a,i) {
          n=split(fl,a,"."); _isRap=0; _name=""; _sample="";
          if (n>=3 && a[1] ~ /^[0-9]+$/ && a[2] ~ /^[0-9]+$/) _isRap=1; else return 0;
          i=2; while (i<=n && a[i] ~ /^[0-9]+$/) i++;
          _name=(i<=n)?a[i]:""; i++; _sample=(i<=n)?a[i]:""; return 1;
      }
      function sortset(c,  m,b,i,j,t) {
          m=split(c,b,","); for(i=1;i<=m;i++)for(j=i+1;j<=m;j++)if(b[j]<b[i]){t=b[i];b[i]=b[j];b[j]=t}
          t=b[1]; for(i=2;i<=m;i++) if(b[i]!=b[i-1]) t=t","b[i]; return t;
      }
      function addset(cur,v) { return (index(","cur",", ","v",")? cur : (cur==""?v:cur","v)) }
      NR==1 { next }
      $1!="legacy" { next }
      {
          nrow++; line[nrow]=$0; rref[nrow]=$3;
          if (parseFlag($17)) {
              isr[nrow]=1; samp[nrow]=_sample; nm[nrow]=(_name!=""?_name:$2);
              if (_sample!="") sampSteps[_sample]=addset(sampSteps[_sample], nm[nrow]);
          } else { isr[nrow]=0; }
      }
      END {
          for (r=1;r<=nrow;r++) {
              if (!isr[r]) { plan[r]="unknown"; continue }
              if (samp[r]!="" && sampSteps[samp[r]]!="") shape=sortset(sampSteps[samp[r]]); else shape=nm[r];
              if (!(shape in shapeId)) { shapeId[shape]="legacyPipe" (++nshape); shapeOf[nshape]=shape }
              plan[r]=shapeId[shape];
          }
          print hdr;
          for (r=1;r<=nrow;r++) {
              n=split(line[r],f,","); f[1]=plan[r];
              s=f[1]; for (i=2;i<=n;i++) s=s","f[i]; print s;
              cnt[plan[r]]++; refs[plan[r]]=addset(refs[plan[r]], rref[r]);
          }
          printf "jrCurate: %d legacy row(s); %d candidate pipeline shape(s)\n", nrow, nshape > "/dev/stderr";
          for (k=1;k<=nshape;k++) { id="legacyPipe" k;
              printf "  %-12s steps={%s} refs={%s} rows=%d\n", id, shapeOf[k], refs[id], cnt[id] > "/dev/stderr" }
          if ("unknown" in cnt)
              printf "  %-12s (standalone / free-form flag -- kept as the shared standalone bucket) rows=%d\n", "unknown", cnt["unknown"] > "/dev/stderr";
      }'
    } > "$out"

    # Pass 2: the manifest template. Header from bash (quoting-safe); one data line
    # per candidate from a quote-free awk over the CSV (sorted by placeholder number).
    {
        printf '#jrCurate-manifest-v1 export-id:%s\n' "$exid"
        printf '# Rename migrated pipelines: edit the newName column ONLY -- replace the "-" with a\n'
        printf '# real pipeline name to keep those rows as that pipeline; leave "-" to keep them as\n'
        printf '# legacy. The other columns are read-only hints. Columns are "|"-separated.\n'
        printf '# Then:  %s import %s %s\n#\n' "$PROG" "$man" "$dir"
        printf '# %-12s| %-16s| %5s| %-30s| %s\n' currentName newName rows steps refs
        awk -F, '
          function sortset(c,  m,b,i,j,t){ m=split(c,b,","); for(i=1;i<=m;i++)for(j=i+1;j<=m;j++)if(b[j]<b[i]){t=b[i];b[i]=b[j];b[j]=t} t=b[1];for(i=2;i<=m;i++)if(b[i]!=b[i-1])t=t","b[i]; return t }
          function addset(cur,v){ return (index(","cur",", ","v",")? cur : (cur==""?v:cur","v)) }
          /^#/ || $1=="pipeline" || $1=="unknown" { next }
          # Sanitize pipe chars out of the read-only hint columns (steps/refs): the
          # manifest is pipe-separated, so a stored value that legitimately contains a
          # pipe (a ref or tool name) would otherwise push the line past 5 fields and
          # trip the NF>5 guard on import, falsely rejecting a line the user never
          # edited. Hints are cosmetic, so mapping pipe to underscore here is harmless;
          # import keys on jobID, not on these columns.
          { ph=$1; s2=$2; s3=$3; gsub(/\|/,"_",s2); gsub(/\|/,"_",s3);
            cnt[ph]++; st[ph]=addset(st[ph],s2); rf[ph]=addset(rf[ph],s3);
            if(!(ph in seen)){ seen[ph]=1; ord[++k]=ph } }
          END{ for(i=1;i<=k;i++)for(j=i+1;j<=k;j++){ na=ord[i];gsub(/[^0-9]/,"",na); nb=ord[j];gsub(/[^0-9]/,"",nb); if(nb+0<na+0){t=ord[i];ord[i]=ord[j];ord[j]=t} }
               for(i=1;i<=k;i++){ ph=ord[i]; printf "%-12s| %-16s| %5d| %-30s| %s\n", ph, "-", cnt[ph], sortset(st[ph]), sortset(rf[ph]) } }
        ' "$out"
    } > "$man"

    echo "$PROG: wrote $out"
    echo "$PROG: wrote $man   <-- EDIT THIS (one line per candidate; fill names after '=')"
    echo "$PROG: then:  $PROG import \"$man\" \"$dir\""
}

# ===========================================================================
# import  -- apply names by jobID; backup + integrity gate. Accepts the manifest
# (recommended) or a directly-edited CSV.
# ===========================================================================
cmd_import() {
    local file="${1:-}"; [ -n "$file" ] || die "usage: $PROG import <manifest|csv> [jobRecordDir]"
    [ -f "$file" ] || die "file not found: $file"
    local dir; dir="$(_resolve_dir "${2:-}")"
    local jr="$dir/jobRecord.txt"
    [ -s "$jr" ] || die "no jobRecord.txt in $dir"

    local map; map="$(mktemp)"     # jobID<TAB>chosen-pipeline

    if head -1 "$file" 2>/dev/null | grep -q '^#jrCurate-manifest'; then
        # MANIFEST MODE: placeholder->name (manifest) + jobID->placeholder (CSV)
        local csv; csv="$(dirname "$file")/jobRecord.curate.csv"
        [ -f "$csv" ] || csv="$dir/jobRecord.curate.csv"
        [ -f "$csv" ] || { rm -f "$map"; die "manifest mode needs jobRecord.curate.csv next to the manifest (re-run export)"; }
        # Bind the manifest to the CSV it was generated from: refuse a manifest edited
        # against a DIFFERENT export. Placeholder numbering (legacyPipeN) can change
        # between exports, so a stale manifest + fresh CSV would silently rename the
        # WRONG rows to a perfectly valid name -- the hardest kind of error to notice.
        local man_id csv_id
        man_id=$(sed -n '1{s/.*export-id:\([^[:space:]]*\).*/\1/p;q}' "$file")
        csv_id=$(sed -n '1{s/^#jrCurate-export-id:\([^[:space:]]*\).*/\1/p;q}' "$csv")
        [ -n "$man_id" ] && [ -n "$csv_id" ] || { rm -f "$map"; die "manifest/CSV predate export-id tagging; re-run '$PROG export' and edit the fresh manifest"; }
        [ "$man_id" = "$csv_id" ] || { rm -f "$map"; die "manifest (export-id $man_id) does not match jobRecord.curate.csv (export-id $csv_id): the manifest is from a different export. Re-run '$PROG export' and re-edit, or pair the manifest with its own CSV."; }
        # manifest columns: currentName | newName | rows | steps | refs ("|"-separated).
        # Edit only newName; "-" or empty means "leave as legacy" (skip). A currentName
        # must appear at most once -- otherwise a duplicate is silently last-wins.
        local phmap; phmap="$(mktemp)"
        awk -F'|' '
            /^#/ || /^[[:space:]]*$/ { next }
            # The template has exactly 5 "|"-columns. More than 5 means a value (almost
            # always an edited newName) contains the "|" separator, which would silently
            # truncate the name -- reject it rather than apply the fragment.
            NF>5 { printf "jrCurate: malformed manifest line (a value contains the '"'"'|'"'"' separator): %s\n", $0 > "/dev/stderr"; bad=1; next }
            { cn=$1; nn=$2;
              gsub(/^[[:space:]]+/,"",cn); gsub(/[[:space:]]+$/,"",cn);
              gsub(/^[[:space:]]+/,"",nn); gsub(/[[:space:]]+$/,"",nn);
              if(cn=="") next;
              if(cn in seen){ printf "jrCurate: duplicate currentName in manifest: %s\n", cn > "/dev/stderr"; dup=1; next }
              seen[cn]=1;
              if(nn=="" || nn=="-") next;
              print cn "\t" nn }
            END{ if(dup||bad) exit 7 }
          ' "$file" > "$phmap" || { rm -f "$map" "$phmap"; die "manifest has duplicate or malformed lines (see above); each candidate appears once, and no value may contain '|'"; }
        # CSV read skips the leading '#export-id' comment; the data header
        # (pipeline,step,...) never matches a placeholder name, so it is ignored too.
        awk -F, -v pf="$phmap" '
            BEGIN{ while((getline l<pf)>0){ split(l,a,"\t"); nm[a[1]]=a[2] } }
            /^#/ { next } ($1 in nm){ print $5 "\t" nm[$1] }
          ' "$csv" > "$map"
        rm -f "$phmap"
    else
        # CSV MODE: the edited pipeline column directly (jobID col5, pipeline col1).
        # Skip the leading '#export-id' comment, then the first real line is the header.
        awk -F, '
            /^#/ { next }
            !hdr { hdr=1; if ($1!="pipeline" || $5!="jobID") { print "BADHEADER" > "/dev/stderr"; exit 3 } ; next }
            NF<17 { next }
            { p=$1; j=$5;
              if (j !~ /^[0-9]+$/) next;
              if (p=="" || p=="legacy" || p=="unknown") next;
              if (p ~ /^legacyPipe[0-9]+$/) next;
              print j "\t" p }
          ' "$file" > "$map" || { rm -f "$map"; die "invalid CSV (bad header)"; }
    fi

    # Count assignments. `grep -c` prints 0 AND exits non-zero on an empty map, so a
    # `|| echo 0` fallback would append a second "0" and make $n a two-line string that
    # breaks the numeric test below; capture stdout only and default an empty result.
    local n; n=$(grep -c . "$map" 2>/dev/null); n=${n:-0}
    [ "${n:-0}" -gt 0 ] || { rm -f "$map"; echo "$PROG: nothing to import (no candidate was given a real pipeline name)."; return 0; }

    # validate every chosen name (flag/CSV-safe, has a non-digit, not reserved)
    local nm r bad=""
    while IFS= read -r nm; do r=$(_badname "$nm"); [ -n "$r" ] && { bad="'$nm' $r"; break; }; done \
        < <(awk -F'\t' '{print $2}' "$map" | sort -u)
    [ -z "$bad" ] || { rm -f "$map"; die "invalid pipeline name: $bad"; }

    # integrity: every assigned jobID must be a legacy row in the store
    local miss
    miss=$(jrEmitV2 "$jr" 2>/dev/null | awk -F, -v mf="$map" '
        BEGIN{ while((getline l < mf)>0){ split(l,a,"\t"); want[a[1]]=1 } }
        NR>1 && $1=="legacy" { have[$5]=1 }
        END{ m=0; for(j in want) if(!(j in have)) m++; print m+0 }')
    [ "${miss:-0}" = 0 ] || { rm -f "$map"; die "$miss assigned jobID(s) are not legacy rows in the store (refusing to touch non-legacy history)"; }

    local ts bak; ts=$(date +%Y-%m-%d_%H-%M-%S); bak="$jr.back.$ts.$$.txt"
    cp -p "$jr" "$bak" || { rm -f "$map"; die "backup failed, aborting"; }
    if jrRewrite "$jr" '
            if(!_ld){ while((getline _l < mapfile)>0){ split(_l,_a,"\t"); want[_a[1]]=_a[2] } close(mapfile); _ld=1 }
            if($(col["pipeline"])=="legacy" && ($(col["jobID"]) in want)) $(col["pipeline"])=want[$(col["jobID"])];
            print
        ' -v mapfile="$map"; then
        # Post-apply integrity gate: a 0 exit from jrRewrite is necessary but NOT
        # sufficient. Re-read the store and confirm EVERY assigned jobID now carries its
        # chosen name -- none still reads "legacy", none vanished.
        # Scope the gate to LEGACY-ORIGIN rows only. We relabel legacy rows by jobID,
        # but a non-legacy row can legitimately SHARE a jobID (the v1 counter resets, or
        # a merged/multi-cluster store) and must keep its own name. So, per assigned
        # jobID, success means: no row for it still reads "legacy" (every legacy row got
        # relabeled) AND some row for it now carries the new name (nothing vanished).
        # Keying "$1 must equal the new name for EVERY row with this jobID" would wrongly
        # flag that unrelated same-jobID row and fail a correct import.
        local bad
        bad=$(jrEmitV2 "$jr" 2>/dev/null | awk -F, -v mf="$map" '
            BEGIN{ while((getline l<mf)>0){ split(l,a,"\t"); want[a[1]]=a[2] } }
            NR>1 && ($5 in want){
                if($1=="legacy")      stillLegacy[$5]=1    # an assigned legacy row was NOT relabeled
                if($1==want[$5])      gotName[$5]=1        # some row for this jobID now carries the new name
            }
            END{ bad=0; for(j in want){ if(j in stillLegacy) bad++; else if(!(j in gotName)) bad++ } print bad+0 }')
        if [ "${bad:-1}" != 0 ]; then
            # The relabel committed but verification could not confirm it. Do NOT
            # auto-restore: a blind "cp backup -> store" would clobber any row a
            # concurrent (lock-free) cleanUp append committed after our backup, and would
            # also undo the migration jrRewrite may have performed. Leave the store as-is,
            # keep the pre-curation backup, and fail loudly so a human inspects.
            rm -f "$map"
            die "post-apply check failed: $bad assigned row(s) not correctly relabeled; store left as-is -- inspect it; pre-curation backup: $bak"
        fi
        echo "$PROG: assigned $n legacy row(s) to a pipeline (backup: $bak)"
        awk -F'\t' '{c[$2]++} END{for(p in c) printf "  %-20s %d row(s)\n", p, c[p]}' "$map"
    else
        # jrRewrite returns non-zero ONLY when it did not commit (lock timeout, awk error,
        # or a non-v2 store). It is atomic (temp+rename), so it NEVER leaves a partial
        # store -- there is nothing to restore. Critically we must NOT "cp backup ->
        # store" here: a concurrent, lock-free cleanUp append (jrAppendLine appends even
        # when it cannot take the lock, by design) may have committed a record after our
        # backup was taken; clobbering with the stale backup would silently delete it.
        # Drop the backup only when it is byte-identical to the live store; otherwise the
        # store moved under us, so keep the backup and name it.
        cmp -s "$jr" "$bak" && rm -f "$bak"
        rm -f "$map"
        [ -f "$bak" ] && die "curation did not apply (lock contention or error); store left as-is (pre-curation backup: $bak)" \
                      || die "curation did not apply (lock contention or error); store left as-is"
    fi
    rm -f "$map"
}

case "$sub" in
    export) cmd_export "$@" ;;
    import) cmd_import "$@" ;;
    *) usage 2 ;;
esac
