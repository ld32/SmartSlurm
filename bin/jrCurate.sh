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
        *[,\ ]*|*"	"*) echo "contains a space/comma/tab" ;;
        legacy|unknown)   echo "is a reserved marker" ;;
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

    local legacy_n
    legacy_n=$(jrEmitV2 "$jr" 2>/dev/null | awk -F, 'NR>1 && $1=="legacy"{c++} END{print c+0}')
    [ "${legacy_n:-0}" -gt 0 ] || { echo "$PROG: no legacy rows to curate in $jr (nothing migrated, or already curated)."; return 0; }

    # Pass 1: cluster + write the CSV (pipeline col = placeholder) + a stderr summary.
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
      }' > "$out"

    # Pass 2: the manifest template. Header from bash (quoting-safe); one data line
    # per candidate from a quote-free awk over the CSV (sorted by placeholder number).
    {
        printf '#jrCurate-manifest-v1\n'
        printf '# Rename migrated pipelines: edit the newName column ONLY -- replace the "-" with a\n'
        printf '# real pipeline name to keep those rows as that pipeline; leave "-" to keep them as\n'
        printf '# legacy. The other columns are read-only hints. Columns are "|"-separated.\n'
        printf '# Then:  %s import %s %s\n#\n' "$PROG" "$man" "$dir"
        printf '# %-12s| %-16s| %5s| %-30s| %s\n' currentName newName rows steps refs
        awk -F, '
          function sortset(c,  m,b,i,j,t){ m=split(c,b,","); for(i=1;i<=m;i++)for(j=i+1;j<=m;j++)if(b[j]<b[i]){t=b[i];b[i]=b[j];b[j]=t} t=b[1];for(i=2;i<=m;i++)if(b[i]!=b[i-1])t=t","b[i]; return t }
          function addset(cur,v){ return (index(","cur",", ","v",")? cur : (cur==""?v:cur","v)) }
          NR==1 || $1=="unknown" { next }
          { ph=$1; cnt[ph]++; st[ph]=addset(st[ph],$2); rf[ph]=addset(rf[ph],$3);
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
        # manifest columns: currentName | newName | rows | steps | refs ("|"-separated).
        # Edit only newName; "-" or empty means "leave as legacy" (skip).
        local phmap; phmap="$(mktemp)"
        awk -F'|' '
            /^#/ || /^[[:space:]]*$/ { next }
            { cn=$1; nn=$2;
              gsub(/^[[:space:]]+/,"",cn); gsub(/[[:space:]]+$/,"",cn);
              gsub(/^[[:space:]]+/,"",nn); gsub(/[[:space:]]+$/,"",nn);
              if(cn=="" || nn=="" || nn=="-") next;
              print cn "\t" nn }
          ' "$file" > "$phmap"
        awk -F, -v pf="$phmap" '
            BEGIN{ while((getline l<pf)>0){ split(l,a,"\t"); nm[a[1]]=a[2] } }
            NR==1{next} ($1 in nm){ print $5 "\t" nm[$1] }
          ' "$csv" > "$map"
        rm -f "$phmap"
    else
        # CSV MODE: the edited pipeline column directly (jobID col5, pipeline col1)
        awk -F, '
            NR==1 { if ($1!="pipeline" || $5!="jobID") { print "BADHEADER" > "/dev/stderr"; exit 3 } ; next }
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

    local ts bak; ts=$(date +%Y-%m-%d_%H-%M-%S); bak="$jr.back.$ts.txt"
    cp -p "$jr" "$bak" || { rm -f "$map"; die "backup failed, aborting"; }
    if jrRewrite "$jr" '
            if(!_ld){ while((getline _l < mapfile)>0){ split(_l,_a,"\t"); want[_a[1]]=_a[2] } close(mapfile); _ld=1 }
            if($(col["pipeline"])=="legacy" && ($(col["jobID"]) in want)) $(col["pipeline"])=want[$(col["jobID"])];
            print
        ' -v mapfile="$map"; then
        echo "$PROG: assigned $n legacy row(s) to a pipeline (backup: $bak)"
        awk -F'\t' '{c[$2]++} END{for(p in c) printf "  %-20s %d row(s)\n", p, c[p]}' "$map"
    else
        cp -p "$bak" "$jr"; rm -f "$map"; die "rewrite failed; store restored from backup"
    fi
    rm -f "$map"
}

case "$sub" in
    export) cmd_export "$@" ;;
    import) cmd_import "$@" ;;
    *) usage 2 ;;
esac
