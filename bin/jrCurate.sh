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
#   jrCurate.sh export [jobRecordDir]        # write an editable, clustered CSV of
#                                            # the legacy rows (pipeline pre-filled
#                                            # with placeholders to rename). Read-only.
#   jrCurate.sh import <file> [jobRecordDir] # merge your edited CSV back into the
#                                            # store, assigning each row's pipeline
#                                            # by jobID. Backs up + integrity-gates.
#
# GUIDED CLUSTERING (export): rows whose flag is rAP-structured
# (step#.dep#.name[.sample...]) are pipeline steps; rows with a free-form flag are
# treated as standalone and left as "unknown" (they keep feeding the standalone
# estimation bucket). Pipeline steps are clustered into candidate pipelines by
# SAMPLE-CHAIN SHAPE: the set of step names that share a sample is one pipeline's
# shape, so two pipelines that merely share a step name (e.g. "mapping") stay
# separate as long as their overall shapes differ. Each distinct shape gets a
# placeholder id (legacyPipeN); you rename the placeholders to real names, using the
# flag/ref/sample/date columns as context, then re-import. The clustering is a guide,
# not an oracle -- you have the final say.
# ---------------------------------------------------------------------------
set -uo pipefail
. "$(dirname "$0")/jobRecordSchema.sh"

PROG=jrCurate
die()  { echo "$PROG: $*" >&2; exit 1; }
warn() { echo "$PROG: $*" >&2; }
usage() { sed -n '2,27p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }

sub="${1:-}"; [ -n "$sub" ] || usage 2
case "$sub" in -h|--help|help) usage 0 ;; esac
shift

# --- resolve the store dir (same convention as jrRepair) ----------------------
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

# ===========================================================================
# export  -- read-only; emits $dir/jobRecord.curate.csv + a summary on stderr.
# v2 column order is canonical (pipeline=1 step=2 ref=3 sample=4 jobID=5 ... flag=17),
# guaranteed by jrEmitV2, so the awks below address columns by that fixed position.
# ===========================================================================
cmd_export() {
    local dir; dir="$(_resolve_dir "${1:-}")"
    local jr="$dir/jobRecord.txt"
    [ -s "$jr" ] || die "no jobRecord.txt in $dir"
    local out="$dir/jobRecord.curate.csv"

    local legacy_n
    legacy_n=$(jrEmitV2 "$jr" 2>/dev/null | awk -F, 'NR>1 && $1=="legacy"{c++} END{print c+0}')
    [ "${legacy_n:-0}" -gt 0 ] || { echo "$PROG: no legacy rows to curate in $jr (nothing migrated, or already curated)."; return 0; }

    jrEmitV2 "$jr" 2>/dev/null | awk -F, -v OFS=, -v hdr="$JR_HEADER_V2" '
      # rAP flag grammar: step#.dep#[.dep#...].name[.sample...]. isRap = starts with
      # num.num; stepName = first non-numeric dotted field; sample = the field after.
      function parseFlag(fl,  n,a,i) {
          n=split(fl,a,"."); _isRap=0; _name=""; _sample="";
          if (n>=3 && a[1] ~ /^[0-9]+$/ && a[2] ~ /^[0-9]+$/) _isRap=1; else return 0;
          i=2; while (i<=n && a[i] ~ /^[0-9]+$/) i++;     # skip step# + numeric deps
          _name=(i<=n)?a[i]:""; i++; _sample=(i<=n)?a[i]:""; return 1;
      }
      function sortset(csv,  m,b,i,j,t) {                 # sort a comma-set, dedup
          m=split(csv,b,","); for(i=1;i<=m;i++)for(j=i+1;j<=m;j++)if(b[j]<b[i]){t=b[i];b[i]=b[j];b[j]=t}
          t=b[1]; for(i=2;i<=m;i++) if(b[i]!=b[i-1]) t=t","b[i]; return t;
      }
      NR==1 { next }                                      # jrEmitV2 header
      $1!="legacy" { next }                               # only legacy rows
      {
          nrow++; line[nrow]=$0; rref[nrow]=$3;
          if (parseFlag($17)) {
              isr[nrow]=1; samp[nrow]=_sample; nm[nrow]=(_name!=""?_name:$2);
              if (_sample!="") {
                  if (index(","sampSteps[_sample]",", ","nm[nrow]",")==0)
                      sampSteps[_sample]=(sampSteps[_sample]==""?nm[nrow]:sampSteps[_sample]","nm[nrow]);
              }
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
              cnt[plan[r]]++;
              if (index(","refs[plan[r]]",", ","rref[r]",")==0)
                  refs[plan[r]]=(refs[plan[r]]==""?rref[r]:refs[plan[r]]","rref[r]);
          }
          printf "jrCurate: %d legacy row(s); %d candidate pipeline shape(s)\n", nrow, nshape > "/dev/stderr";
          for (k=1;k<=nshape;k++) { id="legacyPipe" k;
              printf "  %-12s steps={%s} refs={%s} rows=%d\n", id, shapeOf[k], refs[id], cnt[id] > "/dev/stderr" }
          if ("unknown" in cnt)
              printf "  %-12s (standalone / free-form flag -- kept as the shared standalone bucket) rows=%d\n", "unknown", cnt["unknown"] > "/dev/stderr";
      }' > "$out"

    echo "$PROG: wrote $out"
    echo "$PROG: rename the legacyPipeN placeholders in the 'pipeline' column to real"
    echo "         names (leave 'unknown' rows as-is to keep them standalone), then:"
    echo "         $PROG import \"$out\" \"$dir\""
}

# ===========================================================================
# import  -- merge chosen pipeline names back by jobID; backup + integrity gate.
# ===========================================================================
cmd_import() {
    local file="${1:-}"; [ -n "$file" ] || die "usage: $PROG import <file> [jobRecordDir]"
    [ -f "$file" ] || die "curation file not found: $file"
    local dir; dir="$(_resolve_dir "${2:-}")"
    local jr="$dir/jobRecord.txt"
    [ -s "$jr" ] || die "no jobRecord.txt in $dir"

    # Build jobID<TAB>pipeline from the edited file. Skip the header, rows left as a
    # placeholder / legacy / unknown / empty, and validate the header + name safety.
    local map; map="$(mktemp)"
    awk -F, '
        NR==1 { if ($1!="pipeline" || $5!="jobID") { print "BADHEADER" > "/dev/stderr"; exit 3 } ; next }
        NF<17 { next }
        { p=$1; j=$5;
          if (j !~ /^[0-9]+$/) next;
          if (p=="" || p=="legacy" || p=="unknown") next;
          if (p ~ /^legacyPipe[0-9]+$/) next;                 # placeholder untouched
          if (p ~ /[, \t]/)  { print "BADNAME " p > "/dev/stderr"; exit 4 }
          print j "\t" p }
      ' "$file" > "$map" || { rm -f "$map"; die "invalid curation file (bad header, or unsafe pipeline name)"; }

    local n; n=$(grep -c . "$map" 2>/dev/null || echo 0)
    [ "${n:-0}" -gt 0 ] || { rm -f "$map"; echo "$PROG: nothing to import (no rows assigned a real pipeline name)."; return 0; }

    # Integrity: every assigned jobID must exist AND be a legacy row (we never touch
    # already-named history). Refuse the whole import otherwise.
    local miss
    miss=$(jrEmitV2 "$jr" 2>/dev/null | awk -F, -v mf="$map" '
        BEGIN{ while((getline l < mf)>0){ split(l,a,"\t"); want[a[1]]=1 } }
        NR>1 && $1=="legacy" { have[$5]=1 }
        END{ m=0; for(j in want) if(!(j in have)) m++; print m+0 }')
    [ "${miss:-0}" = 0 ] || { rm -f "$map"; die "$miss assigned jobID(s) are not legacy rows in the store (refusing to touch non-legacy history)"; }

    # Backup, then reassign via the locked rewrite API. The map is loaded once inside
    # the per-row program (jrRewrite supplies col[] and runs the program per row).
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
