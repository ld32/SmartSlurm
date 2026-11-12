#!/bin/bash
# ---------------------------------------------------------------------------
# extract_dag.sh -- print the job DAG that the wf-mode mock sbatch recorded.
#
# Precondition: the user's submit script was run under ss_wf_activate
# (register-don't-run), so $MOCK_SLURM_STATE/wf/<id>.{name,deps,run,array} exist.
# We read that state instead of statically parsing the submit script, so dynamic
# job-id capture, loops, conditionals and arrays are all handled -- the submission
# logic actually ran (no payloads did).
#
# Output (TSV, ids resolved to step names):
#   NODE <name> <single|array[1,2,3]> <command>
#   EDGE <name> <deptype(AND|OR)> <depname>
# A step whose name is reused (fan-out loop) appears once per registered job.
# ---------------------------------------------------------------------------
set -uo pipefail
wf="${MOCK_SLURM_STATE:?set MOCK_SLURM_STATE (activate the harness + run the submit script in wf mode first)}/wf"
[ -d "$wf" ] && ls "$wf"/*.name >/dev/null 2>&1 || {
  echo "extract_dag: no registered workflow in $wf (did you ss_wf_activate and run the submit script?)" >&2
  exit 1
}

ids=$(ls "$wf"/*.name | sed 's#.*/##; s/\.name$//' | sort -n)

# id -> name map (for resolving dependency ids)
declare -A NAME
for id in $ids; do NAME[$id]="$(cat "$wf/$id.name" 2>/dev/null)"; done

for id in $ids; do
  name="${NAME[$id]:-job$id}"
  if [ -f "$wf/$id.array" ]; then kind="array[$(tr ' ' ',' < "$wf/$id.array")]"; else kind="single"; fi
  cmd="$(tail -1 "$wf/$id.run" 2>/dev/null)"        # the payload line (after the cd)
  printf 'NODE\t%s\t%s\t%s\n' "$name" "$kind" "$cmd"

  deps="$(cat "$wf/$id.deps" 2>/dev/null)"
  [ -z "$deps" ] && continue
  if [[ "$deps" == *'?'* ]]; then conj=OR; IFS='?' read -ra terms <<< "$deps"
  else                            conj=AND; IFS=',' read -ra terms <<< "$deps"; fi
  for term in "${terms[@]}"; do
    [ -z "$term" ] && continue
    type="${term%%:*}"
    if [ "$type" = singleton ]; then
      printf 'EDGE\t%s\t%s\t%s\n' "$name" "singleton($conj)" "-"; continue
    fi
    for d in $(echo "${term#*:}" | tr ':' ' '); do
      d="${d%%+*}"; [[ "$d" =~ ^[0-9]+$ ]] || continue
      printf 'EDGE\t%s\t%s\t%s\n' "$name" "$type($conj)" "${NAME[$d]:-job$d}"
    done
  done
done
