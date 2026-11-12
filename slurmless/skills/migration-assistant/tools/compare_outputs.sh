#!/bin/bash
# ---------------------------------------------------------------------------
# compare_outputs.sh <origDir> <rapDir> <contract>
#
# Differential test for a port: compare the AUTHORITATIVE outputs of two runs
# (the original workflow vs its rAP port) under an identity contract. The contract
# is what the agent proposes and the user approves — which files are authoritative,
# and how to normalize volatile fields (timestamps, job ids, absolute paths,
# ordering) that legitimately differ between runs.
#
# Contract lines (see references/identity-contract.md):
#   FILE <glob>        an authoritative output to compare (repeatable; globs ok)
#   NORMALIZE <sed>    an ERE sed expression applied to every file before diffing
#                      (repeatable), e.g.  s/job [0-9]+/job <ID>/g
#   SORT <glob>        compare these files order-insensitively (repeatable)
#
# Prints "MATCH <file>" / "DIFF <file>" (+ the diff) per authoritative file, then
# "RESULT: MATCH" or "RESULT: MISMATCH". Exit 0 iff everything matches. A file
# present in only one run is a DIFF.
# ---------------------------------------------------------------------------
set -uo pipefail
od="${1:?usage: compare_outputs.sh <origDir> <rapDir> <contract>}"
rd="${2:?rapDir}"; contract="${3:?contract}"
[ -d "$od" ] && [ -d "$rd" ] && [ -f "$contract" ] || { echo "compare_outputs: bad args" >&2; exit 2; }

files=(); norms=(); sorts=()
while IFS= read -r line || [ -n "$line" ]; do
  case "$line" in
    ''|\#*)          continue ;;
    "FILE "*)        files+=("${line#FILE }") ;;
    "NORMALIZE "*)   norms+=("${line#NORMALIZE }") ;;
    "SORT "*)        sorts+=("${line#SORT }") ;;
  esac
done < "$contract"
[ ${#files[@]} -gt 0 ] || { echo "compare_outputs: contract lists no FILE" >&2; exit 2; }

# resolve the FILE globs to a de-duplicated relative-path list across BOTH dirs
declare -A SEEN; list=()
for g in "${files[@]}"; do
  for d in "$od" "$rd"; do
    for f in "$d"/$g; do
      [ -e "$f" ] || continue
      rel="${f#"$d"/}"
      [ -n "${SEEN[$rel]:-}" ] && continue; SEEN[$rel]=1; list+=("$rel")
    done
  done
done
[ ${#list[@]} -gt 0 ] || { echo "compare_outputs: no files matched the contract in either run" >&2; exit 2; }

sedargs=(); for nrm in "${norms[@]}"; do sedargs+=(-e "$nrm"); done
_want_sort() { local rel="$1" g; for g in "${sorts[@]}"; do case "$rel" in $g) return 0 ;; esac; done; return 1; }
_prep() { # <dir> <rel> -> normalized (+ optionally sorted) content; a sentinel if missing
  local f="$1/$2"
  [ -e "$f" ] || { printf '<<MISSING in %s>>\n' "$1"; return; }
  { if [ ${#sedargs[@]} -gt 0 ]; then sed -E "${sedargs[@]}" "$f"; else cat "$f"; fi; } \
    | { if _want_sort "$2"; then sort; else cat; fi; }
}

tmp="$(mktemp)"; fail=0
for rel in "${list[@]}"; do
  if diff <(_prep "$od" "$rel") <(_prep "$rd" "$rel") > "$tmp" 2>&1; then
    echo "MATCH $rel"
  else
    echo "DIFF  $rel"; sed 's/^/    /' "$tmp"; fail=1
  fi
done
rm -f "$tmp"
echo "RESULT: $([ "$fail" -eq 0 ] && echo MATCH || echo MISMATCH)"
[ "$fail" -eq 0 ]
