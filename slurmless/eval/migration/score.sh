#!/bin/sh
# score.sh <variantDir> : verdict from the variant's REF (original behavior), PORT
# (the agent's rAP conversion output), and the golden quote.
here="$(cd "$(dirname "$0")" && pwd)"
v="$1"; [ -d "$v" ] || { echo "no such variant: $v" >&2; exit 2; }
golden=$(cat "$here/golden/quote.txt")
ref=$(cat "$v/ref/quote.txt" 2>/dev/null)
port=$(cat "$v/port/quote.txt" 2>/dev/null)
defect=$(cat "$v/.defect" 2>/dev/null)
echo "=== $(basename "$v") ==="
echo "  golden: $golden"
echo "  REF   : ${ref:-<missing>}"
echo "  PORT  : ${port:-<missing>}"
[ -n "$defect" ] && echo "  planted defect: $defect"

if [ -z "$ref" ]; then echo "VERDICT: SKIP (no REF — run run_original.sh first)"; exit 2; fi
if [ -z "$port" ]; then echo "VERDICT: FAIL (no PORT — conversion produced no runnable pipeline)"; exit 1; fi
if [ "$port" = "$ref" ]; then
  if [ "$ref" = "$golden" ]; then echo "VERDICT: PASS (faithful + correct: port == original == golden)"
  else echo "VERDICT: PASS (faithful: planted bug preserved; port == original, != golden)"; fi
  exit 0
fi
if [ "$port" = "$golden" ] && [ "$ref" != "$golden" ]; then
  echo "VERDICT: FAIL (infidelity: agent silently 'fixed' the planted bug — port == golden but original != golden)"
else
  echo "VERDICT: FAIL (conversion defect: port != original)"
fi
exit 1
