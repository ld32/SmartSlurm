#!/bin/bash
# =============================================================================
# sizedLoad.sh <input_file>
# A test "tool" that consumes MEMORY and WALL-TIME in proportion to the number
# of rows in its input, deterministically and with tunable knobs. Use it to
# exercise SmartSlurm's jobRecord fitting + resource auto-increment: run it over
# several differently-sized inputs so rAP can fit a curve, and set FAIL_PROB>0
# to make it overshoot the request and trigger OOM/OOT auto-rerun.
#
# Env knobs (all optional):
#   MEM_KB_PER_ROW  memory held per input row, KB      (default 20)
#   SEC_PER_ROW     seconds slept per input row        (default 0.02)
#   NOISE_PCT       +/- jitter on the totals, percent  (default 5)
#   FAIL_PROB       chance [0..1] of a deliberate overshoot (default 0 = never)
#   FAIL_MODE       oom | oot | both                   (default both)
#   OVERSHOOT_X     multiplier applied on overshoot    (default 4)
# =============================================================================
set -u
in=${1:?usage: sizedLoad.sh <input_file>}
rows=$(wc -l < "$in" 2>/dev/null); rows=${rows:-0}

: "${MEM_KB_PER_ROW:=20}"; : "${SEC_PER_ROW:=0.02}"; : "${NOISE_PCT:=5}"
: "${FAIL_PROB:=0}"; : "${FAIL_MODE:=both}"; : "${OVERSHOOT_X:=4}"

# base targets ~ rows, with +/- NOISE_PCT jitter
vals=$(awk -v r="$rows" -v mk="$MEM_KB_PER_ROW" -v sp="$SEC_PER_ROW" -v n="$NOISE_PCT" 'BEGIN{
  srand(); jm=1+(rand()*2-1)*n/100; jt=1+(rand()*2-1)*n/100;
  printf "%d %.2f", (r*mk*jm)+1, (r*sp*jt) }')
mem_kb=${vals% *}; sec=${vals#* }

# optional deterministic overshoot (for the OOM/OOT rerun test)
if awk -v p="$FAIL_PROB" 'BEGIN{srand(); exit !(rand()<p)}'; then
  case "$FAIL_MODE" in
    oom)  mem_kb=$(awk -v m="$mem_kb" -v x="$OVERSHOOT_X" 'BEGIN{printf "%d", m*x}') ;;
    oot)  sec=$(awk   -v s="$sec"    -v x="$OVERSHOOT_X" 'BEGIN{printf "%.2f", s*x}') ;;
    *)    mem_kb=$(awk -v m="$mem_kb" -v x="$OVERSHOOT_X" 'BEGIN{printf "%d", m*x}')
          sec=$(awk   -v s="$sec"    -v x="$OVERSHOOT_X" 'BEGIN{printf "%.2f", s*x}') ;;
  esac
  echo "sizedLoad: OVERSHOOT fired (mode=$FAIL_MODE x$OVERSHOOT_X)" >&2
fi

echo "sizedLoad: rows=$rows -> hold ~${mem_kb} KB for ~${sec}s" >&2
# hold ~mem_kb KB in an awk array for ~sec seconds (both scale with rows)
awk -v kb="$mem_kb" -v sec="$sec" 'BEGIN{
  n=int(kb);
  for(i=0;i<n;i++) a[i]=sprintf("%1023d",i);  # distinct ~1KB entries -> RSS ~ n KB
  system("sleep " sec) }'                      # hold it for the duration
echo "sizedLoad: done" >&2
