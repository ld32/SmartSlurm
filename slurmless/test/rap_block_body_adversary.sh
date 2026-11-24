#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# rap_block_body_adversary.sh -- regression guards for two SILENT bugs in the
# T2 verbatim #@end body / submit-scope baking path of bin/runAsPipeline
# (commits 915cb0e, 8dd423a on feat/rap-block-body).
#
#   SMARTSLURM_ROOT=/path/to/smartslurm test/rap_block_body_adversary.sh
#
# Both are REGRESSIONS vs ld32/master (master produced the correct result; the
# T2 baking produces a silently wrong one, exit 0). Each guard asserts the
# CORRECT result, so it is RED on the buggy tool and GREEN once fixed.
#
#   G1  submit-scope ARRAY baked as a scalar (only element 0 survives)
#   G2  positional-param bake fabricates `set -- ''` when the pipeline has no
#       args, so $# becomes 1 instead of 0
# ---------------------------------------------------------------------------
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/lib/harness.sh"
ss_activate >/dev/null 2>&1

fails=0
pass() { echo "    PASS - $1"; }
fail() { echo "    FAIL - $1"; [ -n "${2:-}" ] && echo "           $2"; fails=$((fails+1)); }
xfail() { echo "    XFAIL (known limit, pre-existing, deferred) - $1"; [ -n "${2:-}" ] && echo "           $2"; }

_convert_run() {  # $1=workspace tag ; stdin = full pipeline.sh
  ss_new_workspace "$WORK_ADV/$1" >/dev/null 2>&1 || ss_new_workspace "/tmp/ss-adv/$1" >/dev/null 2>&1
  cat > pipeline.sh; chmod +x pipeline.sh
  ss_run_pipeline "./pipeline.sh" "sbatch -p short -t 5" >/dev/null 2>&1 || true
}
: "${WORK_ADV:=/tmp/ss-adv}"

echo "== [rap_block_body_adversary] T2 baking silent-correctness guards =="

# --- G1: submit-scope array must not be silently truncated to element 0 -------
_convert_run g1_array <<'PIPE'
#!/bin/sh
samples=(alpha beta gamma)
#@1,0,t,,,sbatch -p short -t 5
echo "count=${#samples[@]} all=${samples[*]}" > arr.out
printf 'ran\n' > ran.txt
#@end
PIPE
got="$(cat arr.out 2>/dev/null)"
if [ "$got" = "count=3 all=alpha beta gamma" ]; then
  pass "G1: submit-scope array resolves fully on the node"
else
  fail "G1: submit-scope array baked as scalar (silent data loss)" \
       "want 'count=3 all=alpha beta gamma'  got '$got'  (cmd.sh bakes: $(grep -m1 '^samples=' smartSlurmLog/jobs/1.0.t/cmd.sh 2>/dev/null))"
fi

# --- G2: no spurious positional param when the pipeline has no args ----------
_convert_run g2_pos <<'PIPE'
#!/bin/sh
#@1,0,t,,,sbatch -p short -t 5
echo "count=$#" > pos.out
printf 'ran\n' > ran.txt
#@end
PIPE
got="$(cat pos.out 2>/dev/null)"
if [ "$got" = "count=0" ]; then
  pass "G2: no-arg pipeline leaves \$# at 0"
else
  fail "G2: 'set -- ''' fabricates an empty positional param (\$# wrong)" \
       "want 'count=0'  got '$got'  (cmd.sh: $(grep -m1 '^set --' smartSlurmLog/jobs/1.0.t/cmd.sh 2>/dev/null))"
fi

# --- G3: ASSOCIATIVE array must bake fully (re-review residual of Finding 1) ---
# The dense-indexed-array fix uses declare -p, but the bake is gated by
# `[ -z "${!v+x}" ] && continue`, whose indirect set-test only sees index/key 0.
# An associative array has no [0] key, so it is SILENTLY SKIPPED -> empty on node.
_convert_run g3_assoc <<'PIPE'
#!/bin/sh
declare -A m=([red]=apple [green]=lime)
#@1,0,t,,,sbatch -p short -t 5
echo "red=${m[red]} n=${#m[@]}" > m.out
printf 'ran\n' > ran.txt
#@end
PIPE
got="$(cat m.out 2>/dev/null)"
if [ "$got" = "red=apple n=2" ]; then
  pass "G3: submit-scope associative array resolves fully on the node"
else
  fail "G3: associative array silently skipped by \${!v+x} set-test (empty on node)" \
       "want 'red=apple n=2'  got '$got'  (cmd.sh bake lines: $(grep -c 'declare -A m' smartSlurmLog/jobs/1.0.t/cmd.sh 2>/dev/null))"
fi

# --- G4: SPARSE indexed array (no index 0) must bake fully (same residual) -----
_convert_run g4_sparse <<'PIPE'
#!/bin/sh
declare -a sp=([3]=three [7]=seven)
#@1,0,t,,,sbatch -p short -t 5
echo "n=${#sp[@]} e3=${sp[3]}" > sp.out
printf 'ran\n' > ran.txt
#@end
PIPE
got="$(cat sp.out 2>/dev/null)"
if [ "$got" = "n=2 e3=three" ]; then
  pass "G4: submit-scope sparse indexed array resolves fully on the node"
else
  fail "G4: sparse indexed array (no [0]) silently skipped (empty on node)" \
       "want 'n=2 e3=three'  got '$got'  (cmd.sh bake lines: $(grep -c 'declare -a sp' smartSlurmLog/jobs/1.0.t/cmd.sh 2>/dev/null))"
fi

# --- G5: submit-scope var name must not collide with a driver-internal global ---
# The driver's progress helpers use UNSCOPED globals (ssHeader: `d=$(printf ...)`
# with no `local d`; ssStart: `x=`, `grp=`). A submit-scope var of the same name is
# clobbered BEFORE baking, so the body receives the driver's scratch value (a string
# of spaces for `d`) instead of the user's. Silent (exit 0), pre-existing on master
# too. Fix: `local` those helper vars (or namespace all driver internals), and/or add
# them to the bake denylist. NOTE: a denylist-only fix would make `d` resolve on the
# NODE (likely empty), which is still wrong -- the real fix is to stop the driver from
# owning the name `d` at all.
_convert_run g5_dcollide <<'PIPE'
#!/bin/sh
d=myvalue
#@1,0,t,,,sbatch -p short -t 5
echo "d=[$d]" > d.out
printf 'ran\n' > ran.txt
#@end
PIPE
got="$(cat d.out 2>/dev/null)"
if [ "$got" = "d=[myvalue]" ]; then
  pass "G5: submit-scope var 'd' survives (no driver-internal collision)"
else
  fail "G5: submit-scope var 'd' clobbered by driver helper global (silent)" \
       "want 'd=[myvalue]'  got '$got'  (ssHeader's unscoped 'd=' pollutes it)"
fi

# --- G6: indirect ${!ptr} to a submit-scope var resolves on the node (FIXED) -------
# Was a silent REGRESSION vs ld32/master: master's unquoted here-doc expanded ${!ptr}
# at submit time, but the T2 bake baked the POINTER 'ptr' and never the variable it
# names (target), so ${!ptr} was empty on the node while the job exited 0. Fixed --
# the driver now resolves each true-indirect pointer and bakes its target too. This
# guard protects the fix against re-regression.
_convert_run g6_indirect <<'PIPE'
#!/bin/sh
target=secretvalue
ptr=target
#@1,0,t,,,sbatch -p short -t 5
echo "indirect=[${!ptr}]" > ind.out
printf 'ran\n' > ran.txt
#@end
PIPE
got="$(cat ind.out 2>/dev/null)"
if [ "$got" = "indirect=[secretvalue]" ]; then
  pass "G6: indirect \${!ptr} to a submit-scope var resolves on the node"
else
  fail "G6: indirect \${!ptr} silently empty on node (REGRESSION vs master, exit 0)" \
       "want 'indirect=[secretvalue]'  got '$got'  (bake baked 'ptr' but not its target)"
fi

# --- G7: submit-scope var referenced ONLY inside $(( )) -> SILENT WRONG ------------
# NOT a regression (master was equally wrong: it escaped '$(' -> arithmetic ran on the
# node). Pinned here because it is a latent silent bug the T2 bake does not close and
# the fuzz corpus does not cover (05_arrays_arith uses only literal operands). The bake
# regex requires a '$' sigil; a bare arithmetic operand (count in $((count+1))) is not
# detected, so count is unset on the node and $((count+1)) evaluates to 1, exit 0.
_convert_run g7_arith <<'PIPE'
#!/bin/sh
count=5
#@1,0,t,,,sbatch -p short -t 5
echo "arith=$((count+1))" > ar.out
printf 'ran\n' > ran.txt
#@end
PIPE
got="$(cat ar.out 2>/dev/null)"
if [ "$got" = "arith=6" ]; then
  pass "G7: (now fixed!) submit-scope var used only in \$(( )) resolves on the node -- promote this to a hard assert"
else
  xfail "G7: arithmetic-only submit-scope var not baked -> \$((count+1)) == 1 (F2; pre-existing, NOT a regression; deferred by decision)" \
        "got '$got'  (bake regex needs a \$ sigil; a bare \$((count)) operand is missed). Fix = bake bare identifiers inside \$((...))."
fi

echo
if [ "$fails" -eq 0 ]; then echo "rap_block_body_adversary: ALL GREEN"; exit 0
else echo "rap_block_body_adversary: $fails FAILING (silent T2 baking bug present)"; exit 1; fi
