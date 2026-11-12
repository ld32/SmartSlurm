#!/bin/bash
# ===========================================================================
# ma_dryrun.sh  --  repeatable Migration Assistant smoke test / worked example.
#
# Drives a full Scenario-2 conversion with the skill's own tools and the slurmless
# harness: a hand-rolled SLURM workflow (sbatch + afterok) -> rAP port, then a
# differential test proving the port reproduces the original's output. Exits 0 on
# MATCH, non-zero otherwise -- so it doubles as an install smoke test.
#
#   bash .../skills/migration-assistant/examples/ma_dryrun.sh
#
# Location-portable: it finds the harness relative to ITS OWN path, so it works
# the same run from a repo checkout or from an installed plugin cache
# (${CLAUDE_PLUGIN_ROOT}/skills/migration-assistant/examples/ma_dryrun.sh).
#
# Env:
#   SMARTSLURM_ROOT   SmartSlurm checkout to test against (else auto-resolved:
#                     vendored parent -> sibling ../SmartSlurm; set it if neither).
#   SS_SCRATCH        scratch root (else a unique mktemp dir, cleaned on success).
#   KEEP_SCRATCH=1    keep the scratch dir even on success (for inspection).
# ===========================================================================
set -uo pipefail

# --- $SL = the slurmless root, resolved from this script's own location --------
SL="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
TOOLS="$SL/skills/migration-assistant/tools"
[ -f "$SL/test/lib/harness.sh" ] || { echo "ma_dryrun: cannot find harness at \$SL=$SL" >&2; exit 2; }

# --- Resolve a SmartSlurm checkout (mirror bootstrap.sh's chain) ---------------
if [ -z "${SMARTSLURM_ROOT:-}" ]; then
  d="$SL"
  while [ "$d" != "/" ]; do
    [ -e "$d/bin/runAsPipeline" ] && { SMARTSLURM_ROOT="$d"; break; }
    d="$(dirname "$d")"
  done
  : "${SMARTSLURM_ROOT:=$(dirname "$SL")/SmartSlurm}"
fi
if [ ! -e "$SMARTSLURM_ROOT/bin/runAsPipeline" ]; then
  echo "ma_dryrun: no SmartSlurm at SMARTSLURM_ROOT=$SMARTSLURM_ROOT" >&2
  echo "           set SMARTSLURM_ROOT to a SmartSlurm checkout and retry." >&2
  exit 2
fi
export SMARTSLURM_ROOT

OWN_SCRATCH=0
if [ -z "${SS_SCRATCH:-}" ]; then SS_SCRATCH="$(mktemp -d /tmp/ma-dryrun-XXXXXX)"; OWN_SCRATCH=1; fi
export SS_SCRATCH
echo "### SL=$SL"
echo "### SMARTSLURM_ROOT=$SMARTSLURM_ROOT"
echo "### SS_SCRATCH=$SS_SCRATCH"

# --- The USER's hand-rolled SLURM workflow (Scenario 2 input) ------------------
STAGE="$SS_SCRATCH/userpipe"; mkdir -p "$STAGE"
cat > "$STAGE/make.sh"  <<'EOF'
#!/bin/bash
echo "hello slurmless world" > greeting.txt
EOF
cat > "$STAGE/shout.sh" <<'EOF'
#!/bin/bash
tr 'a-z' 'A-Z' < greeting.txt > result.txt
EOF
cat > "$STAGE/submit.sh" <<'EOF'
#!/bin/bash
# produce greeting, then (only on success) shout it -- classic afterok chain
jid=$(sbatch make.sh | awk '{print $NF}')
sbatch --dependency=afterok:$jid shout.sh
EOF
chmod +x "$STAGE"/*.sh
echo; echo "===== INPUT: user's hand-rolled SLURM workflow (submit.sh) ====="; cat "$STAGE/submit.sh"

# --- Harness --------------------------------------------------------------------
# shellcheck disable=SC1090
source "$SL/test/lib/harness.sh" && ss_activate >/dev/null 2>&1

# ===== ARM 1: recover the DAG (register-don't-run), then run the original =======
ss_new_workspace "$SS_SCRATCH/orig" >/dev/null
cp "$STAGE"/*.sh .
ss_wf_activate
bash ./submit.sh >/dev/null 2>&1                 # registers jobs; nothing runs yet
echo; echo "===== STEP 1: recovered DAG (extract_dag.sh) ====="
bash "$TOOLS/extract_dag.sh" | tee dag.tsv
ss_wf_drain >/dev/null 2>&1                       # execute the real DAG
REF="$(cat result.txt 2>/dev/null)"
echo; echo "### ORIGINAL output (REF): [$REF]"

# ===== STEP 2: author the logical plan, emit the rAP port ======================
cat > "$SS_SCRATCH/plan.txt" <<'EOF'
# Two linear steps, no fan-out: make -> (afterok) -> shout
SBATCH  sbatch -p short -c 1 --mem 100M -t 5:00
STEP make  dep=-    fanout=no in=none
    echo "hello slurmless world" > greeting.txt
STEP shout dep=make fanout=no in=none
    tr 'a-z' 'A-Z' < greeting.txt > result.txt
EOF
echo; echo "===== STEP 2: emitted rAP script (emit_rap.sh) ====="
bash "$TOOLS/emit_rap.sh" "$SS_SCRATCH/plan.txt" | tee "$SS_SCRATCH/pipeline.sh"

# ===== ARM 2: run the rAP port in its own workspace ============================
ss_new_workspace "$SS_SCRATCH/rap" >/dev/null
cp "$SS_SCRATCH/pipeline.sh" .
echo; echo "===== STEP 3: run the rAP port (ss_run_pipeline) ====="
ss_run_pipeline "pipeline.sh" "sbatch -p short -c 1 --mem 100M -t 5:00"
PORT="$(cat result.txt 2>/dev/null)"
echo "### PORT output: [$PORT]  (runAsPipeline log: $LAST_RUN_LOG)"

# ===== STEP 4: differential test under an identity contract ====================
printf 'FILE result.txt\n' > "$SS_SCRATCH/contract.txt"
echo; echo "===== STEP 4: differential (compare_outputs.sh) ====="
if bash "$TOOLS/compare_outputs.sh" "$SS_SCRATCH/orig" "$SS_SCRATCH/rap" "$SS_SCRATCH/contract.txt"; then
  verdict=0
else
  verdict=1
fi

echo
if [ "$verdict" -eq 0 ] && [ -n "$REF" ] && [ "$REF" = "$PORT" ]; then
  echo "##### ma_dryrun: PASS (port reproduces the original: [$PORT]) #####"
  [ "$OWN_SCRATCH" -eq 1 ] && [ -z "${KEEP_SCRATCH:-}" ] && rm -rf "$SS_SCRATCH"
  exit 0
else
  echo "##### ma_dryrun: FAIL (REF=[$REF] PORT=[$PORT] compare=$verdict) #####" >&2
  echo "scratch kept for inspection: $SS_SCRATCH" >&2
  exit 1
fi
