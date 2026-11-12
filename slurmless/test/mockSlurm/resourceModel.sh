#!/bin/bash
# ---------------------------------------------------------------------------
# resourceModel.sh  --  ground-truth resource curve for the Slurm-free harness.
#
# Given a model file and an input size, print the resources a job of that size
# "really" needs, optionally with reproducible stochastic noise:
#
#     resourceModel.sh <modelFile> <inputSize> [seed]
#         -> prints:  <memUsed_M> <timeUsed_min> <maxRSS_M>
#
# This is the synthetic stand-in for what memCpuMonitor.sh + sacct would have
# MEASURED on a real cluster. The mock sbatch calls it at submit time, caches
# the result per job id, and the mock memCpuMonitor/sacct replay it so that the
# UNMODIFIED cleanUp.sh writes a realistic row into jobRecord.txt. estimateResource.sh
# then fits a curve against those rows -- the whole learning loop, Slurm-free.
#
# Model file (KEY=VALUE, shell-sourceable). Example (linear mem in input KB):
#     MODEL_PROGRAM=grow          # optional; informational / match key
#     MODEL_REF=lin               # optional
#     MEM_FORM=linear             # linear | exp | poly | const
#     MEM_A=0.5                   # form params (see below)
#     MEM_B=200
#     TIME_FORM=linear
#     TIME_A=0.01
#     TIME_B=5
#     MAXRSS_FACTOR=1.0           # sacct MaxRSS = factor * memUsed (knob for
#                                 #   "sacct under/over-reports vs the cgroup peak")
#     NOISE_MEM=0.0               # fractional gaussian sd (0.05 = 5%)
#     NOISE_TIME=0.0
#     SEED=1234                   # base seed; final seed = SEED + <seed arg>
#
# Forms (x = inputSize in KB):
#     linear : y = A*x + B
#     exp    : y = B * exp(A*x)
#     poly   : y = A*x^2 + B*x + C      (MEM_C/TIME_C default 0)
#     const  : y = B
#
# Noise is multiplicative gaussian: y *= (1 + N(0, NOISE)), clamped to >= 1.
# It is seeded ONLY by (SEED + seed-arg), so a given (model, size, seed) always
# yields the same number -- reruns are byte-for-byte reproducible.
# ---------------------------------------------------------------------------
set -uo pipefail

modelFile="${1:?usage: resourceModel.sh <modelFile> <inputSize> [seed]}"
inputSize="${2:?inputSize required}"
seedArg="${3:-0}"

[ -f "$modelFile" ] || { echo "resourceModel.sh: no such model file: $modelFile" >&2; exit 2; }

# defaults, then override from the model file
MODEL_PROGRAM=""; MODEL_REF=""
MEM_FORM=const;  MEM_A=0; MEM_B=0; MEM_C=0
TIME_FORM=const; TIME_A=0; TIME_B=0; TIME_C=0
MAXRSS_FACTOR=1.0
NOISE_MEM=0; NOISE_TIME=0
SEED=0
# shellcheck disable=SC1090
source "$modelFile"

# non-numeric input size (e.g. "missingInputFile") -> no synthetic curve
case "$inputSize" in
  ''|*[!0-9]*) echo "0 0 0"; exit 0;;
esac

finalSeed=$(( SEED + seedArg ))

awk -v x="$inputSize" \
    -v mform="$MEM_FORM"  -v ma="$MEM_A"  -v mb="$MEM_B"  -v mc="$MEM_C" \
    -v tform="$TIME_FORM" -v ta="$TIME_A" -v tb="$TIME_B" -v tc="$TIME_C" \
    -v nmem="$NOISE_MEM"  -v ntime="$NOISE_TIME" \
    -v mfac="$MAXRSS_FACTOR" -v seed="$finalSeed" '
    function curve(form, a, b, c,    v) {
        if      (form == "linear") v = a*x + b
        else if (form == "exp")    v = b * exp(a*x)
        else if (form == "poly")   v = a*x*x + b*x + c
        else                       v = b            # const
        return v
    }
    # Self-contained LCG (Numerical Recipes constants, modulus 2^32). awk uses
    # doubles; 1664525*state stays < 2^53 so every step is an exact integer.
    # We roll our own because mawk ignores srand(constant) (non-reproducible),
    # and reproducible noise is the whole point of a curve-fitting test.
    function urand(    ) {
        _lcg = (1664525*_lcg + 1013904223) % 4294967296
        return _lcg / 4294967296
    }
    # Box-Muller gaussian, N(0,1)
    function gauss(    u1, u2) {
        u1 = urand(); if (u1 < 1e-12) u1 = 1e-12
        u2 = urand()
        return sqrt(-2.0*log(u1)) * cos(6.283185307179586 * u2)
    }
    BEGIN {
        _lcg = (seed % 4294967296 + 4294967296) % 4294967296
        # Decorrelate sequential seeds. An LCG K-th output is AFFINE in the seed,
        # so a fixed warm-up still leaves sequential seeds tracing a smooth (aliased)
        # curve. Breaking linearity needs a SEED-DEPENDENT number of warm-up steps.
        _warm = 16 + (seed % 97) + (seed % 13) * 7
        for (_w = 0; _w < _warm; _w++) urand()
        mem  = curve(mform, ma, mb, mc)
        # draw time noise AFTER mem noise so the two are distinct but reproducible
        if (nmem  > 0) mem  *= (1.0 + nmem  * gauss())
        tim  = curve(tform, ta, tb, tc)
        if (ntime > 0) tim  *= (1.0 + ntime * gauss())
        if (mem < 1) mem = 1
        if (tim < 1) tim = 1
        maxrss = mem * mfac
        if (maxrss < 1) maxrss = 1
        printf "%d %d %d\n", (mem+0.5), (tim+0.5), (maxrss+0.5)
    }'
