#!/bin/bash
# ---------------------------------------------------------------------------
# parse_rap_dag.sh <rAP driver>  ->  the step DAG of an annotated rAP script.
#
# The read-back companion to emit_rap.sh: parse a driver's runAsPipeline markers
#   #@<num>,<dependsOn(dot-joined)>,<name>,<...>
# and print one line per step with its dependency step-numbers resolved to NAMES:
#
#   <num>  <name>        dep=<depName[.depName...]|->
#
# A dependency-free step prints `dep=-`. Use it to verify an emitted port's
# structure against the plan reviewed at checkpoint 1, or to diff two drivers'
# wiring. It reads the driver text only (no harness, no SmartSlurm needed).
# ---------------------------------------------------------------------------
set -uo pipefail
f="${1:?usage: parse_rap_dag.sh <rAP-driver>}"
[ -r "$f" ] || { echo "parse_rap_dag: cannot read '$f'" >&2; exit 2; }

awk -F',' '
  /^[[:space:]]*#@[0-9]/ {
    line=$0; sub(/^[[:space:]]*#@/,"",line); split(line,a,",");
    num=a[1]; deps=a[2]; name=a[3];
    numName[num]=name; order[++k]=num; dep[num]=deps;
  }
  END {
    if (k==0) { print "parse_rap_dag: no #@ step markers found" > "/dev/stderr"; exit 3 }
    for (i=1;i<=k;i++){ num=order[i]; d=dep[num]; out="-";
      if (d!="" && d!="0"){ m=split(d,b,"."); out="";
        for (j=1;j<=m;j++){ nm=(b[j] in numName)?numName[b[j]]:("#"b[j]);
          out=(out==""?nm:out"."nm); } }
      printf "%2s  %-14s dep=%s\n", num, numName[num], out;
    }
  }' "$f"
