#!/bin/bash
# A step that depends on a tool which may not be installed in the sandbox. Under
# the REAL tier it fails (set -e + missing tool); the tiered harness can STUB it
# so the pipeline's wiring still validates.
set -e
nosuchtool_xyz --do-science
echo "ran toolstep" > ran.toolstep
