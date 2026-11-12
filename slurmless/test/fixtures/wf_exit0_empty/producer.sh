#!/bin/bash
# Silent failure: a tool that EXITS 0 but never writes its promised output
# (out.txt). The scheduler sees success, so an afterok consumer is NOT orphaned —
# it runs and then fails for real on the missing input.
echo "producer: did nothing useful (no out.txt written)" >&2
exit 0
