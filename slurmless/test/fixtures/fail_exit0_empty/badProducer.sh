#!/bin/bash
# Cannot read its input; writes nothing to stdout; STILL exits 0 (the bug class).
echo "warning: input not found, continuing" >&2
exit 0
