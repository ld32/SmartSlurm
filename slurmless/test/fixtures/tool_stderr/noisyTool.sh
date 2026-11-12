#!/bin/bash
echo "noisyTool: starting analysis of $1"
echo "noisyTool: processing records"
echo "WARNING: 3 records skipped" >&2
echo "noisyTool: done, wrote 0 usable records"
echo "ERROR: no usable records - aborting" >&2
exit 4
