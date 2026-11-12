#!/bin/bash
# toolstep needs an uninstalled tool; consumer runs afterok of it.
set -u
A=$(sbatch -J toolstep toolstep.sh | awk '{print $4}')
sbatch -J consumer --dependency=afterok:"$A" job.sh consumer 0 >/dev/null
# a --wrap step calling a missing tool directly (exercises SS_WF_STUB_MISSING)
sbatch -J wrapstep --wrap "nosuchtool_abc hello" >/dev/null
