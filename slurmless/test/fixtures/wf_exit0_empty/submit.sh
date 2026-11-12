#!/bin/bash
# producer "succeeds" (exit 0) but writes nothing; consumer runs afterok of it.
set -u
p=$(sbatch -J producer producer.sh | awk '{print $4}')
sbatch -J consumer --dependency=afterok:"$p" consumer.sh >/dev/null
