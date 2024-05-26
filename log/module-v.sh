#!/bin/bash
date

trap "{ cleanUp.sh module-v module none 0 1 2048 10 2048 10 short  \"rccg\" \"none\" 10 5; }" EXIT
memCpuMonitor.sh module-v 2048 10 2048 1 &

unset SLURM_CPU_BIND
srun -n 1 -A rccg sh -e -o pipefail -c 'module -v; ' && touch log/module-v.success 

#Command used to submit the job: /usr/bin/sbatch --mail-type=FAIL --requeue --parsable -p $myPartition -c 1 --mem $myMem -t $myTime --open-mode=append -o log/module-v.out -e log/module-v.out -J module-v   -A rccg  log/module-v.sh

#myMem=2048 myTime=0-00:10:00
