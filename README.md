# betterSlurm

slurm wrappers to work-around the limitations or bugs

1 Auto choose partition according run time request. Otherwise, if give wrong partition name, sbatch give error.

2 Auto create output folder if not exist. Otherwise if output folder does not exist, job dies silently.
