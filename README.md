# betterSlurm


# features:

1) Auto adjust partition according to run-time request if they does not match up

For example, this command:

sbatch -p medium -t 0-0:0:10 myjob.sh

becomes:

sbatch -p short -t 0-0:0:10 myjob.sh

This command will be not changed:

sbatch -p priority -t 0-0:0:10 myjob.sh

2) Auto check if slurm script exists

For example, this command will fail witout eror if myjob.sh does not exist: 

sbatch -p priority -t 0-0:0:10 myjob.sh 

3) Auto create workding directory if not exist

For example, this command fails if /home/ld/workdir does not exist

sbatch -D /home/ld/workdir -p priority -t 0-0:0:10 myjob.sh 

3) Auto create output and erro folders if not exist

For example, this command fails if folder out or err does not exist

sbatch -p priority -t 0-0:0:10 -o out/out -e err/myjob.sh 
