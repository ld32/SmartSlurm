# SmarterSlurm

# How to use sbatchTop

cd ~

git clone git@github.com:ld32/smarterSlurm.git

export PATH=$HOME/smartSlurm/bin:$PATH

sbatchAndTop <sbatch option1> <sbatch option 2> <sbatch option 3> <...>

Such as: 

sbatchAndTop -p short -c 1 -t 2:0:0 --mem 2G --wrap "my_application para1 para2" # Here -p short is optional, because ssbatch chooses partition according to run time.

or: 

sbatchAndTop job.sh

# How to use ssbatch:

cd ~

git clone git@github.com:ld32/smarterSlurm.git

export PATH=$HOME/smartSlurm/bin:$PATH

ssbatch <sbatch option1> <sbatch option 2> <sbatch option 3> <...>

Such as: 

ssbatch -p short -c 1 -t 2:0:0 --mem 2G --wrap "my_application para1 para2" # Here -p short is optional, because ssbatch chooses partition according to run time.

or: 

ssbatch job.sh

Or if you would like to use this ssbatch command to replace the regular sbatch command (so that you don't have to modify your pipeline): 

sbatch() { $HOME/smartSlurm/bin/ssbatch $@; }  # define a bash function called sbatch 

export -f sbatch                             # enable it

type sbatch   # confirm the new sbatch function overwrites the regular sbatch command

Then you can run slurm jobs as usual. After you finish using ssbatch, run this command to disable it: 

unset sbatch

# ssbatch features:

1) Auto adjust partition according to run-time request if they does not match up

For example, this command:

sbatch -p medium -t 0-0:0:10 myjob.sh.  # notice mediume partition only allows job longer than 12 hours. 

becomes:

sbatch -p short -t 0-0:0:10 myjob.sh.   # medium partition is replaced by short partition

sbatch -t 0-0:0:10 --wrap hostname # notice there is no partition is given with the sbatch command

becomes:

sbatch -p short -t 0-0:0:10 --wrap hostname # notice short partition is chosen for this 10 minute job

This command will be not changed:

sbatch -p priority -t 0-0:0:10 myjob.sh # because prioirty partition allow any time less than 30 days, we keep to use priority partion. 

2) Auto check if slurm script exists

For example, using defaut sbatch, this command will fail witout eror if myjob.sh does not exist: 

sbatch -p priority -t 0-0:0:10 myjob.sh 

3) Auto create workding directory if not exist

For example, using default sbatch, this command fails if /home/ld/workdir does not exist

sbatch -D /home/ld/workdir -p priority -t 0-0:0:10 myjob.sh 

3) Auto create output and erro folders if not exist

For example, using default sbatch, this command fails if folder out or err does not exist

sbatch -p priority -t 0-0:0:10 -o out/out -e err/myjob.sh 
