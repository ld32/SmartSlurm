# SmarterSlurm

ssbath was designed to run https://github.com/ENCODE-DCC/atac-seq-pipeline, so that users don't have to modify the original workflow and sbatch can automatially modify the partitions according user's local cluster partition settings. 

# How to use ssbatch:

cd ~    
git clone git://github.com/ld32/smarterSlurm.git  
export PATH=$HOME/smartSlurm/bin:$PATH  
ssbatch <sbatch option1> <sbatch option 2> <sbatch option 3> <...>

Such as:     
ssbatch -p short -c 1 -t 2:0:0 --mem 2G --wrap "my_application para1 para2" # Here -p short is optional, because ssbatch chooses partition according to run time.   
or:     
ssbatch job.sh

Or if you would like to use this ssbatch command to replace the regular sbatch command (so that you don't have to modify your pipeline, such as https://github.com/ENCODE-DCC/atac-seq-pipeline):    
sbatch() { $HOME/smartSlurm/bin/ssbatch "$@"; }  # define a bash function called sbatch   
export -f sbatch                  # enable the function it    
type sbatch                       # confirm the new sbatch function overwrites the regular sbatch command

Then you can run slurm jobs as usual. After you finish using ssbatch, run this command to disable it:    
unset sbatch

# ssbatch features:

1) Auto adjust partition according to run-time request if they does not match up

For example, this command:  
ssbatch -p medium -t 0-0:0:10 myjob.sh.  # notice mediume partition only allows job longer than 12 hours.    
becomes:    
sbatch -p short -t 0-0:0:10 myjob.sh.   # medium partition is replaced by short partition

ssbatch -t 0-0:0:10 --wrap hostname # notice there is no partition is given with the sbatch command  
becomes:    
sbatch -p short -t 0-0:0:10 --wrap hostname # notice short partition is chosen for this 10 minute job

This command will be not change the partition:   
ssbatch -p priority -t 0-0:0:10 myjob.sh # because prioirty partition allow any time less than 30 days, we keep to use priority partion. 

2) Auto check if slurm script exists    
For example, using defaut sbatch, this command will fail witout eror if myjob.sh does not exist:    
ssbatch -p priority -t 0-0:0:10 myjob.sh 

3) Auto create workding directory if not exist  
For example, using default sbatch, this command fails if /home/ld/workdir does not exist    
ssbatch -D /home/ld/workdir -p priority -t 0-0:0:10 myjob.sh 

3) Auto create output and erro folders if not exist     
For example, using default sbatch, this command fails if folder out or err does not exist       
ssbatch -p priority -t 0-0:0:10 -o out/out -e err/myjob.sh 

# How does it works

config/partitions.txt contains partion time limit and bash function adjustPartition to adjust partion for sbatch jobs: 

\# Genernal partions, ordered by maximum allowed run-time in hours 
partition1Name=short   
partition2Name=medium  
partition3Name=long        
partition4Name=      
partition5Name=     
partition1TimeLimit=12  # run-time > 0 hours and <= 12 hours    
partition2TimeLimit=120 # run-time > 12 hours and <= 5 days     
partition3TimeLimit=720 # run-time > 5 days and <= 30 days  
partition4TimeLimit=    
partition5TimeLimit=    

\# Special pertitions with special restrictions
partition6Name=priority    # only allow two job running at the same time        
partition7Name=highmem     # run-time <= 30 days, special partision     
partition8Name=interactive      
partition9Name=mpi      
partition10Name=        

partition6TimeLimit=720 # run-time <= 30 days   
partition7TimeLimit=720 # run-time <= 30 days   
partition8TimeLimit=12  # run-time <= 12 hours      
partition9TimeLimit=720 # run-time <= 30 days       
partition10TimeLimit=       

\#function 
adjustPartition() {         
    ... # please open the file to see the content         
}       
export -f adjustPartition    

# How to use sbatchAndTop

cd ~    
git clone git@github.com:ld32/smarterSlurm.git  
export PATH=$HOME/smartSlurm/bin:$PATH    
sbatchAndTop <sbatch option1> <sbatch option 2> <sbatch option 3> <...> 

Such as:    
sbatchAndTop -p short -c 1 -t 2:0:0 --mem 2G --wrap "my_application para1 para2" # Here -p short is optional, because ssbatch chooses partition according to run time.  

or:     
sbatchAndTop job.sh 

# sbatchAndTop features:

1) Submit slurm job using ssbatch (scrool down to see ssbatch features) and run scontrol top on the job
