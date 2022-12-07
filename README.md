=============================
# SmartSlurm
- [ssbatch](#ssbatch)
    - [ssbatch features](#ssbatch-features)
    - [Quick Start](#quick-start)
    - [How to Use ssbatch in Details]($how-to-use-ssbatch-in-details) 
    - [How does It Works]($how-does-it-works)
- [sbatchAndTop](#sbatchAndTop)
# ssbatch
ssbath was originally designed to run https://github.com/ENCODE-DCC/atac-seq-pipeline, so that users don't have to modify the original workflow and sbatch can automatially modify the partitions according user's local cluster partition settings. The script was later improved to have more features.
## ssbatch features:
1) Auto adjust memory and run-time according to the input size based on statistics from earlier jobs
2) Auto choose partition according to run-time request
3) Auto re-run failed OOM (out of memory) and OOT (out of run-time) jobs
4) Get good emails: by default Slurm emails only have a subject. ssbatch attaches the content of the sbatch script, the output and error log to email
## Quick Start
``` bash
# Install and setup
cd $HOME 
git clone git://github.com/ld32/smarterSlurm.git  
export PATH=$HOME/smartSlurm/bin:$PATH  
sbatch() { $HOME/smartSlurm/bin/ssbatch "$@"; }; export -f sbatch                                 

# Create some text files for testing
createBigTextFiles.sh

# Notice, you don't have to run this section, because I have run it and save the statistics in $HOME/smartSlurm
# Run 5 jobs to get memory and run-time statistics
for i in {1..5}; do
    export SSBATCH_I=bigText$i.txt # This is to tell ssbatch the input file to calculate input file size
    sbatch -t 0:30:0 --mem 2000M --wrap="~/smartSlurm/bin/useSomeMemTimeAccordingInputSize.sh bigText$i.txt"
done

# Auto adjust memory and run-time according input file size
export SSBATCH_I=bigText1.txt,bigText2.txt # This is to tell ssbatch the input file to calculate input file size 
sbatch -t 0:30:0 --mem 2000M --wrap="~/smartSlurm/bin/useSomeMemTimeWithInput.sh bigText1.txt bigText$2.txt 3"

# Notice, you don't have to run this section, because I have run it and save the statistics in $HOME/smartSlurm
# Run 3 jobs to get memory and run-time statistics
unset SSBATCH_I
for i in {1..3}; do
    sbatch -t 0:30:0 --mem 2000M --wrap="~/smartSlurm/bin/useSomeMemTimeNoInput.sh $i $i"
done

# Auto adjust memory and run-time so that 90% jobs can finish successfully
sbatch -t 0:30:0 --mem 2000M --wrap="~/smartSlurm/bin/useSomeMemTimeNoInput.sh bigText1.txt 1 2"

# After you finish using ssbatch, run this command to disable it:    
unset sbatch
unset SSBATCH_I
```

## How to use ssbatch in details
``` bash
cd $HOME 
git clone git://github.com/ld32/smarterSlurm.git  
export PATH=$HOME/smartSlurm/bin:$PATH  
ssbatch <sbatch option1> <sbatch option 2> <sbatch option 3> <...>

# Such as:     
ssbatch -p short -c 1 -t 2:0:0 --mem 2G --wrap "my_application para1 para2"
# Notice above, '-p short' is optional because ssbatch can automatically choose partition according to run time.   

#or:     
ssbatch job.sh

# Or if you would like to use this ssbatch command to replace the regular sbatch command (so that you 
# don't have to modify your pipeline, such as https://github.com/ENCODE-DCC/atac-seq-pipeline). 
# Please run these two commands before submitting sbatch jobs:
sbatch() { $HOME/smartSlurm/bin/ssbatch "$@"; }  # define a bash function called sbatch   
export -f sbatch                                 # enable the function it    

#Then you can run sbatch jobs as usual. After you finish using ssbatch, run this command to disable it:    
unset sbatch
```



## How does it works
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



# sbatchAndTop
## How to use sbatchAndTop
```
cd ~    
git clone git@github.com:ld32/smarterSlurm.git  
export PATH=$HOME/smartSlurm/bin:$PATH    
sbatchAndTop <sbatch option1> <sbatch option 2> <sbatch option 3> <...> 

# Such as:    
sbatchAndTop -p short -c 1 -t 2:0:0 --mem 2G --wrap "my_application para1 para2" 
# Here -p short is optional, because ssbatch chooses partition according to run time.  

# or:     
sbatchAndTop job.sh 

## sbatchAndTop features:
1) Submit slurm job using ssbatch (scrool up to see ssbatch features) and run scontrol top on the job
