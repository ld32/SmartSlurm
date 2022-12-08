=============================
# SmartSlurm
- [ssbatch](#ssbatch)
    - [ssbatch features](#ssbatch-features)
    - [Quick Start](#quick-start)
    - [How to Use ssbatch]($how-to-use-ssbatch) 
    - [How does It Works]($how-does-it-works)
- [sbatchAndTop](#sbatchAndTop)
# ssbatch
ssbath was originally designed to run https://github.com/ENCODE-DCC/atac-seq-pipeline, so that users don't have to modify the original workflow and sbatch can automatially modify the partitions according user's local cluster partition settings. The script was later improved to have more features.
## ssbatch features:
1) Auto adjust memory and run-time according to statistics from earlier jobs
2) Auto choose partition according to run-time request
3) Auto re-run failed OOM (out of memory) and OOT (out of run-time) jobs
4) Get good emails: by default Slurm emails only have a subject. ssbatch attaches the content of the sbatch script, the output and error log to email
## Quick Start
``` bash
# Download
cd $HOME 
git clone git://github.com/ld32/smartSlurm.git  

# Setup path
export PATH=$HOME/smartSlurm/bin:$PATH  

# Set up a function so that ssbatch is called when running sbatch
sbatch() { $HOME/smartSlurm/bin/ssbatch "$@"; }; export -f sbatch                                 

# Create some text files for testing
createBigTextFiles.sh

# Run 3 jobs to get memory and run-time statistics for useSomeMemTimeNoInput.sh
# Notice 1: you don't have to run this section, because I have run it and save the statistics in $HOME/smartSlurm
# Notice 2: Slurm will submit three jobs to short partition, each reserves 2G memory and 2 hour run-time 
export SSBATCH_S=useSomeMemTimeNoInput.sh # This is optional because the software name is the same as the script
for i in {1..3}; do
    sbatch --mem 2G -t 2:0:0 --wrap="useSomeMemTimeNoInput.sh $i"
done

# After the 3 jobs finish, auto adjust memory and run-time so that 90% jobs can finish successfully
# Notice: Slurm will submit five jobs to short partition, and reserved 19M memory and 7 minute run-time 
sbatch --mem 2G -t 2:0:0 --mem 2G --wrap="useSomeMemTimeNoInput.sh bigText1.txt"

# After you finish using ssbatch, run this command to disable it:    
unset sbatch
unset SSBATCH_S
```

## How to use ssbatch
``` bash
# Download
cd $HOME 
git clone git://github.com/ld32/smartSlurm.git  

# Setup path
export PATH=$HOME/smartSlurm/bin:$PATH  

# Set up a function so that ssbatch is called when running sbatch
sbatch() { $HOME/smartSlurm/bin/ssbatch "$@"; }; export -f sbatch                                 

# Create some text files for testing
createBigTextFiles.sh

# Run 3 jobs to get memory and run-time statistics for useSomeMemTimeNoInput.sh
# Notice 1: you don't have to run this section, because I have run it and save the statistics in $HOME/smartSlurm
# Notice 2: Slurm will submit three jobs to short partition, each reserves 2G memory and 2 hour run-time 
export SSBATCH_S=useSomeMemTimeNoInput.sh # This is optional because the software name is the same as the script
for i in {1..3}; do
    sbatch --mem 2G -t 2:0:0 --wrap="useSomeMemTimeNoInput.sh $i"
done

# After the 3 jobs finish, auto adjust memory and run-time so that 90% jobs can finish successfully
# Notice: Slurm will submit five jobs to short partition, and reserved 19M memory and 7 minute run-time 
sbatch --mem 2G -t 2:0:0 --mem 2G --wrap="useSomeMemTimeNoInput.sh bigText1.txt"

# Run 5 jobs to get memory and run-time statistics for useSomeMemTimeAccordingInputSize.sh
# Notice 1: you don't have to run this section, because I have run it and save the statistics in $HOME/smartSlurm
# Notice 2: Slurm will submit five jobs to short partition, each reserves 2G memory and 2 hour run-time 
export SSBATCH_S=useSomeMemTimeAccordingInputSize.sh # This is optional because the software name is the same as the script
for i in {1..5}; do
    export SSBATCH_I=bigText$i.txt # This is to tell ssbatch the input file to calculate input file size
    sbatch -t 2:0:0 --mem 2G --wrap="useSomeMemTimeAccordingInputSize.sh bigText$i.txt"
done

# After the 5 jobs finish, auto adjust memory and run-time according input file size
# Notice: Slurm will submit one job to short partition, and reserves 21M memory and 13 minute run-time 
export SSBATCH_S=useSomeMemTimeAccordingInputSize.sh # This is optional because the software name is the same as the script
export SSBATCH_I=bigText1.txt,bigText2.txt # This is to tell ssbatch the input file to calculate input file size 
sbatch -t 2:0:0 --mem 2G --wrap="useSomeMemTimeAccordingInputSize.sh bigText1.txt bigText$2.txt"

# The second way to tell the input file name and software name: 
sbatch --comment="SSBATCH_S=useSomeMemTimeAccordingInputSize.sh SSBATCH_I=bigText1.txt,bigText2.txt" \
    -t 2:0:0 --mem 2G --wrap="useSomeMemTimeAccordingInputSize.sh bigText1.txt bigText$2.txt"

# The third way to tell the input file name and software name: 
sbatch -t 2:0:0 --mem 2G job.sh

cat job.sh
#!/bin/bash
#SBATCH --commen="SSBATCH_S=useSomeMemTimeAccordingInputSize.sh SSBATCH_I=bigText1.txt,bigText2.txt"
useSomeMemTimeAccordingInputSize.sh bigText1.txt bigText$2.txt

# After you finish using ssbatch, run this command to disable it:    
unset sbatch
unset SSBATCH_S
unset SSBATCH_I
```

## How does it works

1) Auto adjust memory and run-time according to statistics from earlier jobs
2) 
~/smartSlurm/jobRecord.txt contains job memory and run-time records. There are three important columns: 
   
   1 2rd colume is input size,
   
   2 7th column is final memory usage
   
   3 8th column is final time usage
   
The data from the three columns are plotted and static  
__________________________________________________________________________________________________________________   
1jobID,2inputSize,3mem,4time,5mem,6time,7mem,8time,9status,10useID,11path,12software,13reference,14output,15script,16error,17cpu,18node,19date,20command
46531,1465,4G,2:0:0,4G,0-2:0:0,3.52,1,COMPLETED,ld32,,useSomeMemTimeAccordingInputSize.sh,none,slurm-%j.out slurm-YRTrRAYA.sh slurm-%j.err,1,compute-a-16-21,slurm-46531.err,Tue Dec 6 15:29:20 EST 2022,"ssbatch -p short -t 2:0:0 --mem=4G --wrap useSomeMemTimeAccordingInputSize.sh bigText1.txt run"

46535,2930,4G,2:0:0,4G,0-2:0:0,6.38,2,COMPLETED,ld32,,useSomeMemTimeAccordingInputSize.sh,none,slurm-%j.out slurm-oT42tyEE.sh slurm-%j.err,1,compute-a-16-21,slurm-46535.err,Tue Dec 6 15:30:46 EST 2022,"ssbatch -p short -t 2:0:0 --mem=4G --wrap useSomeMemTimeAccordingInputSize.sh bigText2.txt run"

46534,4395,4G,2:0:0,4G,0-2:0:0,9.24,4,COMPLETED,ld32,,useSomeMemTimeAccordingInputSize.sh,none,slurm-%j.out slurm-TQyBOQ5f.sh slurm-%j.err,1,compute-a-16-21,slurm-46534.err,Tue Dec 6 15:32:40 EST 2022,"ssbatch -p short -t 2:0:0 --mem=4G --wrap useSomeMemTimeAccordingInputSize.sh bigText3.txt run"


3) Auto choose partition according to run-time request

~/smartSlurm/config/partitions.txt contains partion time limit and bash function adjustPartition to adjust partion for sbatch jobs: 

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

4) Auto re-run failed OOM (out of memory) and OOT (out of run-time) jobs


6) Get good emails: by default Slurm emails only have a subject. ssbatch attaches the content of the sbatch script, the output and error log to email





# sbatchAndTop
## How to use sbatchAndTop
```
cd ~    
git clone git@github.com:ld32/smartSlurm.git  
export PATH=$HOME/smartSlurm/bin:$PATH    
sbatchAndTop <sbatch option1> <sbatch option 2> <sbatch option 3> <...> 

# Such as:    
sbatchAndTop -p short -c 1 -t 2:0:0 --mem 2G --wrap "my_application para1 para2" 
# Here -p short is optional, because ssbatch chooses partition according to run time.  

# or:     
sbatchAndTop job.sh 

## sbatchAndTop features:
1) Submit slurm job using ssbatch (scrool up to see ssbatch features) and run scontrol top on the job
