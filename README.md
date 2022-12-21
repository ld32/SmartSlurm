=============================
# SmartSlurm
- [Smart sbatch](#smart-sbatch)
    - [Smart sbatch features](#smart-sbatch-features)
    - [Quick start](#ssbatch-quick-start)
    - [How to use ssbatch](#how-to-use-ssbatch)
    - [How does smart sbatch work]($how-does-smart-sbatch-work)

- [Run bash script as smart pipeline using smart sbatch](#Run-bash-script-as-smart-pipeline-using-smart-sbatch)
    - [Smart pipeline features](#smart-pipeline-features)
    - [smart pipline Quick start](#quick-start)
    - [How to use smart pipeline](#how-to-use-smart-pipeline)
    - [How does smart pipeline work]($how-does-smart-pipeline-work)

- [sbatchAndTop](#sbatchAndTop)


# Smart sbatch
ssbath was originally designed to run https://github.com/ENCODE-DCC/atac-seq-pipeline, so that users don't have to modify the original workflow and sbatch can automatially modify the partitions according user's local cluster partition settings. The script was later improved to have more features.

![](https://github.com/ld32/smartSlurm/blob/main/stats/useSomeMemTimeAccordingInputSize.sh.none.time.png)

## ssbatch features:
1) Auto adjust memory and run-time according to statistics from earlier jobs
2) Auto choose partition according to run-time request
3) Auto re-run failed OOM (out of memory) and OOT (out of run-time) jobs
4) Get good emails: by default Slurm emails only have a subject. ssbatch attaches the content of the sbatch script, the output and error log to email
## ssbatch Quick Start
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

# run sbatch as usual 
# Notice: Slurm will submit a jobs to short partition, and reserved 21M memory and 7 minute run-time 
sbatch --mem 2G -t 2:0:0 --mem 2G --wrap="useSomeMemTimeNoInput.sh 1"

# After you finish using ssbatch, run this command to disable it:    
unset sbatch
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

# After the 3 jobs finish, when submitting more jobs, ssbatch auto adjusts memory and run-time so that 90% jobs can finish successfully
# Notice: Slurm will submit five jobs to short partition, and reserved 19M memory and 7 minute run-time 
sbatch --mem 2G -t 2:0:0 --mem 2G --wrap="useSomeMemTimeNoInput.sh 1"

# Run 5 jobs to get memory and run-time statistics for useSomeMemTimeAccordingInputSize.sh
# Notice 1: you don't have to run this section, because I have run it and save the statistics in $HOME/smartSlurm
# Notice 2: Slurm will submit five jobs to short partition, each reserves 2G memory and 2 hour run-time 
export SSBATCH_S=useSomeMemTimeAccordingInputSize.sh # This is optional because the software name is the same as the script
for i in {1..5}; do
    export SSBATCH_I=bigText$i.txt # This is to tell ssbatch the input file to calculate input file size
    sbatch -t 2:0:0 --mem 2G --wrap="useSomeMemTimeAccordingInputSize.sh bigText$i.txt"
done

# After the 5 jobs finish, when submitting more jobs, ssbatch auto adjusts memory and run-time according input file size
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

## How does it work

1) Auto adjust memory and run-time according to statistics from earlier jobs

~/smartSlurm/jobRecord.txt contains job memory and run-time records. There are three important columns: 
   
   2rd colume is input size,
   
   7th column is final memory usage
   
   8th column is final time usage
   
The data from the three columns are plotted and statistics  
__________________________________________________________________________________________________________________   
1jobID,2inputSize,3mem,4time,5mem,6time,7mem,8time,9status,10useID,11path,12software,13reference,14output,15script,16error,17cpu,18node,19date,20command
46531,1465,4G,2:0:0,4G,0-2:0:0,3.52,1,COMPLETED,ld32,,useSomeMemTimeAccordingInputSize.sh,none,slurm-%j.out slurm-YRTrRAYA.sh slurm-%j.err,1,compute-a-16-21,slurm-46531.err,Tue Dec 6 15:29:20 EST 2022,"ssbatch -p short -t 2:0:0 --mem=4G --wrap useSomeMemTimeAccordingInputSize.sh bigText1.txt run"

46535,2930,4G,2:0:0,4G,0-2:0:0,6.38,2,COMPLETED,ld32,,useSomeMemTimeAccordingInputSize.sh,none,slurm-%j.out slurm-oT42tyEE.sh slurm-%j.err,1,compute-a-16-21,slurm-46535.err,Tue Dec 6 15:30:46 EST 2022,"ssbatch -p short -t 2:0:0 --mem=4G --wrap useSomeMemTimeAccordingInputSize.sh bigText2.txt run"

46534,4395,4G,2:0:0,4G,0-2:0:0,9.24,4,COMPLETED,ld32,,useSomeMemTimeAccordingInputSize.sh,none,slurm-%j.out slurm-TQyBOQ5f.sh slurm-%j.err,1,compute-a-16-21,slurm-46534.err,Tue Dec 6 15:32:40 EST 2022,"ssbatch -p short -t 2:0:0 --mem=4G --wrap useSomeMemTimeAccordingInputSize.sh bigText3.txt run"

\#Here is the input size vs memory plot for useSomeMemTimeAccordingInputSize.sh: 

![](https://github.com/ld32/smartSlurm/blob/main/stats/useSomeMemTimeAccordingInputSize.sh.none.mem.png)

\#Here is the input size vs run-time plot for useSomeMemTimeAccordingInputSize.sh: 

![](https://github.com/ld32/smartSlurm/blob/main/stats/useSomeMemTimeAccordingInputSize.sh.none.time.png)

\#Here is the run-time vs memory plot for useSomeMemTimeNoInput.sh: 

![](https://github.com/ld32/smartSlurm/blob/main/stats/useSomeMemTimeNoInput.sh.none.stat.noInput.png)

2) Auto choose partition according to run-time request

~/smartSlurm/config/partitions.txt contains partion time limit and bash function adjustPartition to adjust partion for sbatch jobs: 

\# Genernal partions, ordered by maximum allowed run-time in hours 

partition1Name=short   
partition1TimeLimit=12  # run-time > 0 hours and <= 12 hours 
partition2Name=medium  
partition2TimeLimit=120 # run-time > 12 hours and <= 5 days
partition3Name=long        
partition3TimeLimit=720 # run-time > 5 days and <= 30 days
partition4Name=      
partition4TimeLimit=  
partition5Name=     
partition5TimeLimit=    

\# Special pertitions with special restrictions

partition6Name=priority    # only allow two job running at the same time        
partition6TimeLimit=720 # run-time <= 30 days
partition7Name=highmem     # run-time <= 30 days, special partision     
partition7TimeLimit=720 # run-time <= 30 days 
partition8Name=interactive      
partition8TimeLimit=12  # run-time <= 12 hours
partition9Name=mpi      
partition9TimeLimit=720 # run-time <= 30 days 
partition10Name=        
partition10TimeLimit=       

\#function 

adjustPartition() {         
    ... # please open the file to see the content         
} ; export -f adjustPartition 

3) Auto re-run failed OOM (out of memory) and OOT (out of run-time) jobs
    
    At end of the job, ~/smartSlurm/bin/cleanUp.sh checkes memory and time usage, save the data in to log ~/smartSlurm/myJobRecord.txt. If job fails, ssbatch re-submit with double memory or double time, clear up the statistic fomular, so that later jobs will re-caculate statistics, 

6) Get good emails: by default Slurm emails only have a subject. ssbatch attaches the content of the sbatch script, the output and error log to email

    ~/smartSlurm/bin/cleanUp.sh also sends a email to user. The email contains the content of the Slurm script, the sbatch command used, and also the content of the output and error log files.


# Smart pipeline
Smart pipeline was originally designed to run bash script as pipelie on Slurm cluster. We added dynamic memory/run-time feature to it and now call it Smart pipine. The runAsPipeline script converts an input bash script to a pipeline that easily submits jobs to the Slurm scheduler for you.

## smart pipeline features:
1) Submit each step as a cluster job using ssbatch, which auto adjusts memory and run-time according to statistics from earlier jobs, and re-run OOM/OOT jobs with doubled memory/run-time
2) Automatically arrange dependencies among jobs
3) Email notifications are sent when each job fails or succeeds
4) If a job fails, all its downstream jobs automatically are killed
5) When re-running the pipeline on the same data folder, if there are any unfinished jobs, the user is asked to kill them or not
6) When re-running the pipeline on the same data folder, the user is asked to confirm to re-run or not if a job or a step was done successfully earlier
7) For re-run, if the script is not changed, runAsPipeline does not re-process the bash script and directly use old one
8) If user has more than one Slurm account, adding -A or —account= to command line to let all jobs to use that Slurm account
9) When adding new input data and re-run the workflow, affected successfully finished jobs will be auto re-run.Run bash script as smart slurm pipeline

## Smart Pipeline Quick Start
``` bash
# Download
cd $HOME 
git clone git://github.com/ld32/smartSlurm.git  

# Setup path
export PATH=$HOME/smartSlurm/bin:$PATH  

# Take a look at the exmplar bash script
cat ~/smartSlurm/bin/bashScriptV1.sh

# Below is the content of bashScriptV1.sh
1 #!/bin/sh
2 
3 for i in {1..1}; do
4    input=bigText$i.txt
5    output=1234.$i.txt
6    useSomeMemTimeAccordingInputSize.sh $input; grep 1234 $input > $output
7
8    output=5678.$i.txt
9    useSomeMemTimeAccordingInputSize.sh $input; grep 5678 $input > $output
10 done
11
12 input=bigText1.txt
13 output=all.txt
14 useSomeMemTimeAccordingInputSize.sh $input; cat 1234.*.txt 5678.*.txt > $output

#Notes about bashScriptV1.sh: 
The script first finds 1234 in file bigText1.txt in row 6, then finds 5678 in bigText1.txt in row 9, then merges the results into all.txt in orow 14 
In order tell smart pipeline which step/command we want to submit as Slurm jobs, we add comments above the commands also some helping commands:  

cat ~/smartSlurm/bin/bashScriptV2.sh

# below is the content of bashScriptV1.sh
1 #!/bin/sh
2 
3 outputs=""
4 for i in {1..1}; do
5    input=bigText$i.txt
6    output=1234.$i.txt
7    #@1,0,useSomeMemTimeAccordingInputSize.sh,,input,sbatch -p short -c 1 --mem 2G -t 50:0
8    useSomeMemTimeAccordingInputSize.sh $input; grep 1234 $input > $output
9    outputs=$outputs,$output
10
11    output=5678.$i.txt
12    #@2,0,useSomeMemTimeAccordingInputSize.sh,,input,sbatch -p short -c 1 --mem 2G -t 50:0
13    useSomeMemTimeAccordingInputSize.sh $input; grep 5678 $input > $output
14    outputs=$outputs,$output
15 done
16 
17 input1=bigText1.txt
18 output=all.txt
19 #@3,1.2,useSomeMemTimeAccordingInputSize.sh,,input
20 useSomeMemTimeAccordingInputSize.sh $input; cat 1234.*.txt 5678.*.txt > $output    

# Notice that there are a few things added to the script here:

    Step 1 is denoted by #@1,0,useSomeMemTimeAccordingInputSize.sh,,input,sbatch -p short -c 1 --mem 2G -t 50:0 (line 7 above), which means this is step 1 that depends on no other step, run software useSomeMemTimeAccordingInputSize.sh, does not use any reference files, and file $input is the input file, needs to be copied to the /tmp directory if user want to use /tmp. The sbatch command tells the pipeline runner the sbatch parameters to run this step.

    Step 2 is denoted by #@2,0,useSomeMemTimeAccordingInputSize.sh,,input,sbatch -p short -c 1 --mem 2G -t 50:0 (line 12 above), which means this is step2 that depends on no other step, run software useSomeMemTimeAccordingInputSize.sh, does not use any reference file, and file $input is the input file, needs be copy to /tmp directory if user wants to use /tmp. The sbatch command tells the pipeline runner the sbatch parameters to run this step.  

    Step 3 is denoted by #@3,1.2,useSomeMemTimeAccordingInputSize.sh,,input (line 19), which means that this is step3 that depends on step1 and step2, and the step runs software useSomeMemTimeAccordingInputSize.sh with on reference file, and use $input as input file. Notice, there is no sbatch here,  so the pipeline runner will use default sbatch command from command line (see below).   

Notice the format of step annotation is #@stepID,dependIDs,stepName,reference,input,sbatchOptions. Reference is optional, which allows the pipeline runner to copy data (file or folder) to local /tmp folder on the computing node to speed up the software. Input is optional, which is used to estimate memory/run-time for the job. sbatchOptions is also optional, and when it is missing, the pipeline runner will use the default sbatch command given from command line (see below).

Here are two more examples:

#@4,1.3,map,,in,sbatch -p short -c 1 -t 50:0   Means step4 depends on step1 and step3, this step is named map, there is no reference data to copy, there is input $in and submits this step with sbatch -p short -c 1 -t 50:0

#@3,1.2,align,db1.db2   Means step3 depends on step1 and step2, this step is named align, $db1 and $db2 are reference data to be copied to /tmp , there is no input and submit with the default sbatch command (see below).



# run sbatch as usual 
# Notice: Slurm will submit a jobs to short partition, and reserved 21M memory and 7 minute run-time 
runAsPipeline ~/smartSlurm/bin/bashScriptV3.sh "sbatch -A rccg -p short -t 10:0 -c 1" noTmp run


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

# After the 3 jobs finish, when submitting more jobs, ssbatch auto adjusts memory and run-time so that 90% jobs can finish successfully
# Notice: Slurm will submit five jobs to short partition, and reserved 19M memory and 7 minute run-time 
sbatch --mem 2G -t 2:0:0 --mem 2G --wrap="useSomeMemTimeNoInput.sh 1"

# Run 5 jobs to get memory and run-time statistics for useSomeMemTimeAccordingInputSize.sh
# Notice 1: you don't have to run this section, because I have run it and save the statistics in $HOME/smartSlurm
# Notice 2: Slurm will submit five jobs to short partition, each reserves 2G memory and 2 hour run-time 
export SSBATCH_S=useSomeMemTimeAccordingInputSize.sh # This is optional because the software name is the same as the script
for i in {1..5}; do
    export SSBATCH_I=bigText$i.txt # This is to tell ssbatch the input file to calculate input file size
    sbatch -t 2:0:0 --mem 2G --wrap="useSomeMemTimeAccordingInputSize.sh bigText$i.txt"
done

# After the 5 jobs finish, when submitting more jobs, ssbatch auto adjusts memory and run-time according input file size
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

## How does it work

1) Auto adjust memory and run-time according to statistics from earlier jobs

~/smartSlurm/jobRecord.txt contains job memory and run-time records. There are three important columns: 
   
   2rd colume is input size,
   
   7th column is final memory usage
   
   8th column is final time usage
   
The data from the three columns are plotted and statistics  
__________________________________________________________________________________________________________________   
1jobID,2inputSize,3mem,4time,5mem,6time,7mem,8time,9status,10useID,11path,12software,13reference,14output,15script,16error,17cpu,18node,19date,20command
46531,1465,4G,2:0:0,4G,0-2:0:0,3.52,1,COMPLETED,ld32,,useSomeMemTimeAccordingInputSize.sh,none,slurm-%j.out slurm-YRTrRAYA.sh slurm-%j.err,1,compute-a-16-21,slurm-46531.err,Tue Dec 6 15:29:20 EST 2022,"ssbatch -p short -t 2:0:0 --mem=4G --wrap useSomeMemTimeAccordingInputSize.sh bigText1.txt run"

46535,2930,4G,2:0:0,4G,0-2:0:0,6.38,2,COMPLETED,ld32,,useSomeMemTimeAccordingInputSize.sh,none,slurm-%j.out slurm-oT42tyEE.sh slurm-%j.err,1,compute-a-16-21,slurm-46535.err,Tue Dec 6 15:30:46 EST 2022,"ssbatch -p short -t 2:0:0 --mem=4G --wrap useSomeMemTimeAccordingInputSize.sh bigText2.txt run"

46534,4395,4G,2:0:0,4G,0-2:0:0,9.24,4,COMPLETED,ld32,,useSomeMemTimeAccordingInputSize.sh,none,slurm-%j.out slurm-TQyBOQ5f.sh slurm-%j.err,1,compute-a-16-21,slurm-46534.err,Tue Dec 6 15:32:40 EST 2022,"ssbatch -p short -t 2:0:0 --mem=4G --wrap useSomeMemTimeAccordingInputSize.sh bigText3.txt run"

\#Here is the input size vs memory plot for useSomeMemTimeAccordingInputSize.sh: 

![](https://github.com/ld32/smartSlurm/blob/main/stats/useSomeMemTimeAccordingInputSize.sh.none.mem.png)

\#Here is the input size vs run-time plot for useSomeMemTimeAccordingInputSize.sh: 

![](https://github.com/ld32/smartSlurm/blob/main/stats/useSomeMemTimeAccordingInputSize.sh.none.time.png)

\#Here is the run-time vs memory plot for useSomeMemTimeNoInput.sh: 

![](https://github.com/ld32/smartSlurm/blob/main/stats/useSomeMemTimeNoInput.sh.none.stat.noInput.png)

2) Auto choose partition according to run-time request

~/smartSlurm/config/partitions.txt contains partion time limit and bash function adjustPartition to adjust partion for sbatch jobs: 

\# Genernal partions, ordered by maximum allowed run-time in hours 

partition1Name=short   
partition1TimeLimit=12  # run-time > 0 hours and <= 12 hours 
partition2Name=medium  
partition2TimeLimit=120 # run-time > 12 hours and <= 5 days
partition3Name=long        
partition3TimeLimit=720 # run-time > 5 days and <= 30 days
partition4Name=      
partition4TimeLimit=  
partition5Name=     
partition5TimeLimit=    

\# Special pertitions with special restrictions

partition6Name=priority    # only allow two job running at the same time        
partition6TimeLimit=720 # run-time <= 30 days
partition7Name=highmem     # run-time <= 30 days, special partision     
partition7TimeLimit=720 # run-time <= 30 days 
partition8Name=interactive      
partition8TimeLimit=12  # run-time <= 12 hours
partition9Name=mpi      
partition9TimeLimit=720 # run-time <= 30 days 
partition10Name=        
partition10TimeLimit=       

\#function 

adjustPartition() {         
    ... # please open the file to see the content         
} ; export -f adjustPartition 

3) Auto re-run failed OOM (out of memory) and OOT (out of run-time) jobs
    
    At end of the job, ~/smartSlurm/bin/cleanUp.sh checkes memory and time usage, save the data in to log ~/smartSlurm/myJobRecord.txt. If job fails, ssbatch re-submit with double memory or double time, clear up the statistic fomular, so that later jobs will re-caculate statistics, 

6) Get good emails: by default Slurm emails only have a subject. ssbatch attaches the content of the sbatch script, the output and error log to email

    ~/smartSlurm/bin/cleanUp.sh also sends a email to user. The email contains the content of the Slurm script, the sbatch command used, and also the content of the output and error log files.



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
