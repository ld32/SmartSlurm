=============================
# SmartSlurm

- [Smart sbatch](#smart-sbatch)
    - [ssbatch features](#ssbatch-features)
    - [How to use ssbatch](#how-to-use-ssbatch)
    - [How does ssbatch work](#how-does-ssbatch-work) 
      
- [Use ssbatch in Snakemake pipeline](#Use-ssbatch-in-Snakemake-pipeline)

- [Use ssbatch in Cromwell pipeline](#Use-ssbatch-in-Cromwell-pipeline])

- [Use ssbatch in Nextflow pipeline](#Use-ssbatch-in-Nextflow-pipeline)

- [Run bash script as smart pipeline using smart sbatch](#Run-bash-script-as-smart-pipeline-using-smart-sbatch)
    - [Smart pipeline features](#smart-pipeline-features)
    - [How to use smart pipeline](#how-to-use-smart-pipeline)
    - [How does smart pipeline work](#how-does-smart-pipeline-work)

- [sbatchAndTop](#sbatchAndTop)


# Smart sbatch
[Back to top](#SmartSlurm)

ssbath was originally designed to run https://github.com/ENCODE-DCC/atac-seq-pipeline, so that users don't have to modify the original workflow and ssbatch can automatially modify the partitions according user's local cluster partition settings. The script was later improved to have more features.

![](https://github.com/ld32/smartSlurm/blob/main/stats/back/useMemTimeWithInput.none.time.png)

##As the figure shown above, the memory usage is roughly co-related to the input size. We can use input size to allocate memory when submit new jobs.

![](https://github.com/ld32/smartSlurm/blob/main/stats/back/barchartMem.png)

##As the figure shown above, Smart Slurm runs the first 5 jobs, it use default memory, then based on the first five jobs, it estimates memory for future jobs. The wasted memory is dramatially decreased for the future jobs.

![](https://github.com/ld32/smartSlurm/blob/main/stats/back/barchartTime.png)

##As the figure shown above, Smart Slurm run the first 5 jobs, it use default time, then based on the first five jobs, it estimates time future jobs. The wasted time is dramatially decreased for the future jobs.

## ssbatch features:
[Back to top](#SmartSlurm)

1) Auto adjust memory and run-time according to statistics from earlier jobs
2) Auto choose partition according to run-time request
3) Auto re-run failed OOM (out of memory) and OOT (out of run-time) jobs
4) (Optional) Generate checkpoint before job runs out of time or memory, and use the checkpoint to re-run jobs.
5) Get good emails: by default Slurm emails only have a subject. ssbatch attaches the content of the sbatch script, the output and error log to email

## How to use ssbatch
[Back to top](#SmartSlurm)

``` bash
# Download 
git clone https://github.com/ld32/smartSlurm.git  

# Setup path
export PATH=$PWD/smartSlurm/bin:$PATH  

# Create some text files for testing
createBigTextFiles.sh

# Run 3 jobs to get memory and run-time statistics for useMemTimeNoInput
for i in {1..3}; do
    ssbatch --mem 2G -t 2:0:0 -S useMemTimeNoInput --wrap="useMemTimeNoInput.sh $i"
done

# After the 3 jobs finish, when submitting more jobs, ssbatch auto adjusts memory 
# and run-time so that 90% jobs can finish successfully
# Notice: this command submits this job to short partition, and reserves 19M memory and 7 minute run-time 
ssbatch --mem 2G -t 2:0:0 --mem 2G -S useMemTimeNoInput --wrap="useMemTimeNoInput.sh 1"

# Run 3 jobs to get memory and run-time statistics for useMemTimeWithInput
for i in {1..3}; do
    ssbatch -t 2:0:0 --mem 2G -S useMemTimeWithInput -I bigText$i.txt \
        --wrap="useMemTimeWithInput.sh bigText$i.txt"
done

# After the 5 jobs finish, when submitting more jobs, ssbatch auto adjusts memory and run-time according input file size
# Notice: this command submits the job to short partition, and reserves 21M memory and 13 minute run-time 
ssbatch -t 2:0:0 --mem 2G -S useMemTimeWithInput \
    -I "bigText1.txt,bigText2.txt" --wrap="useMemTimeWithInput.sh bigText1.txt bigText2.txt"

# The second way to tell the input file name: 
sbatch -t 2:0:0 --mem 2G job.sh

cat job.sh
#!/bin/bash
#SBATCH --commen="S=useMemTimeWithInput I=bigText1.txt,bigText2.txt"
useMemTimeWithInput.sh bigText1.txt bigText$2.txt

```

## How does ssbatch work    
[Back to top](#SmartSlurm)

1) Auto adjust memory and run-time according to statistics from earlier jobs

$smartSlurmJobRecordDir/jobRecord.txt contains job memory and run-time records. There are three important columns: 
   
   1st colume is job ID
   
   2rd colume is input size
   
   7th column is actual memory usage
   
   8th column is actual time usage
   
The data from the three columns are plotted and statistics  
__________________________________________________________________________________________________________________   
1jobID,2inputSize,3mem,4time,5mem,6time,7mem,8time,9status,10useID,11path,12software,13reference

46531,1465,4G,2:0:0,4G,0-2:0:0,3.52,1,COMPLETED,ld32,,useMemTimeWithInput,none

46535,2930,4G,2:0:0,4G,0-2:0:0,6.38,2,COMPLETED,ld32,,useMemTimeWithInput,none

46534,4395,4G,2:0:0,4G,0-2:0:0,9.24,4,COMPLETED,ld32,,useMemTimeWithInput,none

\#Here is the input size vs memory plot for useMemTimeWithInput: 

![](https://github.com/ld32/smartSlurm/blob/main/stats/back/useMemTimeWithInput.none.mem.png)

\#Here is the input size vs run-time plot for useMemTimeWithInput: 

![](https://github.com/ld32/smartSlurm/blob/main/stats/back/useMemTimeWithInput.none.time.png)

\#Here is the run-time vs memory plot for useMemTimeNoInput: 

![](https://github.com/ld32/smartSlurm/blob/main/stats/back/useMemTimeNoInput.none.stat.noInput.png)

2) Auto choose partition according to run-time request

smartSlrm/config/config.txt contains partion time limit and bash function adjustPartition to adjust partion for sbatch jobs: 

\# Genernal partions, ordered by maximum allowed run-time in hours 

partition1Name=short; partition1TimeLimit=12  # run-time > 0 hours and <= 12 hours 

partition2Name=medium; partition2TimeLimit=120 # run-time > 12 hours and <= 5 days

partition3Name=long; partition3TimeLimit=720 # run-time > 5 days and <= 30 days

...        

\#function 

adjustPartition() {         
    ... # please open the file to see the content         
} ; export -f adjustPartition 

3) Auto re-run failed OOM (out of memory) and OOT (out of run-time) jobs
    
    At end of the job, $smartSlurmJobRecordDir/bin/cleanUp.sh checkes memory and time usage, save the data in to log $smartSlurmJobRecordDir/myJobRecord.txt. If job fails, ssbatch re-submit with double memory or double time, clear up the statistic fomular, so that later jobs will re-caculate statistics, 

4) Checkpoint
    
    If checkpoint is enabled, before the job run out of memory or time, ssbatch generate checkpoint and resubmit the job.

5) Get good emails: by default Slurm emails only have a subject. ssbatch attaches the content of the sbatch script, the output and error log to email

    $smartSlurmJobRecordDir/bin/cleanUp.sh also sends a email to user. The email contains the content of the Slurm script, the sbatch command used, and also the content of the output and error log files.


# Use ssbatch in Snakemake pipeline
[Back to top](#SmartSlurm)

``` bash
# Download smartSlurm if it is not done yet 
git clone https://github.com/ld32/smartSlurm.git  

# Download snakemake tutorial (from: https://snakemake.readthedocs.io/en/v3.11.0/tutorial/setup.html)
wget https://bitbucket.org/snakemake/snakemake-tutorial/get/v3.9.0-1.tar.bz2
tar -xf v3.9.0-1.tar.bz2 --strip 1
cp $PWD/smartSlurm/bin/Snakefile .
cp $PWD/smartSlurm/bin/config.yaml .

# Create Snakemake conda env (from: https://snakemake.readthedocs.io/en/v3.11.0/tutorial/setup.html)
module load miniconda3
mamba env create --name snakemakeEnv --file $PWD/smartSlurm/config/snakemakeEnv.yaml

# Review Snakefile, activate the snakemake env and run test
module load miniconda3
source activate snakemakeEnv
export PATH=$PWD/smartSlurm/bin:$PATH
cat Snakefile
snakemake -p -j 999 --latency-wait=80 --cluster "ssbatch -t 100 --mem 1G -p short"

If you have multiple Slurm account:
snakemake -p -j 999 --latency-wait=80 --cluster "ssbatch -A mySlurmAccount -t 100 --mem 1G"

```

# Use ssbatch in Cromwell pipeline
[Back to top](#SmartSlurm)

``` bash
# Download 
git clone https://github.com/ld32/smartSlurm.git  

# Setup path
export PATH=$PWD/smartSlurm/bin:$PATH  

# Create some text files for testing
createBigTextFiles.sh

# Use ssbatch to replace regular sbatch, set up a fuction. 
# This way, whenever you run sbatch, ssbatch is called. 
sbatch() { $HOME/smartSlurm/bin/ssbatch "$@"; }; export -f sbatch                                 

# Run 3 jobs to get memory and run-time statistics for useMemTimeNoInput
for i in {1..3}; do
    sbatch --mem 2G -t 2:0:0 --commen="S=useMemTimeNoInput" \
        --wrap="useMemTimeNoInput.sh $i"
done

# After the 3 jobs finish, when submitting more jobs, ssbatch auto adjusts memory 
# and run-time so that 90% jobs can finish successfully
# Notice: this command submits this job to short partition, and reserves 19M memory and 7 minute run-time 
sbatch --mem 2G -t 2:0:0 --mem 2G --commen="S=useMemTimeNoInput" --wrap="useMemTimeNoInput.sh 1"

# Run 3 jobs to get memory and run-time statistics for useMemTimeWithInput
for i in {1..3}; do
    sbatch -t 2:0:0 --mem 2G --commen="S=useMemTimeWithInput I=bigText$i.txt" \
        --wrap="useMemTimeWithInput.sh bigText$i.txt"
done

# After the 5 jobs finish, when submitting more jobs, ssbatch auto adjusts memory and run-time according input file size
# Notice: this command submits the job to short partition, and reserves 21M memory and 13 minute run-time 
sbatch -t 2:0:0 --mem 2G --commen="S=useMemTimeWithInput \
    I=bigText1.txt,bigText2.txt" --wrap="useMemTimeWithInput.sh bigText1.txt bigText2.txt"

# The second way to tell the input file name: 
sbatch -t 2:0:0 --mem 2G job.sh

cat job.sh
#!/bin/bash
#SBATCH --commen="S=useMemTimeWithInput I=bigText1.txt,bigText2.txt"
useMemTimeWithInput.sh bigText1.txt bigText$2.txt

# After you finish using ssbatch, run these command to disable ssbatch:    
unset sbatch

```

# Use ssbatch in Nextflow pipeline
[Back to top](#SmartSlurm)

``` bash
# Download smartSlurm if it is not done yet 
git clone https://github.com/ld32/smartSlurm.git  

# Create Nextflow conda env
module load miniconda3
mamba create -n  nextflowEnv -c bioconda -y nextflow

# Review nextflow file, activate the nextflow env and run test
module load miniconda3
export PATH=$PWD/smartSlurm/bin:$PATH  
source activate nextflowEnv
cp $PWD/smartSlurm/bin/nextflow.nf .
cp $PWD/smartSlurm/config/nextflow.config .

# If you have multiple Slurm account, modify the config file:
nano nextflow.config

#change:
//process.clusterOptions = '--account=mySlurmAcc'
to 
process.clusterOptions = '--account=mySlurmAcc'

# save the file 

# Ready to run:
nextflow run nextflow.nf -profile slurm

```
# Run bash script as smart pipeline using smart sbatch
[Back to top](#SmartSlurm)

Smart pipeline was originally designed to run bash script as pipelie on Slurm cluster. We added dynamic memory/run-time feature to it and now call it Smart pipeline. The runAsPipeline script converts an input bash script to a pipeline that easily submits jobs to the Slurm scheduler for you.

\#Here is the memeory usage by the optimzed workflow: The orignal pipeline has 11 steps. Most of the steps only need less than 10G memory to run. But one of the step needs 140G. Because the original pipeline is submitted as a single huge job, 140G is reserved for all the steps. (Each compute node in the cluster has 256 GB RAM.) By submitting each step as a separate job, most steps only need to reserve 10G, which decreases the memory usage dramatically. (The pink part of the graph below shows these savings.) Another optimization is to dynamically allocate memory based on the reference genome size and input sequencing data size. (This in shown in the yellow part of the graph.)
Because of the decreased resource demand, the jobs can start earlier, and in turn increase the overall throughput.

![](https://github.com/ld32/smartSlurm/blob/main/stats/back/barchartMemSaved.png)

## smart pipeline features:
[Back to top](#SmartSlurm)

1) Submit each step as a cluster job using ssbatch, which auto adjusts memory and run-time according to statistics from earlier jobs, and re-run OOM/OOT jobs with doubled memory/run-time
2) Automatically arrange dependencies among jobs
3) Email notifications are sent when each job fails or succeeds
4) If a job fails, all its downstream jobs automatically are killed
5) When re-running the pipeline on the same data folder, if there are any unfinished jobs, the user is asked to kill them or not
6) When re-running the pipeline on the same data folder, the user is asked to confirm to re-run or not if a job or a step was done successfully earlier
7) For re-run, if the script is not changed, runAsPipeline does not re-process the bash script and directly use old one
8) If user has more than one Slurm account, adding -A or —account= to command line to let all jobs to use that Slurm account
9) When adding new input data and re-run the workflow, affected successfully finished jobs will be auto re-run.Run bash script as smart slurm pipeline

## How to use smart pipeline
[Back to top](#SmartSlurm)

``` bash
# Download
git clone https://github.com/ld32/smartSlurm.git  

# Setup path
export PATH=$PWD/smartSlurm/bin:$PATH  

# Take a look at a regular exmaple bash script
cat $PWD/smartSlurm/scripts/bashScriptV1.sh

# Below is the content of a regular bashScriptV1.sh 
1 #!/bin/sh
2
3 for i in {1..1}; do
4    input=bigText$i.txt
5    output=1234.$i.txt
6    useMemTimeWithInput.sh $input; grep 1234 $input > $output
7
8    output=5678.$i.txt
9    useMemTimeWithInput.sh $input; grep 5678 $input > $output
10 done
11
12 input=bigText1.txt
13 output=all.txt
14 useMemTimeWithInput.sh $input; cat 1234.*.txt 5678.*.txt > $output

# Notes about bashScriptV1.sh: 
The script first finds 1234 in file bigText1.txt in row 6, then finds 
5678 in bigText1.txt in row 9, then merges the results into all.txt in orow 14 

# In order tell smart pipeline which step/command we want to submit as Slurm jobs, 
# we add comments above the commands also some helping commands:  
cat $smartSlurmJobRecordDir/scripts/bashScriptV2.sh

# below is the content of bashScriptV2.sh
1 #!/bin/sh
2 
3 outputs=""
4
5 for i in {1..1}; do
6    input=bigText$i.txt
7    output=1234.$i.txt
8    #@1,0,useMemTimeWithInput,,input,sbatch -p short -c 1 --mem 2G -t 2:0:0
9    useMemTimeWithInput.sh $input; grep 1234 $input > $output
10   outputs=$outputs,$output
11
12   output=5678.$i.txt
13   #@2,0,useMemTimeWithInput,,input,sbatch -p short -c 1 --mem 2G -t 2:0:0
14   useMemTimeWithInput.sh $input; grep 5678 $input > $output
15   outputs=$outputs,$output
16 done
17 
18 input=bigText1.txt
19 output=all.txt
20 #@3,1.2,useMemTimeWithInput,,input
21 useMemTimeWithInput.sh $input; cat 1234.*.txt 5678.*.txt > $output    
```

## Notice that there are a few things added to the script here:
[Back to top](#SmartSlurm)

Before the for loop start, there is #loopStart:i, which means all the steps inside the loop use $i as part of unique job identifier.

Step 1 is denoted by #@1,0,useMemTimeWithInput,,input,sbatch -p short -c 1 --mem 2G -t 2:0:0 (line 7 above), which means this is step 1 that depends on no other step, run software useMemTimeWithInput, use the value of $i as unique job identifier for this this step, does not use any reference files, and file $input is the input file, needs to be copied to the /tmp directory if user want to use /tmp. The sbatch command tells the pipeline runner the sbatch parameters to run this step.

Step 2 is denoted by #@2,0,useMemTimeWithInput,,input,sbatch -p short -c 1 --mem 2G -t 2:0:0 (line 12 above), which means this is step2 that depends on no other step, run software useMemTimeWithInput, use the value of $i as unique job identifier for this step, does not use any reference file, and file $input is the input file, needs be copy to /tmp directory if user wants to use /tmp. The sbatch command tells the pipeline runner the sbatch parameters to run this step.  

Step 3 is denoted by #@3,1.2,useMemTimeWithInput,,input (line 19), which means that this is step3 that depends on step1 and step2, and the step runs software useMemTimeWithInput with no reference file, does not need unique identifier because there is only one job in the step, and use $input as input file. Notice, there is no sbatch here,  so the pipeline runner will use default sbatch command from command line (see below).   

Notice the format of step annotation is #@stepID,dependIDs,sofwareName,reference,input,sbatchOptions. Reference is optional, which allows the pipeline runner to copy data (file or folder) to local /tmp folder on the computing node to speed up the software. Input is optional, which is used to estimate memory/run-time for the job. sbatchOptions is also optional, and when it is missing, the pipeline runner will use the default sbatch command given from command line (see below).

Here are two more examples:

#@4,1.3,map,,in,sbatch -p short -c 1 -t 2:0:0  #Means step4 depends on step1 and step3, this step run software 'map', there is no reference data to copy, there is input $in and submits this step with sbatch -p short -c 1 -t 2:0:0

#@3,1.2,align,db1.db2   # Means step3 depends on step1 and step2, this step run software 'align', $db1 and $db2 are reference data to be copied to /tmp , there is no input and submit with the default sbatch command (see below).

# Test run the modified bash script as a pipeline
[Back to top](#SmartSlurm)

runAsPipeline bashScriptV2.sh "sbatch -p short -t 10:0 -c 1" useTmp

This command will generate new bash script of the form slurmPipeLine.checksum.sh in log folder. The checksum portion of the filename will have a MD5 hash that represents the file contents. We include the checksum in the filename to detect when script contents have been updated. If it is not changed, we don not re-create the pipeline script.

This runAsPipeline command will run a test of the script, meaning does not really submit jobs. It will only show a fake job ids like 1234 for each step. If you were to append run at the end of the command, the pipeline would actually be submitted to the Slurm scheduler.

Ideally, with useTmp, the software should run faster using local /tmp disk space for database/reference than the network storage. For this small query, the difference is small, or even slower if you use local /tmp. If you don't need /tmp, you can use noTmp.

With useTmp, the pipeline runner copy related data to /tmp and all file paths will be automatically updated to reflect a file's location in /tmp when using the useTmp option. 
Sample output from the test run

Note that only step 2 used -t 2:0:0, and all other steps used the default -t 10:0. The default walltime limit was set in the runAsPipeline command, and the walltime parameter for step 2 was set in the bash_script_v2.sh script.
runAsPipeline bashScriptV2.sh "sbatch -p short -t 10:0 -c 1" useTmp


# here is the outputs:

Wed Dec 21 15:50:43 EST 2022
Running: /home/ld32/smartSlurm/bin/runAsPipeline /home/ld32/smartSlurm/scripts/bashScriptV2.sh sbatch -p short -t 10:0 -c 1 noTmp

Currently Loaded Modules:
  1) gcc/6.2.0


converting /home/ld32/smartSlurm/scripts/bashScriptV2.sh to $smartSlurmLogDir/slurmPipeLine.6f93dc8953b9c1d1f96b4fabd657446a.sh

find loop start: for i in {1..1}; do

find job marker:
#@1,0,useMemTimeWithInput,,input,sbatch -p short -c 1 --mem 2G -t 2:0:0
sbatch options: sbatch -p short -c 1 --mem 2G -t 2:0:0

find job:
useMemTimeWithInput.sh $input; grep 1234 $input > $output

find job marker:
#@2,0,useMemTimeWithInput,,input,sbatch -p short -c 1 --mem 2G -t 2:0:0
sbatch options: sbatch -p short -c 1 --mem 2G -t 2:0:0

find job:
useMemTimeWithInput.sh $input; grep 5678 $input > $output
find loop end: done

find job marker:
#@3,1.2,useMemTimeWithInput,,input

find job:
useMemTimeWithInput.sh $input; cat 1234.*.txt 5678.*.txt > $output
smartSlurmLog/slurmPipeLine.6f93dc8953b9c1d1f96b4fabd657446a.sh /home/ld32/smartSlurm/scripts/bashScriptV2.sh is ready to run. Starting to run ...
Running $smartSlurmLogDir/slurmPipeLine.6f93dc8953b9c1d1f96b4fabd657446a.sh /home/ld32/smartSlurm/scripts/bashScriptV2.sh

Currently Loaded Modules:
  1) gcc/6.2.0

---------------------------------------------------------

step: 1, depends on: 0, job name: useMemTimeWithInput, flag: useMemTimeWithInput.1
Running:
ssbatch -L /home/ld32/smartSlurm -S useMemTimeWithInput -R none -F 1.0.useMemTimeWithInput.1 -I ,bigText1.txt -D null -p short -c 1 --me
m 2G -t 2:0:0 --wrap "set -e; useMemTimeWithInput.sh bigText1.txt; grep 1234 bigText1.txt > 1234.1.txt"

Parsing result from sbatch commandline:
sbatch options: partition: short time: 2:0:0 mem: 2G mem-per-cpu: task: core: 1 node: out: err: dep:
wrapCMD: set -e; useMemTimeWithInput.sh bigText1.txt; grep 1234 bigText1.txt > 1234.1.txt
additional sbatch parameter: -c 1
test or run: set -e; useMemTimeWithInput.sh bigText1.txt; grep 1234 bigText1.txt > 1234.1.txt
depend on no job

Check if there input file list and this job does not depend on other jobs
inputSize: 1465
Running: estimateMemTime.sh useMemTimeWithInput none 1465
Estimating mem:
Finala: 0.001952218430034 Finalb: 0.660000000000000 Maximum: 7325.000000000000000
mem formula: ( 0.001952218430034 x 1465 + 0.660000000000000 ) x 1.0

Estimating time:
Finala: 0.001023890784415 Finalb: -0.699999996948272 Maximum: 7325.000000000000000
time formula: ( 0.001023890784415 x 1465 + -0.699999996948272 ) x 1.0
Got 3.519999999999810 .800000002219703
Got estimation inputsize: 1465 mem: 9M time: 6

0 0, 0 hour, 6 min, 0 sec

Building new sbatch command ...
New slurmScirpt is ready. The content is:
#!/bin/bash
trap "{ cleanUp.sh \"/home/ld32/smartSlurm\" "useMemTimeWithInput" "none" \"1.0.useMemTimeWithInput.1\" "1465" "1" "2G" "2:0:0" "9M" "0-0:6:0" "short"  \"\" \"/home/ld32/smartSlurm/bin/smartSbatch -L /home/ld32/smartSlurm -S useMemTimeWithInput -R none -F 1.0.useMemTimeWithInput.1 -I
,bigText1.txt -D null -p short -c 1 --mem 2G -t 2:0:0 --wrap set -e; useMemTimeWithInput.sh bigText1.txt; grep 1234 bigText1.txt > 1234.1.txt\"; }" EXIT
srun -n 1 bash -e -c "{ set -e; useMemTimeWithInput.sh bigText1.txt; grep 1234 bigText1.txt > 1234.1.txt; } && touch /home/ld32/smartSlurm/smartSlurmLog/1.0.useMemTimeWithInput.1.success"
New sbatch command to submit job:
/usr/bin/sbatch --mail-type=FAIL --requeue --parsable -p short --mem 9M -t 0-0:6:0 --open-mode=append -o /home/ld32/smartSlurm/smartSlurmLog/1.0.useMemTimeWithInput.1.out
 -e /home/ld32/smartSlurm/smartSlurmLog/1.0.useMemTimeWithInput.1.err -J 1.0.useMemTimeWithInput.1 -c 1 /home/ld32/smartSlurm/smartSlurmLog/1.0.useMemTimeWithInput.1.sh
This is a testing, not really running a job...

step: 2, depends on: 0, job name: useMemTimeWithInput, flag: useMemTimeWithInput.1
Running:
ssbatch -L /home/ld32/smartSlurm -S useMemTimeWithInput -R none -F 2.0.useMemTimeWithInput.1 -I ,bigText1.txt -D null -p short -c 1 --me
m 2G -t 2:0:0 --wrap "set -e; useMemTimeWithInput.sh bigText1.txt; grep 5678 bigText1.txt > 5678.1.txt"

Parsing result from sbatch commandline:
sbatch options: partition: short time: 2:0:0 mem: 2G mem-per-cpu: task: core: 1 node: out: err: dep:
wrapCMD: set -e; useMemTimeWithInput.sh bigText1.txt; grep 5678 bigText1.txt > 5678.1.txt
additional sbatch parameter: -c 1 -A 
test or run: set -e; useMemTimeWithInput.sh bigText1.txt; grep 5678 bigText1.txt > 5678.1.txt
depend on no job

Check if there input file list and this job does not depend on other jobs
inputSize: 1465
Running: estimateMemTime.sh useMemTimeWithInput none 1465
Estimating mem:
Finala: 0.001952218430034 Finalb: 0.660000000000000 Maximum: 7325.000000000000000
mem formula: ( 0.001952218430034 x 1465 + 0.660000000000000 ) x 1.0
Estimating time:
Finala: 0.001023890784415 Finalb: -0.699999996948272 Maximum: 7325.000000000000000
time formula: ( 0.001023890784415 x 1465 + -0.699999996948272 ) x 1.0
Got 3.519999999999810 .800000002219703
Got estimation inputsize: 1465 mem: 9M time: 6

0 0, 0 hour, 6 min, 0 sec

Building new sbatch command ...
New slurmScirpt is ready. The content is:
#!/bin/bash
trap "{ cleanUp.sh \"/home/ld32/smartSlurm\" "useMemTimeWithInput" "none" \"2.0.useMemTimeWithInput.1\" "1465" "1" "2G" "2:0:0" "9M" "0-0:6:0" "s
hort"  \"\" \"/home/ld32/smartSlurm/bin/smartSbatch -L /home/ld32/smartSlurm -S useMemTimeWithInput -R none -F 2.0.useMemTimeWithInput.1 -I
,bigText1.txt -D null -p short -c 1 --mem 2G -t 2:0:0 --wrap set -e; useMemTimeWithInput.sh bigText1.txt; grep 5678 bigText1.txt > 5678.1.txt\"; }" EXIT
srun -n 1 bash -e -c "{ set -e; useMemTimeWithInput.sh bigText1.txt; grep 5678 bigText1.txt > 5678.1.txt; } && touch /home/ld32/smartSlurm/smartSlurmLog/2.0.useM
emTimeWithInput.1.success"

New sbatch command to submit job:
/usr/bin/sbatch --mail-type=FAIL --requeue --parsable -p short --mem 9M -t 0-0:6:0 --open-mode=append -o /home/ld32/smartSlurm/smartSlurmLog/2.0.useMemTimeWithInput.1.out
 -e /home/ld32/smartSlurm/smartSlurmLog/2.0.useMemTimeWithInput.1.err -J 2.0.useMemTimeWithInput.1 -c 1 /home/ld32/smartSlurm/smartSlurmLog/2.0.useMemTimeWithInput.1.sh
This is a testing, not really running a job...

step: 3, depends on: 1.2, job name: useMemTimeWithInput, flag: useMemTimeWithInput
Running:
ssbatch -L /home/ld32/smartSlurm -S useMemTimeWithInput -R none -F 3.1.2.useMemTimeWithInput -I ,bigText1.txt -D ..123..123 -p s
hort -t 10:0 -c 1 --wrap "set -e; useMemTimeWithInput.sh bigText1.txt; cat 1234.1.txt 1234.2.txt 5678.1.txt 5678.2.txt > all.txt"

Parsing result from sbatch commandline:
sbatch options: partition: short time: 10:0 mem: mem-per-cpu: task: core: 1 node: out: err: dep:
wrapCMD: set -e; useMemTimeWithInput.sh bigText1.txt; cat 1234.1.txt 1234.2.txt 5678.1.txt 5678.2.txt > all.txt
additional sbatch parameter: -c 1
test or run: set -e; useMemTimeWithInput.sh bigText1.txt; cat 1234.1.txt 1234.2.txt 5678.1.txt 5678.2.txt > all.txt
depend on multiple jobs
working on 123
working on 123

Check if there input file list and this job does not depend on other jobs
inputSize: 1465
Running: estimateMemTime.sh useMemTimeWithInput none 1465
Estimating mem:
Finala: 0.001952218430034 Finalb: 0.660000000000000 Maximum: 7325.000000000000000
mem formula: ( 0.001952218430034 x 1465 + 0.660000000000000 ) x 1.0

Estimating time:
Finala: 0.001023890784415 Finalb: -0.699999996948272 Maximum: 7325.000000000000000
time formula: ( 0.001023890784415 x 1465 + -0.699999996948272 ) x 1.0
Got 3.519999999999810 .800000002219703
Got estimation inputsize: 1465 mem: 9M time: 6

0 0, 0 hour, 6 min, 0 sec

Building new sbatch command ...
New slurmScirpt is ready. The content is:
#!/bin/bash
trap "{ cleanUp.sh \"/home/ld32/smartSlurm\" "useMemTimeWithInput" "none" \"3.1.2.useMemTimeWithInput\" "1465" "1" "2G" "10:0" "9M" "0-0:6:0" "short"  \"\" \"/home/ld32/smartSlurm/bin/smartSbatch -L /home/ld32/smartSlurm -S useMemTimeWithInput -R none -F 3.1.2.useMemTimeWithInput -I ,bigText1.txt -D ..123..123 -p short -t 10:0 -c 1 --wrap set -e; useMemTimeWithInput.sh bigText1.txt; cat 1234.*.txt 5678.*.txt > all.txt\"; }" EXIT
srun -n 1 bash -e -c "{ set -e; useMemTimeWithInput.sh bigText1.txt; cat 1234.1.txt 1234.2.txt 5678.1.txt 5678.2.txt > all.txt; } && touch /home/ld32/smartSlurm/smartSlurmLog/3.1.2.useMemTimeWithInput.success"

New sbatch command to submit job:
/usr/bin/sbatch --mail-type=FAIL --requeue --parsable -p short --mem 9M -t 0-0:6:0 --open-mode=append -o /home/ld32/smartSlurm/smartSlurmLog/3.1.2.useMemTimeWithInput.out -e /home/ld32/smartSlurm/smartSlurmLog/3.1.2.useMemTimeWithInput.err -J 3.1.2.useMemTimeWithInput --dependency=afterok:123:123 -c 1 /home/ld32/smartSlurm/smartSlurmLog/3.1.2.useMemTimeWithInput
This is a testing, not really running a job...

All submitted jobs:
job_id       depend_on              job_flag
123         null                  1.0.useMemTimeWithInput.1
124         null                  2.0.useMemTimeWithInput.1
125         ..123..124            3.1.2.useMemTimeWithInput
---------------------------------------------------------
Note: This is just a test run, so no job is actually submitted. In real run it should submit jobs and report as above.


Run the modified bash script as a pipeline

Thus far in the example, we have not actually submitted any jobs to the scheduler. To submit the pipeline, you will need to append the run parameter to the command. If run is not specified, test mode will be used, which does not submit jobs and gives the placeholder of 1234for jobids in the command's output. 
runAsPipeline ~/smartSbatch/scripts/bashScriptV2.sh "sbatch -p short -t 10:0 -c 1" useTmp run

# Below is the output
[Back to top](#SmartSlurm)


Wed Dec 21 16:02:47 EST 2022
Running: /home/ld32/smartSlurm/bin/runAsPipeline /home/ld32/smartSlurm/scripts/bashScriptV2.sh sbatch -p short -t 10:0 -c 1 noTmp run

Currently Loaded Modules:
  1) gcc/6.2.0


This is a re-run with the same command and script is not changed, no need to convert the script. Using the old one: $smartSlurmLogDir/slurmPipeLine.6f93dc8953b9c1d1f96b4fabd657446a.run.sh
Running $smartSlurmLogDir/slurmPipeLine.6f93dc8953b9c1d1f96b4fabd657446a.run.sh /home/ld32/smartSlurm/scripts/bashScriptV2.sh

Currently Loaded Modules:
  1) gcc/6.2.0

---------------------------------------------------------

step: 1, depends on: 0, job name: useMemTimeWithInput, flag: useMemTimeWithInput.1
Running:
ssbatch -L /home/ld32/smartSlurm -S useMemTimeWithInput -R none -F 1.0.useMemTimeWithInput.1 -I ,bigText1.txt -D null -p short -c 1 --mem 2G -t 2:0:0 --wrap "set -e; useMemTimeWithInput bigText1.txt; grep 1234 bigText1.txt > 1234.1.txt" run

sbatch options: partition: short time: 2:0:0 mem: 2G mem-per-cpu: task: core: 1 node: out: err: dep:
wrapCMD: set -e; useMemTimeWithInput.sh bigText1.txt; grep 1234 bigText1.txt > 1234.1.txt
additional sbatch parameter: -c 1 -A 
test or run: run
depend on no job

Check if there input file list and this job does not depend on other jobs
inputSize: 1465
Running: estimateMemTime.sh useMemTimeWithInput none 1465
Estimating mem:
Finala: 0.001952218430034 Finalb: 0.660000000000000 Maximum: 7325.000000000000000
mem formula: ( 0.001952218430034 x 1465 + 0.660000000000000 ) x 1.0

Estimating time:
Finala: 0.001023890784415 Finalb: -0.699999996948272 Maximum: 7325.000000000000000
time formula: ( 0.001023890784415 x 1465 + -0.699999996948272 ) x 1.0
Got 3.519999999999810 .800000002219703
Got estimation inputsize: 1465 mem: 9M time: 6

0 0, 0 hour, 6 min, 0 sec

Building new sbatch command ...
New slurmScirpt is ready. The content is:
#!/bin/bash
trap "{ cleanUp.sh \"/home/ld32/smartSlurm\" "useMemTimeWithInput" "none" \"1.0.useMemTimeWithInput.1\" "1465" "1" "2G" "2:0:0" "9M" "0-0:6:0" "short"  \"\" \"/home/ld32/smartSlurm/bin/smartSbatch -L /home/ld32/smartSlurm -S useMemTimeWithInput -R none -F 1.0.useMemTimeWithInput.1 -I ,bigText1.txt -D null -p short -c 1 --mem 2G -t 2:0:0 --wrap set -e; useMemTimeWithInput bigText1.txt; grep 1234 bigText1.txt > 1234.1.txt run\"; }" EXIT
srun -n 1 bash -e -c "{ set -e; useMemTimeWithInput.sh bigText1.txt; grep 1234 bigText1.txt > 1234.1.txt; } && touch /home/ld32/smartSlurm/smartSlurmLog/1.0.useMemTimeWithInput.1.success"

New sbatch command to submit job:
/usr/bin/sbatch --mail-type=FAIL --requeue --parsable -p short --mem 9M -t 0-0:6:0 --open-mode=append -o /home/ld32/smartSlurm/smartSlurmLog/1.0.useMemTimeWithInput.1.out -e /home/ld32/smartSlurm/smartSlurmLog/1.0.useMemTimeWithInput.1.err -J 1.0.useMemTimeWithInput.1 -c 1 /home/ld32/smartSlurm/smartSlurmLog/1.0.useMemTimeWithInput.1.sh
Start submtting job...

step: 2, depends on: 0, job name: useMemTimeWithInput, flag: useMemTimeWithInput.1
Running:
ssbatch -L /home/ld32/smartSlurm -S useMemTimeWithInput -R none -F 2.0.useMemTimeWithInput.1 -I ,bigText1.txt -D null -p short -c 1 --mem 2G -t 2:0:0 --wrap "set -e; useMemTimeWithInput bigText1.txt; grep 5678 bigText1.txt > 5678.1.txt" run

Parsing result from sbatch commandline:
sbatch options: partition: short time: 2:0:0 mem: 2G mem-per-cpu: task: core: 1 node: out: err: dep:
wrapCMD: set -e; useMemTimeWithInput.sh bigText1.txt; grep 5678 bigText1.txt > 5678.1.txt
additional sbatch parameter: -c 1 -A 
test or run: run
depend on no job

Check if there input file list and this job does not depend on other jobs
inputSize: 1465
Running: estimateMemTime.sh useMemTimeWithInput none 1465
Estimating mem:
Finala: 0.001952218430034 Finalb: 0.660000000000000 Maximum: 7325.000000000000000
mem formula: ( 0.001952218430034 x 1465 + 0.660000000000000 ) x 1.0

Estimating time:
Finala: 0.001023890784415 Finalb: -0.699999996948272 Maximum: 7325.000000000000000
time formula: ( 0.001023890784415 x 1465 + -0.699999996948272 ) x 1.0
Got 3.519999999999810 .800000002219703
Got estimation inputsize: 1465 mem: 9M time: 6

0 0, 0 hour, 6 min, 0 sec

Building new sbatch command ...
New slurmScirpt is ready. The content is:
#!/bin/bash
trap "{ cleanUp.sh \"/home/ld32/smartSlurm\" "useMemTimeWithInput" "none" \"2.0.useMemTimeWithInput.1\" "1465" "1" "2G" "2:0:0" "9M" "0-0:6:0" "short"  \"\" \"/home/ld32/smartSlurm/bin/smartSbatch -L /home/ld32/smartSlurm -S useMemTimeWithInput -R none -F 2.0.useMemTimeWithInput.1 -I ,bigText1.txt -D null -p short -c 1 --mem 2G -t 2:0:0 --wrap set -e; useMemTimeWithInput.sh bigText1.txt; grep 5678 bigText1.txt > 5678.1.txt run\"; }" EXIT
srun -n 1 bash -e -c "{ set -e; useMemTimeWithInput.sh bigText1.txt; grep 5678 bigText1.txt > 5678.1.txt; } && touch /home/ld32/smartSlurm/smartSlurmLog/2.0.useMemTimeWithInput.1.success"

New sbatch command to submit job:
/usr/bin/sbatch --mail-type=FAIL --requeue --parsable -p short --mem 9M -t 0-0:6:0 --open-mode=append -o /home/ld32/smartSlurm/smartSlurmLog/2.0.useMemTimeWithInput.1.out -e /home/ld32/smartSlurm/smartSlurmLog/2.0.useMemTimeWithInput.1.err -J 2.0.useMemTimeWithInput.1 -c 1 /home/ld32/smartSlurm/smartSlurmLog/2.0.useMemTimeWithInput.1.sh
Start submtting job...

step: 3, depends on: 1.2, job name: useMemTimeWithInput, flag: useMemTimeWithInput
Running:

ssbatch -L /home/ld32/smartSlurm -S useMemTimeWithInput -R none -F 3.1.2.useMemTimeWithInput -I ,bigText1.txt -D ..46631..46632 -p short -t 10:0 -c 1 --wrap "set -e; useMemTimeWithInput bigText1.txt; cat 1234.1.txt 1234.2.txt 5678.1.txt 5678.2.txt > all.txt" run

Parsing result from sbatch commandline:
sbatch options: partition: short time: 10:0 mem: mem-per-cpu: task: core: 1 node: out: err: dep:
wrapCMD: set -e; useMemTimeWithInput.sh bigText1.txt; cat 1234.1.txt 1234.2.txt 5678.1.txt 5678.2.txt > all.txt
additional sbatch parameter: -c 1
test or run: run
depend on multiple jobs
working on 46631
working on 46632

Check if there input file list and this job does not depend on other jobs
inputSize: 1465
Running: estimateMemTime.sh useMemTimeWithInput none 1465
Estimating mem:
Finala: 0.001952218430034 Finalb: 0.660000000000000 Maximum: 7325.000000000000000
mem formula: ( 0.001952218430034 x 1465 + 0.660000000000000 ) x 1.0

Estimating time:
Finala: 0.001023890784415 Finalb: -0.699999996948272 Maximum: 7325.000000000000000
time formula: ( 0.001023890784415 x 1465 + -0.699999996948272 ) x 1.0
Got 3.519999999999810 .800000002219703
Got estimation inputsize: 1465 mem: 9M time: 6

0 0, 0 hour, 6 min, 0 sec

Building new sbatch command ...

New slurmScirpt is ready. The content is:
#!/bin/bash
trap "{ cleanUp.sh \"/home/ld32/smartSlurm\" "useMemTimeWithInput" "none" \"3.1.2.useMemTimeWithInput\" "1465" "1" "2G" "10:0" "9M" "0-0:6:0" "short"  \"\" \"/home/ld32/smartSlurm/bin/smartSbatch -L /home/ld32/smartSlurm -S useMemTimeWithInput -R none -F 3.1.2.useMemTimeWithInput -I ,bigText1.txt -D ..46631..46632 -p short -t 10:0 -c 1 --wrap set -e; useMemTimeWithInput.sh bigText1.txt; cat 1234.*.txt 5678.*.txt > all.txt run\"; }" EXIT
srun -n 1 bash -e -c "{ set -e; useMemTimeWithInput.sh bigText1.txt; cat 1234.1.txt 1234.2.txt 5678.1.txt 5678.2.txt > all.txt; } && touch /home/ld32/smartSlurm/smartSlurmLog/3.1.2.useMemTimeWithInput.success"

New sbatch command to submit job:
/usr/bin/sbatch --mail-type=FAIL --requeue --parsable -p short --mem 9M -t 0-0:6:0 --open-mode=append -o /home/ld32/smartSlurm/smartSlurmLog/3.1.2.useMemTimeWithInput.out -e /home/ld32/smartSlurm/smartSlurmLog/3.1.2.useMemTimeWithInput.err -J 3.1.2.useMemTimeWithInput --dependency=afterok:46631:46632 -c 1 /home/ld32/smartSlurm/smartSlurmLog/3.1.2.useMemTimeWithInput
Start submtting job...

All submitted jobs:
job_id       depend_on              job_flag
46631       null                  1.0.useMemTimeWithInput.1
46632       null                  2.0.useMemTimeWithInput.1
46633       ..46631..46632        3.1.2.useMemTimeWithInput

Monitoring the jobs

You can use the command:
squeue -u $USER --Format=jobid:10,username:6,partition:14,name:35,state:14,timeused:10,timeleft:10,timelimit:180,starttime:18,nodelist:18,numcpus:8,minmemory:30

To see the job status (running, pending, etc.). You also get two emails for each step, one at the start of the step, one at the end of the step.
Successful job email
Email subject: Success: job id:46631 name:1.0.useMemTimeWithInput.1

Email content:

Job script content:
#!/bin/bash
trap "{ cleanUp.sh \"/home/ld32/smartSlurm\" "useMemTimeWithInput" "none" \"1.0.useMemTimeWithInput.1\" "1465" "1" "2G" "2:0:0" "9M" "0-0:6:0" "short"  \"\" \"/home/ld32/smartSlurm/bin/smartSbatch -L /home/ld32/smartSlurm -S useMemTimeWithInput -R none -F 1.0.useMemTimeWithInput.1 -I ,bigText1.txt -D null -p short -c 1 --mem 2G -t 2:0:0 --wrap set -e; useMemTimeWithInput.sh bigText1.txt; grep 1234 bigText1.txt > 1234.1.txt run\"; }" EXIT
srun -n 1 bash -e -c "{ set -e; useMemTimeWithInput.sh bigText1.txt; grep 1234 bigText1.txt > 1234.1.txt; } && touch /home/ld32/smartSlurm/smartSlurmLog/1.0.useMemTimeWithInput.1.success"

#Command used to submit the job:
#/usr/bin/sbatch --mail-type=FAIL --requeue --parsable -p short --mem 9M -t 0-0:6:0 --open-mode=append -o /home/ld32/smartSlurm/smartSlurmLog/1.0.useMemTimeWithInput.1.out -e /home/ld32/smartSlurm/smartSlurmLog/1.0.useMemTimeWithInput.1.err -J 1.0.useMemTimeWithInput.1     -c 1    /home/ld32/smartSlurm/smartSlurmLog/1.0.useMemTimeWithInput.1.sh

#Sbatch command output:
#Submitted batch job 46631

Job output:
Begin allocating memory...
...end allocating memory. Begin sleeping for 60 seconds...
Done
Running /home/ld32/smartSlurm/bin/cleanUp.sh /home/ld32/smartSlurm useMemTimeWithInput none 1.0.useMemTimeWithInput.1 1465 1 2G 2:0:0 9M 0-0:6:0 short /home/ld32/smartSlurm/bin/smartSbatch -L /home/ld32/smartSlurm -S useMemTimeWithInput -R none -F 1.0.useMemTimeWithInput.1 -I ,bigText1.txt -D null -p short -c 1 --mem 2G -t 2:0:0 --wrap set -e; useMemTimeWithInput.sh bigText1.txt; grep 1234 bigText1.txt > 1234.1.txt run

Job summary:
JobID                     Submit               Start                 End     MaxRSS      State                       NodeList  Partition                        ReqTRES   TotalCPU        Elapsed      Timelimit
------------ ------------------- ------------------- ------------------- ---------- ---------- ------------------------------ ---------- ------------------------------ ---------- -------------- --------------
46631        2022-12-21T16:03:24 2022-12-21T16:03:24             Unknown               RUNNING                  compute-x      short  billing=1,cpu=1,mem=9M,node=1  00:00.076       00:01:06       00:06:00
46631.batch  2022-12-21T16:03:24 2022-12-21T16:03:24             Unknown               RUNNING                  compute-x                                             00:00:00       00:01:06               
46631.extern 2022-12-21T16:03:24 2022-12-21T16:03:24             Unknown               RUNNING                  compute-x                                             00:00:00       00:01:06               
46631.0      2022-12-21T16:03:25 2022-12-21T16:03:25 2022-12-21T16:04:25      3.49M  COMPLETED                  compute-x                                            00:00.076       00:01:00               
*Notice the sacct report above: while the main job is still running for sacct command, user task is completed.
Last row of job summary: 46631.0      2022-12-21T16:03:25 2022-12-21T16:03:25 2022-12-21T16:04:25      3.49M  COMPLETED                  compute-x                                            00:00.076       00:01:00               
start: 1671656605 finish: 1671656665 mem: 3.49M mins: 1
jobStatus: COMPLETED
Added this line to $smartSlurmJobRecordDir/myJobRecord.txt:
46631,1465,2G,2:0:0,9M,0-0:6:0,3.49,1,COMPLETED,ld32,/home/ld32/smartSlurm,useMemTimeWithInput,none,1.0.useMemTimeWithInput.1,1,compute-x,/home/ld32/smartSlurm/smartSlurmLog/1.0.useMemTimeWithInput.1.err,Wed Dec 21 16:04:30 EST 2022,"/home/ld32/smartSlurm/bin/smartSbatch -L /home/ld32/smartSlurm -S useMemTimeWithInput -R none -F 1.0.useMemTimeWithInput.1 -I ,bigText1.txt -D null -p short -c 1 --mem 2G -t 2:0:0 --wrap set -e; useMemTimeWithInput.sh bigText1.txt; grep 1234 bigText1.txt > 1234.1.txt run"
Running: /home/ld32/smartSlurm/bin/adjustDownStreamJobs.sh /home/ld32/smartSlurm/log
Find current job id (flag: 1.0.useMemTimeWithInput.1):
46631

Find all downstream jobs which depend on current job
job idNames:
..46631..46632 3.1.2.useMemTimeWithInput
1working on ..46631..46632 3.1.2.useMemTimeWithInput
2working on 46631
Ignore. It is the current job. It should adjust the mem and time for the downsteam job.
2working on 46632
look for the job flag for 46632
This job was done!
Dependants for 3.1.2.useMemTimeWithInput are all done except for the current job. Ready to adjust mem/runtime
Do not have a formula. Let us build one...
Running
/home/ld32/smartSlurm/bin/jobStatistics.sh 4
Usage: checkJobRecord.sh software reference
Try to build fomular, but it was not successful

Error output:
The key elements are time and memory used.

# Check job logs

You can use the command:
ls -l log

This command list all the logs created by the pipeline runner. *.sh files are the slurm scripts for each step, *.out files are output files for each step, *.success files means job successfully finished for each step and *.failed means job failed for each steps.

You also get two emails for each step, one at the start of the step, one at the end of the step.
Cancel all jobs

You can use the command to cancel running and pending jobs:
cancelAllJobs $smartSlurmLogDir/alljobs.jid
What happens if there is some error? 

You can re-run this command in the same folder. We will delete an input file to see what happens.
# We are intentionally removing an input file to see a "failed job" email message
rm universityB.txt
runAsPipeline bashScriptV2.sh "sbatch -p short -t 10:0 -c 1" useTmp run

# Here is the output
Fri Sep 24 10:00:36 EDT 2021
Running: /home/ld32/smartSlurm/bin/runAsPipeline bashScriptV2.sh sbatch -p short -t 10:0 -c 1 useTmp run

Currently Loaded Modules:
This is a re-run with the same command and script is not changed, no need to convert the script. Using the old one: $smartSlurmLogDir/slurmPipeLine.a855454a70b2198fa5b2643bb1d41762.run.sh
Running $smartSlurmLogDir/slurmPipeLine.a855454a70b2198fa5b2643bb1d41762.run.sh bashScriptV2.sh


Could not find any jobs to cancel.
---------------------------------------------------------

step: 1, depends on: 0, job name: find1, flag: find1.A reference: .u
depend on no job
1.0.find1.A was done before, do you want to re-run it?
y:        To re-run this job, press y, then enter key.
ystep:    To re-run all jobs for step 3: hisatCount, type yall, then press enter key.
yall:     To re-run all jobs, type yallall, then press enter key.
enter:    To not re-run this job, directly press enter key.
nstep:    To not re-run all successful jobs for step 3: hisatCount, type nall, then press enter key.
nall:     To not re-run all successful jobs, type nallall, then press enter key.

\# type enter here to not re-run

step: 2, depends on: 0, job name: find2, flag: find2.A reference: .u
depend on no job
2.0.find2.A was done before, do you want to re-run it?
y:        To re-run this job, press y, then enter key.
ystep:    To re-run all jobs for step 3: hisatCount, type yall, then press enter key.
yall:     To re-run all jobs, type yallall, then press enter key.
enter:    To not re-run this job, directly press enter key.
nstep:    To not re-run all successful jobs for step 3: hisatCount, type nall, then press enter key.
nall:     To not re-run all successful jobs, type nallall, then press enter key.

\# type enter here to not re-run

job 2.0.find2.A is not submitted

step: 1, depends on: 0, job name: find1, flag: find1.B reference: .u
depend on no job
1.0.find1.B was done before, do you want to re-run it?
y:        To re-run this job, press y, then enter key.
ystep:    To re-run all jobs for step 3: hisatCount, type yall, then press enter key.
yall:     To re-run all jobs, type yallall, then press enter key.
enter:    To not re-run this job, directly press enter key.
nstep:    To not re-run all successful jobs for step 3: hisatCount, type nall, then press enter key.
nall:     To not re-run all successful jobs, type nallall, then press enter key.

\# type ‘y’ and enter here to re-run

Will re-run the down stream steps even if they are done before.
sbatch -p short -c 1 -t 2:0:0 --requeue --nodes=1  -J 1.0.find1.B -o /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.B.out -e /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.B.out /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.B.sh
\# Submitted batch job 41209197

step: 2, depends on: 0, job name: find2, flag: find2.B reference: .u
depend on no job
2.0.find2.B was done before, do you want to re-run it?
y:        To re-run this job, press y, then enter key.
ystep:    To re-run all jobs for step 3: hisatCount, type yall, then press enter key.
yall:     To re-run all jobs, type yallall, then press enter key.
enter:    To not re-run this job, directly press enter key.
nstep:    To not re-run all successful jobs for step 3: hisatCount, type nall, then press enter key.
nall:     To not re-run all successful jobs, type nallall, then press enter key.

\# type enter here to not re-run

job 2.0.find2.B is not submitted

step: 3, depends on: 1.2, job name: merge , flag: merge reference:
depend on other jobs
sbatch -p short -t 10:0 -c 1 --requeue --nodes=1 --dependency=afterok:41209197 -J 3.1.2.merge -o /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/3.1.2.merge.out -e /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/3.1.2.merge.out /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/3.1.2.merge.sh
\# Submitted batch job 41209210

\# Notice above, rcbio didn’t ask if user wants to re-run step3 or not and directly re-run it.

All submitted jobs:
job_id       depend_on              job_flag
41209197    null                  1.0.find1.B
41209210    ..41209197.           3.1.2.merge
---------------------------------------------------------

This command will check if the earlier run is finished or not. If not, ask user to kill the running jobs or not, then ask user to rerun the successfully finished steps or not. Click 'y', it will rerun, directly press 'enter' key, it will not rerun. 
Failed job email
Email subject: Failed: job id:41209197 name:1.0.find1.B

Email content:
Job script content:
#!/bin/bash
#Commands:
trap "{ cleanup.sh /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.B; }” EXIT
touch /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.B.start
srun -n 1 bash -e -c "{ set -e; rsyncToTmp  /tmp/rcbio/universityB.txt; grep -H John /tmp/rcbio/universityB.txt >>  John.txt; grep -H Mike /tmp/rcbio/universityB.txt >>  Mike.txt        ; } && touch /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.B.success || touch /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.B.failed"

#sbatch command:
#sbatch -p short -c 1 -t 2:0:0 --requeue --nodes=1  -J 1.0.find1.B -o /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.B.out -e /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.B.out /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.B.sh

\# Submitted batch job 41209197
Job output:
Working to copy: /tmp/rcbio/universityB.txt, waiting lock...
Reference file or folder not exist: /universityB.txt
grep: /tmp/rcbio/universityB.txt: No such file or directory
grep: /tmp/rcbio/universityB.txt: No such file or directory
Job done. Summary:
       JobID              Submit               Start                 End      State  Partition              ReqTRES  Timelimit    CPUTime     MaxRSS                       NodeList
------------ ------------------- ------------------- ------------------- ---------- ---------- -------------------- ---------- ---------- ---------- ------------------------------
41209197     2021-09-24T10:02:43 2021-09-24T10:03:09             Unknown    RUNNING      short billing=1,cpu=1,mem+   00:2:0:0   00:00:09                          compute-x
41209197.ba+ 2021-09-24T10:03:09 2021-09-24T10:03:09             Unknown    RUNNING                                              00:00:09                          compute-x
41209197.ex+ 2021-09-24T10:03:09 2021-09-24T10:03:09             Unknown    RUNNING                                              00:00:09                          compute-x
41209197.0   2021-09-24T10:03:13 2021-09-24T10:03:13 2021-09-24T10:03:13  COMPLETED                                              00:00:00          0               compute-x
*Notice the sacct report above: while the main job is still running for sacct command, user task is completed.

    The key element here is the error message.

    Notice here, step2 job is automatically canceled because this job failed. We deleted universityB.txt, so the job has failed. We don’t get an email from the downstream step3 job. 

Fix the error and re-run the pipeline

You can rerun this command in the same folder
cp universityA.txt universityB.txt
runAsPipeline bashScriptV2.sh "sbatch -p short -t 10:0 -c 1" useTmp run

This command will automatically check if the earlier run is finished. If the run has not finished, the script will ask the user if they want to kill the running jobs or not, then ask user to rerun the successfully finished steps or not. Click 'y', it will rerun, directly press 'enter' key, it will not rerun. 

Notice here, step3 will run by default. It will run without prompting the user for permission.
What happens if we add more input data and re-run the pipeline?

You can rerun this command in the same folder
cp universityA.txt universityC.txt
cp bashScriptV2.sh bashScriptV3.sh 
nano bashScriptV3.sh  
\# change
for i in A B; do
to: 
for i in A B C; do

\# save the file and run:
runAsPipeline bashScriptV3.sh "sbatch -p short -t 10:0 -c 1" useTmp run

\# Here are the output: 
Fri Sep 24 10:56:16 EDT 2021
Running: /home/ld32/smartSlurm/bin/runAsPipeline bashScriptV3.sh sbatch -p short -t 10:0 -c 1 useTmp run

Currently Loaded Modules:
  

converting bashScriptV3.sh to $smartSlurmLogDir/slurmPipeLine.b72e7f91da30d312a2c85d0735896f79.run.sh

find loop start: for i in A B C; do

find job marker:
#@1,0,find1,u,,sbatch -p short -c 1 -t 2:0:0
sbatch options: sbatch -p short -c 1 -t 2:0:0

find job:
grep -H John $u >>  John.txt; grep -H Mike $u >>  Mike.txt

find job marker:
#@2,0,find2,u,,sbatch -p short -c 1 -t 2:0:0
sbatch options: sbatch -p short -c 1 -t 2:0:0

find job:
grep -H Nick $u >>  Nick.txt; grep -H Julia $u >>  Julia.txt
find loop end: done

find job marker:
#@3,1.2,merge

find job:
cat John.txt Mike.txt Nick.txt Julia.txt > all.txt
smartSlurmLog/slurmPipeLine.b72e7f91da30d312a2c85d0735896f79.run.sh bashScriptV3.sh is ready to run. Starting to run ...
Running $smartSlurmLogDir/slurmPipeLine.b72e7f91da30d312a2c85d0735896f79.run.sh bashScriptV3.sh

Currently Loaded Modules:
  


Could not find any jobs to cancel.
---------------------------------------------------------

step: 1, depends on: 0, job name: find1, flag: find1.A reference: .u
depend on no job
1.0.find1.A was done before, do you want to re-run it?
y:        To re-run this job, press y, then enter key.
ystep:    To re-run all jobs for step 3: hisatCount, type yall, then press enter key.
yall:     To re-run all jobs, type yallall, then press enter key.
enter:    To not re-run this job, directly press enter key.
nstep:    To not re-run all successful jobs for step 3: hisatCount, type nall, then press enter key.
nall:     To not re-run all successful jobs, type nallall, then press enter key.

\# type enter here to not re-run

job 1.0.find1.A is not submitted

step: 2, depends on: 0, job name: find2, flag: find2.A reference: .u
depend on no job
2.0.find2.A was done before, do you want to re-run it?
y:        To re-run this job, press y, then enter key.
ystep:    To re-run all jobs for step 3: hisatCount, type yall, then press enter key.
yall:     To re-run all jobs, type yallall, then press enter key.
enter:    To not re-run this job, directly press enter key.
nstep:    To not re-run all successful jobs for step 3: hisatCount, type nall, then press enter key.
nall:     To not re-run all successful jobs, type nallall, then press enter key.

\# type enter here to not re-run

job 2.0.find2.A is not submitted

step: 1, depends on: 0, job name: find1, flag: find1.B reference: .u
depend on no job
1.0.find1.B was done before, do you want to re-run it?
y:        To re-run this job, press y, then enter key.
ystep:    To re-run all jobs for step 3: hisatCount, type yall, then press enter key.
yall:     To re-run all jobs, type yallall, then press enter key.
enter:    To not re-run this job, directly press enter key.
nstep:    To not re-run all successful jobs for step 3: hisatCount, type nall, then press enter key.
nall:     To not re-run all successful jobs, type nallall, then press enter key.

\# type enter here to not re-run

job 1.0.find1.B is not submitted

step: 2, depends on: 0, job name: find2, flag: find2.B reference: .u
depend on no job
2.0.find2.B was done before, do you want to re-run it?
y:        To re-run this job, press y, then enter key.
ystep:    To re-run all jobs for step 3: hisatCount, type yall, then press enter key.
yall:     To re-run all jobs, type yallall, then press enter key.
enter:    To not re-run this job, directly press enter key.
nstep:    To not re-run all successful jobs for step 3: hisatCount, type nall, then press enter key.
nall:     To not re-run all successful jobs, type nallall, then press enter key.

\# type enter here to not re-run

job 2.0.find2.B is not submitted

step: 1, depends on: 0, job name: find1, flag: find1.C reference: .u
depend on no job
sbatch -p short -c 1 -t 2:0:0 --requeue --nodes=1  -J 1.0.find1.C -o /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.C.out -e /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.C.out /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.C.sh
\# Submitted batch job 41211380

step: 2, depends on: 0, job name: find2, flag: find2.C reference: .u
depend on no job
sbatch -p short -c 1 -t 2:0:0 --requeue --nodes=1  -J 2.0.find2.C -o /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/2.0.find2.C.out -e /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/2.0.find2.C.out /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/2.0.find2.C.sh
\# Submitted batch job 41211381

step: 3, depends on: 1.2, job name: merge , flag: merge reference:
depend on multiple jobs
sbatch -p short -t 10:0 -c 1 --requeue --nodes=1 --dependency=afterok:41211380:41211381 -J 3.1.2.merge -o /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/3.1.2.merge.out -e /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/3.1.2.merge.out /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/3.1.2.merge.sh
\# Submitted batch job 41211382

All submitted jobs:
job_id       depend_on              job_flag
41211380    null                  1.0.find1.C
41211381    null                  2.0.find2.C
41211382    ..41211380..41211381  3.1.2.merge
---------------------------------------------------------

This command will check if the earlier run is finished, and will prompt the user if they kill any running jobs. Next, it will then ask the user if they want to rerun any successfully finished steps. Click 'y', it will rerun, directly press 'enter' key, it will not rerun. 

For the new data, RCBio will submit 2 jobs. Step3 will also still automatically run.
Re-run a single job manually
\# /working/directory is a placeholder, replace it with your actual working directory path
cd /working/directory
\# all/related/modules is a placeholder, replace it with the actual other modules/versions you need
module load rcbio/1.3.3 and all/related/modules

\# submit job with proper partition, time, number of cores and memory
sbatch --requeue --mail-type=ALL -p short -t 2:0:0 -c 2 --mem 2G /working/directory/smartSlurmLog/stepID.loopID.stepName.sh

Or:
runSingleJob "module load bowtie/1.2.2; bowtie -x /n/groups/shared_databases/bowtie_indexes/hg19 -p 2 -1 read1.fq -2 read2.fq --sam > out.bam" "sbatch -p short -t 1:0:0 -c 2 -mem 8G"

For details about the second option: Get more informative slurm email notification and logs through rcbio/1.3 
To run your own script as Slurm pipeline

If you have a bash script with multiple steps and you wish to run it as Slurm pipeline, here is how you can do that:

    modify your old script and add the notation to mark the start and end of any loops, and the start of any step for which you want to submit as an sbatch job. 

    use runAsPipeline with your modified bash script, as detailed above. 

How does the runAsPipeline RCBio pipeline runner work?

In case you wonder how it works, here is a simple example to explain.

For each step per loop, the pipeline runner creates a file that looks like the one below. (Here it is named flag.sh): 
#!/bin/bash 
srun -n 1 bash -c "{ echo I am running...; hostname; otherCommands; } && touch flag.success" 
sleep 5 
export SLURM_TIME_FORMAT=relative 
echo Job done. Summary: 
sacct --format=JobID,Submit,Start,End,State,Partition,ReqTRES%30,CPUTime,MaxRSS,NodeList%30 --units=M -j $SLURM_JOBID 
sendJobFinishEmail.sh flag 
[ -f flag.success ] && exit 0 || exit 1 

Your analysis commands will be wrapped in an srun so we can monitor if it completed successfully. If your commands worked (meaning exited in 0 status), then we will create the success file. Next, we will run sacct to get stats for the job step, and will send a job completion email with sendJobFinishEmail.sh. The sendJobFinishEmail.sh script is available in /home/ld32/smartSlurm/bin/, if you are interested in the contents of that script.

Then the job script will be submitted with: 
sbatch -p short -t 10:0 -o flag.out -e flag.out flag.sh

 

Let us know if you have any questions by emailing rchelp@hms.harvard.edu. Please include your working folder and the commands used in your email. Any comments and suggestions are welcome!

We have additional example ready-to-run workflows available, which may be of interest to you.Test run the modified bash script as a pipeline
runAsPipeline bashScriptV2.sh "sbatch -p short -t 10:0 -c 1" useTmp

This command will generate new bash script of the form slurmPipeLine.checksum.sh in flag folder. The checksum portion of the filename will have a MD5 hash that represents the file contents. We include the checksum in the filename to detect when script contents have been updated.

This runAsPipeline command will run a test of the script, meaning does not really submit jobs. It will only show a fake job ids like 1234 for each step. If you were to append run at the end of the command, the pipeline would actually be submitted to the Slurm scheduler.

Ideally, with useTmp, the software should run faster using local /tmp disk space for database/reference than the network storage. For this small query, the difference is small, or even slower if you use local /tmp. If you don't need /tmp, you can use noTmp.

With useTmp, the pipeline runner copy related data to /tmp and all file paths will be automatically updated to reflect a file's location in /tmp when using the useTmp option. 
Sample output from the test run

Note that only step 2 used -t 2:0:0, and all other steps used the default -t 10:0. The default walltime limit was set in the runAsPipeline command, and the walltime parameter for step 2 was set in the bash_script_v2.sh script.
runAsPipeline bashScriptV2.sh "sbatch -p short -t 10:0 -c 1" useTmp

Fri Sep 24 09:46:15 EDT 2021
Running: /home/ld32/smartSlurm/bin/runAsPipeline bashScriptV2.sh sbatch -p short -t 10:0 -c 1 useTmp

Currently Loaded Modules:
  

converting bashScriptV2.sh to $smartSlurmLogDir/slurmPipeLine.a855454a70b2198fa5b2643bb1d41762.sh

find loop start: for i in A B; do

find job marker:
#@1,0,find1,u,,sbatch -p short -c 1 -t 2:0:0
sbatch options: sbatch -p short -c 1 -t 2:0:0

find job:
grep -H John $u >>  John.txt; grep -H Mike $u >>  Mike.txt

find job marker:
#@2,0,find2,u,,sbatch -p short -c 1 -t 2:0:0
sbatch options: sbatch -p short -c 1 -t 2:0:0

find job:
grep -H Nick $u >>  Nick.txt; grep -H Julia $u >>  Julia.txt
find loop end: done

find job marker:
#@3,1.2,merge

find job:
cat John.txt Mike.txt Nick.txt Julia.txt > all.txt
smartSlurmLog/slurmPipeLine.a855454a70b2198fa5b2643bb1d41762.sh bashScriptV2.sh is ready to run. Starting to run ...
Running $smartSlurmLogDir/slurmPipeLine.a855454a70b2198fa5b2643bb1d41762.sh bashScriptV2.sh

Currently Loaded Modules:
  

---------------------------------------------------------

step: 1, depends on: 0, job name: find1, flag: find1.A reference: .u
depend on no job
sbatch -p short -c 1 -t 2:0:0 --requeue --nodes=1  -J 1.0.find1.A -o /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.A.out -e /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.A.out /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.A.sh
\# This is testing, so no job is submitted. In real run it should submit job such as: Submitted batch job 1349

step: 2, depends on: 0, job name: find2, flag: find2.A reference: .u
depend on no job
sbatch -p short -c 1 -t 2:0:0 --requeue --nodes=1  -J 2.0.find2.A -o /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/2.0.find2.A.out -e /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/2.0.find2.A.out /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/2.0.find2.A.sh
\# This is testing, so no job is submitted. In real run it should submit job such as: Submitted batch job 1560

step: 1, depends on: 0, job name: find1, flag: find1.B reference: .u
depend on no job
sbatch -p short -c 1 -t 2:0:0 --requeue --nodes=1  -J 1.0.find1.B -o /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.B.out -e /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.B.out /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.B.sh
\# This is testing, so no job is submitted. In real run it should submit job such as: Submitted batch job 1766

step: 2, depends on: 0, job name: find2, flag: find2.B reference: .u
depend on no job
sbatch -p short -c 1 -t 2:0:0 --requeue --nodes=1  -J 2.0.find2.B -o /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/2.0.find2.B.out -e /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/2.0.find2.B.out /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/2.0.find2.B.sh
\# This is testing, so no job is submitted. In real run it should submit job such as: Submitted batch job 1970

step: 3, depends on: 1.2, job name: merge , flag: merge reference:
depend on multiple jobs
sbatch -p short -t 10:0 -c 1 --requeue --nodes=1 --dependency=afterok:1349:1766:1560:1970 -J 3.1.2.merge -o /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/3.1.2.merge.out -e /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/3.1.2.merge.out /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/3.1.2.merge.sh
\# This is testing, so no job is submitted. In real run it should submit job such as: Submitted batch job 2172

All submitted jobs:
job_id       depend_on              job_flag
1349        null                  1.0.find1.A
1560        null                  2.0.find2.A
1766        null                  1.0.find1.B
1970        null                  2.0.find2.B
2172        ..1349.1766..1560.1970  3.1.2.merge
---------------------------------------------------------
Note: This is just a test run, so no job is actually submitted. In real run it should submit jobs and report as above.
Run the modified bash script as a pipeline

Thus far in the example, we have not actually submitted any jobs to the scheduler. To submit the pipeline, you will need to append the run parameter to the command. If run is not specified, test mode will be used, which does not submit jobs and gives the placeholder of 1234for jobids in the command's output. 
runAsPipeline bashScriptV2.sh "sbatch -p short -t 10:0 -c 1" useTmp run

\# Below is the output
Fri Sep 24 09:48:12 EDT 2021
Running: /home/ld32/smartSlurm/bin/runAsPipeline bashScriptV2.sh sbatch -p short -t 10:0 -c 1 useTmp run

Currently Loaded Modules:
  

converting bashScriptV2.sh to $smartSlurmLogDir/slurmPipeLine.a855454a70b2198fa5b2643bb1d41762.run.sh

find loop start: for i in A B; do

find job marker:
#@1,0,find1,u,,sbatch -p short -c 1 -t 2:0:0
sbatch options: sbatch -p short -c 1 -t 2:0:0

find job:
grep -H John $u >>  John.txt; grep -H Mike $u >>  Mike.txt

find job marker:
#@2,0,find2,u,,sbatch -p short -c 1 -t 2:0:0
sbatch options: sbatch -p short -c 1 -t 2:0:0

find job:
grep -H Nick $u >>  Nick.txt; grep -H Julia $u >>  Julia.txt
find loop end: done

find job marker:
#@3,1.2,merge

find job:
cat John.txt Mike.txt Nick.txt Julia.txt > all.txt
smartSlurmLog/slurmPipeLine.a855454a70b2198fa5b2643bb1d41762.run.sh bashScriptV2.sh is ready to run. Starting to run ...
Running $smartSlurmLogDir/slurmPipeLine.a855454a70b2198fa5b2643bb1d41762.run.sh bashScriptV2.sh

Currently Loaded Modules:
  

Could not find any jobs to cancel.
---------------------------------------------------------

step: 1, depends on: 0, job name: find1, flag: find1.A reference: .u
depend on no job
sbatch -p short -c 1 -t 2:0:0 --requeue --nodes=1  -J 1.0.find1.A -o /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.A.out -e /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.A.out /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.A.sh
\# Submitted batch job 41208893

step: 2, depends on: 0, job name: find2, flag: find2.A reference: .u
depend on no job
sbatch -p short -c 1 -t 2:0:0 --requeue --nodes=1  -J 2.0.find2.A -o /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/2.0.find2.A.out -e /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/2.0.find2.A.out /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/2.0.find2.A.sh
\# Submitted batch job 41208894

step: 1, depends on: 0, job name: find1, flag: find1.B reference: .u
depend on no job
sbatch -p short -c 1 -t 2:0:0 --requeue --nodes=1  -J 1.0.find1.B -o /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.B.out -e /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.B.out /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.B.sh
\# Submitted batch job 41208895

step: 2, depends on: 0, job name: find2, flag: find2.B reference: .u
depend on no job
sbatch -p short -c 1 -t 2:0:0 --requeue --nodes=1  -J 2.0.find2.B -o /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/2.0.find2.B.out -e /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/2.0.find2.B.out /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/2.0.find2.B.sh
\# Submitted batch job 41208898

step: 3, depends on: 1.2, job name: merge , flag: merge reference:
depend on multiple jobs
sbatch -p short -t 10:0 -c 1 --requeue --nodes=1 --dependency=afterok:41208893:41208895:41208894:41208898 -J 3.1.2.merge -o /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/3.1.2.merge.out -e /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/3.1.2.merge.out /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/3.1.2.merge.sh
\# Submitted batch job 41208899

All submitted jobs:
job_id       depend_on              job_flag
41208893    null                  1.0.find1.A
41208894    null                  2.0.find2.A
41208895    null                  1.0.find1.B
41208898    null                  2.0.find2.B
41208899    ..41208893.41208895..41208894.41208898  3.1.2.merge
---------------------------------------------------------
Monitoring the jobs

You can use the command:
O2squeue -u $USER

To see the job status (running, pending, etc.). You also get two emails for each step, one at the start of the step, one at the end of the step.
Successful job email
Email subject: Success: job id:41208893 name:1.0.find1.A

Email content:

Job script content:
#!/bin/bash
#Commands:
trap "{ cleanup.sh /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.A; }” EXIT
touch /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.A.start
srun -n 1 bash -e -c "{ set -e; rsyncToTmp  /tmp/rcbio/universityA.txt; grep -H John /tmp/rcbio/universityA.txt >>  John.txt; grep -H Mike /tmp/rcbio/universityA.txt >>  Mike.txt        ; } && touch /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.A.success || touch /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.A.failed"

#sbatch command:
#sbatch -p short -c 1 -t 2:0:0 --requeue --nodes=1  -J 1.0.find1.A -o /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.A.out -e /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.A.out /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.A.sh

\# Submitted batch job 41208893
Job output:
Working to copy: /tmp/rcbio/universityA.txt, waiting lock...
Got lock: /tmp/-tmp-rcbio-universityA.txt. Copying data to: /tmp/rcbio/universityA.txt
Copying is done for /tmp/rcbio/universityA.txt
Job done. Summary:
       JobID              Submit               Start                 End      State  Partition              ReqTRES  Timelimit    CPUTime     MaxRSS                       NodeList
------------ ------------------- ------------------- ------------------- ---------- ---------- -------------------- ---------- ---------- ---------- ------------------------------
41208893     2021-09-24T09:48:13 2021-09-24T09:48:24             Unknown    RUNNING      short billing=1,cpu=1,mem+   00:2:0:0   00:00:10                          compute-x
41208893.ba+ 2021-09-24T09:48:24 2021-09-24T09:48:24             Unknown    RUNNING                                              00:00:10                          compute-x
41208893.ex+ 2021-09-24T09:48:24 2021-09-24T09:48:24             Unknown    RUNNING                                              00:00:10                          compute-x
41208893.0   2021-09-24T09:48:29 2021-09-24T09:48:29 2021-09-24T09:48:29  COMPLETED                                              00:00:00          0               compute-x
*Notice the sacct report above: while the main job is still running for sacct command, user task is completed.

The key elements are time and memory used.
Check job logs

You can use the command:
ls -l flag

This command list all the logs created by the pipeline runner. *.sh files are the slurm scripts for each step, *.out files are output files for each step, *.success files means job successfully finished for each step and *.failed means job failed for each steps.

You also get two emails for each step, one at the start of the step, one at the end of the step.
Cancel all jobs

You can use the command to cancel running and pending jobs:
cancelAllJobs $smartSlurmLogDir/alljobs.jid
What happens if there is some error? 

You can re-run this command in the same folder. We will delete an input file to see what happens.
\# We are intentionally removing an input file to see a "failed job" email message
rm universityB.txt
runAsPipeline bashScriptV2.sh "sbatch -p short -t 10:0 -c 1" useTmp run

\# Here is the output
Fri Sep 24 10:00:36 EDT 2021
Running: /home/ld32/smartSlurm/bin/runAsPipeline bashScriptV2.sh sbatch -p short -t 10:0 -c 1 useTmp run

Currently Loaded Modules:
  

This is a re-run with the same command and script is not changed, no need to convert the script. Using the old one: $smartSlurmLogDir/slurmPipeLine.a855454a70b2198fa5b2643bb1d41762.run.sh
Running $smartSlurmLogDir/slurmPipeLine.a855454a70b2198fa5b2643bb1d41762.run.sh bashScriptV2.sh

Currently Loaded Modules:
  

Could not find any jobs to cancel.
---------------------------------------------------------

step: 1, depends on: 0, job name: find1, flag: find1.A reference: .u
depend on no job
1.0.find1.A was done before, do you want to re-run it?
y:        To re-run this job, press y, then enter key.
ystep:    To re-run all jobs for step 3: hisatCount, type yall, then press enter key.
yall:     To re-run all jobs, type yallall, then press enter key.
enter:    To not re-run this job, directly press enter key.
nstep:    To not re-run all successful jobs for step 3: hisatCount, type nall, then press enter key.
nall:     To not re-run all successful jobs, type nallall, then press enter key.

\# type enter here to not re-run

step: 2, depends on: 0, job name: find2, flag: find2.A reference: .u
depend on no job
2.0.find2.A was done before, do you want to re-run it?
y:        To re-run this job, press y, then enter key.
ystep:    To re-run all jobs for step 3: hisatCount, type yall, then press enter key.
yall:     To re-run all jobs, type yallall, then press enter key.
enter:    To not re-run this job, directly press enter key.
nstep:    To not re-run all successful jobs for step 3: hisatCount, type nall, then press enter key.
nall:     To not re-run all successful jobs, type nallall, then press enter key.

\# type enter here to not re-run

job 2.0.find2.A is not submitted

step: 1, depends on: 0, job name: find1, flag: find1.B reference: .u
depend on no job
1.0.find1.B was done before, do you want to re-run it?
y:        To re-run this job, press y, then enter key.
ystep:    To re-run all jobs for step 3: hisatCount, type yall, then press enter key.
yall:     To re-run all jobs, type yallall, then press enter key.
enter:    To not re-run this job, directly press enter key.
nstep:    To not re-run all successful jobs for step 3: hisatCount, type nall, then press enter key.
nall:     To not re-run all successful jobs, type nallall, then press enter key.

\# type ‘y’ and enter here to re-run

Will re-run the down stream steps even if they are done before.
sbatch -p short -c 1 -t 2:0:0 --requeue --nodes=1  -J 1.0.find1.B -o /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.B.out -e /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.B.out /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.B.sh
\# Submitted batch job 41209197

step: 2, depends on: 0, job name: find2, flag: find2.B reference: .u
depend on no job
2.0.find2.B was done before, do you want to re-run it?
y:        To re-run this job, press y, then enter key.
ystep:    To re-run all jobs for step 3: hisatCount, type yall, then press enter key.
yall:     To re-run all jobs, type yallall, then press enter key.
enter:    To not re-run this job, directly press enter key.
nstep:    To not re-run all successful jobs for step 3: hisatCount, type nall, then press enter key.
nall:     To not re-run all successful jobs, type nallall, then press enter key.

\# type enter here to not re-run

job 2.0.find2.B is not submitted

step: 3, depends on: 1.2, job name: merge , flag: merge reference:
depend on other jobs
sbatch -p short -t 10:0 -c 1 --requeue --nodes=1 --dependency=afterok:41209197 -J 3.1.2.merge -o /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/3.1.2.merge.out -e /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/3.1.2.merge.out /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/3.1.2.merge.sh
\# Submitted batch job 41209210

\# Notice above, rcbio didn’t ask if user wants to re-run step3 or not and directly re-run it.

All submitted jobs:
job_id       depend_on              job_flag
41209197    null                  1.0.find1.B
41209210    ..41209197.           3.1.2.merge
---------------------------------------------------------

This command will check if the earlier run is finished or not. If not, ask user to kill the running jobs or not, then ask user to rerun the successfully finished steps or not. Click 'y', it will rerun, directly press 'enter' key, it will not rerun. 
Failed job email
Email subject: Failed: job id:41209197 name:1.0.find1.B

Email content:
Job script content:
#!/bin/bash
#Commands:
trap "{ cleanup.sh /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.B; }” EXIT
touch /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.B.start
srun -n 1 bash -e -c "{ set -e; rsyncToTmp  /tmp/rcbio/universityB.txt; grep -H John /tmp/rcbio/universityB.txt >>  John.txt; grep -H Mike /tmp/rcbio/universityB.txt >>  Mike.txt        ; } && touch /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.B.success || touch /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.B.failed"

#sbatch command:
#sbatch -p short -c 1 -t 2:0:0 --requeue --nodes=1  -J 1.0.find1.B -o /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.B.out -e /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.B.out /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.B.sh

# Submitted batch job 41209197
Job output:
Working to copy: /tmp/rcbio/universityB.txt, waiting lock...
Reference file or folder not exist: /universityB.txt
grep: /tmp/rcbio/universityB.txt: No such file or directory
grep: /tmp/rcbio/universityB.txt: No such file or directory
Job done. Summary:
       JobID              Submit               Start                 End      State  Partition              ReqTRES  Timelimit    CPUTime     MaxRSS                       NodeList
------------ ------------------- ------------------- ------------------- ---------- ---------- -------------------- ---------- ---------- ---------- ------------------------------
41209197     2021-09-24T10:02:43 2021-09-24T10:03:09             Unknown    RUNNING      short billing=1,cpu=1,mem+   00:2:0:0   00:00:09                          compute-x
41209197.ba+ 2021-09-24T10:03:09 2021-09-24T10:03:09             Unknown    RUNNING                                              00:00:09                          compute-x
41209197.ex+ 2021-09-24T10:03:09 2021-09-24T10:03:09             Unknown    RUNNING                                              00:00:09                          compute-x
41209197.0   2021-09-24T10:03:13 2021-09-24T10:03:13 2021-09-24T10:03:13  COMPLETED                                              00:00:00          0               compute-x
*Notice the sacct report above: while the main job is still running for sacct command, user task is completed.

    The key element here is the error message.

    Notice here, step2 job is automatically canceled because this job failed. We deleted universityB.txt, so the job has failed. We don’t get an email from the downstream step3 job. 

Fix the error and re-run the pipeline

You can rerun this command in the same folder
cp universityA.txt universityB.txt
runAsPipeline bashScriptV2.sh "sbatch -p short -t 10:0 -c 1" useTmp run

This command will automatically check if the earlier run is finished. If the run has not finished, the script will ask the user if they want to kill the running jobs or not, then ask user to rerun the successfully finished steps or not. Click 'y', it will rerun, directly press 'enter' key, it will not rerun. 

Notice here, step3 will run by default. It will run without prompting the user for permission.
What happens if we add more input data and re-run the pipeline?

You can rerun this command in the same folder
cp universityA.txt universityC.txt
cp bashScriptV2.sh bashScriptV3.sh 
nano bashScriptV3.sh  
\# change
for i in A B; do
to: 
for i in A B C; do

\# save the file and run:
runAsPipeline bashScriptV3.sh "sbatch -p short -t 10:0 -c 1" useTmp run

\# Here are the output: 
Fri Sep 24 10:56:16 EDT 2021
Running: /home/ld32/smartSlurm/bin/runAsPipeline bashScriptV3.sh sbatch -p short -t 10:0 -c 1 useTmp run

Currently Loaded Modules:
  

converting bashScriptV3.sh to $smartSlurmLogDir/slurmPipeLine.b72e7f91da30d312a2c85d0735896f79.run.sh

find loop start: for i in A B C; do

find job marker:
\#@1,0,find1,u,,sbatch -p short -c 1 -t 2:0:0
sbatch options: sbatch -p short -c 1 -t 2:0:0

find job:
grep -H John $u >>  John.txt; grep -H Mike $u >>  Mike.txt

find job marker:
\#@2,0,find2,,u,,sbatch -p short -c 1 -t 2:0:0
sbatch options: sbatch -p short -c 1 -t 2:0:0

find job:
grep -H Nick $u >>  Nick.txt; grep -H Julia $u >>  Julia.txt
find loop end: done

find job marker:
\#@3,1.2,merge

find job:
cat John.txt Mike.txt Nick.txt Julia.txt > all.txt
smartSlurmLog/slurmPipeLine.b72e7f91da30d312a2c85d0735896f79.run.sh bashScriptV3.sh is ready to run. Starting to run ...
Running $smartSlurmLogDir/slurmPipeLine.b72e7f91da30d312a2c85d0735896f79.run.sh bashScriptV3.sh

Currently Loaded Modules:
  


Could not find any jobs to cancel.
---------------------------------------------------------

step: 1, depends on: 0, job name: find1, flag: find1.A reference: .u
depend on no job
1.0.find1.A was done before, do you want to re-run it?
y:        To re-run this job, press y, then enter key.
ystep:    To re-run all jobs for step 3: hisatCount, type yall, then press enter key.
yall:     To re-run all jobs, type yallall, then press enter key.
enter:    To not re-run this job, directly press enter key.
nstep:    To not re-run all successful jobs for step 3: hisatCount, type nall, then press enter key.
nall:     To not re-run all successful jobs, type nallall, then press enter key.

\# type enter here to not re-run

job 1.0.find1.A is not submitted

step: 2, depends on: 0, job name: find2, flag: find2.A reference: .u
depend on no job
2.0.find2.A was done before, do you want to re-run it?
y:        To re-run this job, press y, then enter key.
ystep:    To re-run all jobs for step 3: hisatCount, type yall, then press enter key.
yall:     To re-run all jobs, type yallall, then press enter key.
enter:    To not re-run this job, directly press enter key.
nstep:    To not re-run all successful jobs for step 3: hisatCount, type nall, then press enter key.
nall:     To not re-run all successful jobs, type nallall, then press enter key.

\# type enter here to not re-run

job 2.0.find2.A is not submitted

step: 1, depends on: 0, job name: find1, flag: find1.B reference: .u
depend on no job
1.0.find1.B was done before, do you want to re-run it?
y:        To re-run this job, press y, then enter key.
ystep:    To re-run all jobs for step 3: hisatCount, type yall, then press enter key.
yall:     To re-run all jobs, type yallall, then press enter key.
enter:    To not re-run this job, directly press enter key.
nstep:    To not re-run all successful jobs for step 3: hisatCount, type nall, then press enter key.
nall:     To not re-run all successful jobs, type nallall, then press enter key.

\# type enter here to not re-run

job 1.0.find1.B is not submitted

step: 2, depends on: 0, job name: find2, flag: find2.B reference: .u
depend on no job
2.0.find2.B was done before, do you want to re-run it?
y:        To re-run this job, press y, then enter key.
ystep:    To re-run all jobs for step 3: hisatCount, type yall, then press enter key.
yall:     To re-run all jobs, type yallall, then press enter key.
enter:    To not re-run this job, directly press enter key.
nstep:    To not re-run all successful jobs for step 3: hisatCount, type nall, then press enter key.
nall:     To not re-run all successful jobs, type nallall, then press enter key.

\# type enter here to not re-run

job 2.0.find2.B is not submitted

step: 1, depends on: 0, job name: find1, flag: find1.C reference: .u
depend on no job
sbatch -p short -c 1 -t 2:0:0 --requeue --nodes=1  -J 1.0.find1.C -o /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.C.out -e /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.C.out /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/1.0.find1.C.sh
\# Submitted batch job 41211380

step: 2, depends on: 0, job name: find2, flag: find2.C reference: .u
depend on no job
sbatch -p short -c 1 -t 2:0:0 --requeue --nodes=1  -J 2.0.find2.C -o /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/2.0.find2.C.out -e /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/2.0.find2.C.out /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/2.0.find2.C.sh
\# Submitted batch job 41211381

step: 3, depends on: 1.2, job name: merge , flag: merge reference:
depend on multiple jobs
sbatch -p short -t 10:0 -c 1 --requeue --nodes=1 --dependency=afterok:41211380:41211381 -J 3.1.2.merge -o /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/3.1.2.merge.out -e /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/3.1.2.merge.out /home/ld32/testRunBashScriptAsSlurmPipeline/smartSlurmLog/3.1.2.merge.sh
\# Submitted batch job 41211382

All submitted jobs:
job_id       depend_on              job_flag
41211380    null                  1.0.find1.C
41211381    null                  2.0.find2.C
41211382    ..41211380..41211381  3.1.2.merge
---------------------------------------------------------

This command will check if the earlier run is finished, and will prompt the user if they kill any running jobs. Next, it will then ask the user if they want to rerun any successfully finished steps. Click 'y', it will rerun, directly press 'enter' key, it will not rerun. 

For the new data, RCBio will submit 2 jobs. Step3 will also still automatically run.
Re-run a single job manually
\# /working/directory is a placeholder, replace it with your actual working directory path
cd /working/directory
\# all/related/modules is a placeholder, replace it with the actual other modules/versions you need
module load rcbio/1.3.3 and all/related/modules

\# submit job with proper partition, time, number of cores and memory
sbatch --requeue --mail-type=ALL -p short -t 2:0:0 -c 2 --mem 2G /working/directory/smartSlurmLog/stepID.loopID.stepName.sh

Or:
runSingleJob "module load bowtie/1.2.2; bowtie -x /n/groups/shared_databases/bowtie_indexes/hg19 -p 2 -1 read1.fq -2 read2.fq --sam > out.bam" "sbatch -p short -t 1:0:0 -c 2 -mem 8G"

For details about the second option: Get more informative slurm email notification and logs through rcbio/1.3 
To run your own script as Slurm pipeline

If you have a bash script with multiple steps and you wish to run it as Slurm pipeline, here is how you can do that:

    modify your old script and add the notation to mark the start and end of any loops, and the start of any step for which you want to submit as an sbatch job. 

    use runAsPipeline with your modified bash script, as detailed above. 

How does the runAsPipeline RCBio pipeline runner work?

In case you wonder how it works, here is a simple example to explain.


## How does smart pipeline work
[Back to top](#SmartSlurm)

1) Auto adjust memory and run-time according to statistics from earlier jobs

$smartSlurmJobRecordDir/jobRecord.txt contains job memory and run-time records. There are three important columns: 
   
   2rd colume is input size,
   
   7th column is final memory usage
   
   8th column is final time usage
   
The data from the three columns are plotted and statistics  
__________________________________________________________________________________________________________________   
1jobID,2inputSize,3mem,4time,5mem,6time,7mem,8time,9status,10useID,11path,12software,13reference,14output,15script,16error,17cpu,18node,19date,20command
46531,1465,4G,2:0:0,4G,0-2:0:0,3.52,1,COMPLETED,ld32,,useMemTimeWithInput,none,slurm-%j.out slurm-YRTrRAYA.sh slurm-%j.err,1,compute-x,slurm-46531.err,Tue Dec 6 15:29:20 EST 2022,"ssbatch -p short -t 2:0:0 --mem=4G --wrap useMemTimeWithInput.sh bigText1.txt run"

46535,2930,4G,2:0:0,4G,0-2:0:0,6.38,2,COMPLETED,ld32,,useMemTimeWithInput,none,slurm-%j.out slurm-oT42tyEE.sh slurm-%j.err,1,compute-x,slurm-46535.err,Tue Dec 6 15:30:46 EST 2022,"ssbatch -p short -t 2:0:0 --mem=4G --wrap useMemTimeWithInput.sh bigText2.txt run"

46534,4395,4G,2:0:0,4G,0-2:0:0,9.24,4,COMPLETED,ld32,,useMemTimeWithInput,none,slurm-%j.out slurm-TQyBOQ5f.sh slurm-%j.err,1,compute-x,slurm-46534.err,Tue Dec 6 15:32:40 EST 2022,"ssbatch -p short -t 2:0:0 --mem=4G --wrap useMemTimeWithInput.sh bigText3.txt run"

\#Here is the input size vs memory plot for useMemTimeWithInput: 

![](https://github.com/ld32/smartSlurm/blob/main/stats/useMemTimeWithInput.none.mem.png)

\#Here is the input size vs run-time plot for useMemTimeWithInput: 

![](https://github.com/ld32/smartSlurm/blob/main/stats/useMemTimeWithInput.none.time.png)

\#Here is the run-time vs memory plot for useMemTimeNoInput: 

![](https://github.com/ld32/smartSlurm/blob/main/stats/useMemTimeNoInput.none.stat.noInput.png)

2) Auto choose partition according to run-time request

smartSlrm/config/config.txt contains partion time limit and bash function adjustPartition to adjust partion for sbatch jobs: 

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

...   

\#function 

adjustPartition() {         
    ... # please open the file to see the content         
} ; export -f adjustPartition 

3) Auto re-run failed OOM (out of memory) and OOT (out of run-time) jobs
    
    At end of the job, $smartSlurmJobRecordDir/bin/cleanUp.sh checkes memory and time usage, save the data in to log $smartSlurmJobRecordDir/myJobRecord.txt. If job fails, ssbatch re-submit with double memory or double time, clear up the statistic fomular, so that later jobs will re-caculate statistics, 

6) Get good emails: by default Slurm emails only have a subject. ssbatch attaches the content of the sbatch script, the output and error log to email

    $smartSlurmJobRecordDir/bin/cleanUp.sh also sends a email to user. The email contains the content of the Slurm script, the sbatch command used, and also the content of the output and error log files.



=======

# sbatchAndTop
## How to use sbatchAndTop
[Back to top](#SmartSlurm)

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
