#!/bin/sh

Usage="Usage: $0 [jobRecordFile, optional, for example: flag/allJobs.txt.old]\nNote: this script will go through job name list in a file (eg: flag/alljobs), to see if the jobs finish successfully or not." 

if [ ! -z "$1" ]; then 
    [ -f "$1" ] || { echo $Usage; exit 1; }
     
    echo According to $1:
    names=`tail -n +1 $1 | awk '{print $3}' | tr "\n" " "`
    for name in $names; do       
       if [ -f flag/$name.start ]; then 
           [ -f flag/$name.success ] && echo Success! $name  || { echo "Failed! $name #############" && [ -f flag/$name.err ] && echo -e "Here is the content of flag/$name.err :\n" &&  tail -n 10 flag/$name.err; }
       else 
          echo Didn not run! $name 
       fi   
    done
    exit 
fi       

[ -f flag/allJobs.txt.first ] || { echo Job list file not exist: flag/allJobs.txt.first; exit 1; }

echo According to flag/allJobs.txt.first:
names=`tail -n +1 flag/allJobs.txt.first | awk '{print $3}' | tr "\n" " "`
for name in $names; do       
        if [ -f flag/$name.start ]; then 
           [ -f flag/$name.success ] && echo Success! $name  || { echo "Failed! $name #############" && [ -f flag/$name.err ] && echo -e "Here is the content of flag/$name.err :\n" &&  tail -n 10 flag/$name.err; }
       else 
          echo Didn not run! $name 
       fi   
done
echo 
echo According to flag/allJobs.txt: 
names=`tail -n +1 flag/allJobs.txt | awk '{print $3}' | tr "\n" " "`
for name in $names; do       
       if [ -f flag/$name.start ]; then 
           [ -f flag/$name.success ] && echo Success! $name  || { echo "Failed! $name #############" && [ -f flag/$name.err ] && echo -e "Here is the content of flag/$name.err :\n" &&  tail -n 10 flag/$name.err; }
       else 
          echo Didn not run! $name 
       fi   
 
done 

echo 
echo According to flag/allJobs.txt: 
names=`tail -n +1 flag/allJobs.txt | awk '{print $3}' | tr "\n" " "`
for name in $names; do       
        if [ -f flag/$name.start ]; then 
           [ -f flag/$name.success ] && echo Success! $name  || { echo "Failed! $name #############" && [ -f flag/$name.err ] && echo -e "Here is the content of flag/$name.err :\n" &&  tail -n 10 flag/$name.err; }
       else 
          echo Didn not run! $name 
       fi   
done 

# check running jobs 
checkJobsSlurm  flag/allJobs.txt
