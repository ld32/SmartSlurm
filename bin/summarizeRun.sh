#!/bin/sh

Usage="Usage: $0 [jobRecordFile, optional, for example: log/allJobs.txt.old]\nNote: this script will go through job name list in a file (eg: log/alljobs), to see if the jobs finish successfully or not." 

if [ ! -z "$1" ]; then 
    [ -f "$1" ] || { echo $Usage; exit 1; }
     
    echo According to $1:
    names=`tail -n +1 $1 | awk '{print $3}' | tr "\n" " "`
    for name in $names; do       
       if [ -f log/$name.start ]; then 
           [ -f log/$name.success ] && echo Success! $name  || { echo "Failed! $name #############" && [ -f log/$name.err ] && echo -e "Here is the content of log/$name.err :\n" &&  tail -n 10 log/$name.err; }
       else 
          echo Didn not run! $name 
       fi   
    done
    exit 
fi       

[ -f log/allJobs.txt.first ] || { echo Job list file not exist: log/allJobs.txt.first; exit 1; }

echo According to log/allJobs.txt.first:
names=`tail -n +1 log/allJobs.txt.first | awk '{print $3}' | tr "\n" " "`
for name in $names; do       
        if [ -f log/$name.start ]; then 
           [ -f log/$name.success ] && echo Success! $name  || { echo "Failed! $name #############" && [ -f log/$name.err ] && echo -e "Here is the content of log/$name.err :\n" &&  tail -n 10 log/$name.err; }
       else 
          echo Didn not run! $name 
       fi   
done
echo 
echo According to log/allJobs.txt: 
names=`tail -n +1 log/allJobs.txt | awk '{print $3}' | tr "\n" " "`
for name in $names; do       
       if [ -f log/$name.start ]; then 
           [ -f log/$name.success ] && echo Success! $name  || { echo "Failed! $name #############" && [ -f log/$name.err ] && echo -e "Here is the content of log/$name.err :\n" &&  tail -n 10 log/$name.err; }
       else 
          echo Didn not run! $name 
       fi   
 
done 

echo 
echo According to log/allJobs.txt: 
names=`tail -n +1 log/allJobs.txt | awk '{print $3}' | tr "\n" " "`
for name in $names; do       
        if [ -f log/$name.start ]; then 
           [ -f log/$name.success ] && echo Success! $name  || { echo "Failed! $name #############" && [ -f log/$name.err ] && echo -e "Here is the content of log/$name.err :\n" &&  tail -n 10 log/$name.err; }
       else 
          echo Didn not run! $name 
       fi   
done 

# check running jobs 
checkJobsSlurm  log/allJobs.txt
