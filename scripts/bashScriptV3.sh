#!/bin/sh

studentName=$1

[ -z "$studentName" ] && echo -e "Error: student name is missing.\nUsage: bashScript <student name>" && exit 1

for i in {1..3}; do
    
    input=university$i.txt
    
    #@1,0,findStudent,,input,sbatch -p short -c 1 --mem 2G -t 50:0 
    findStudent.sh John $input > $studentName.$i.txt
  
done

#@2,1,mergeStudent,,,sbatch -p short -c 1 --mem 2G -t 50:0 
cat $studentName.*.txt > all$studentName.txt