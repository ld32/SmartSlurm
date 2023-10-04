#!/bin/sh

#loopStart:i
for i in {1..3}; do
    #@1,0,step1,,input,sbatch -p short -c 1 --mem 2G -t 50:0
    echo loop i=$i step 1
   
    #loopStart:j
    for j in {1..2}; do 
        #@2,1,step2,,input,sbatch -p short -c 1 --mem 2G -t 50:0
        echo loop j=$j step2  
    
    done

    #@3,2,step3,,input,sbatch -p short -c 1 --mem 2G -t 50:0
    echo loop i=$i step3
done

#@4,3,step4,,input
echo step 4
