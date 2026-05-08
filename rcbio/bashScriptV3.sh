#!/bin/sh

for i in A B C; do            
    
    u=university$i.txt    
    
    #@1,0,find1,u,sbatch -p short -c 1 -t 50:0   
    grep -H John $u >>  John.txt; grep -H Mike $u >>  Mike.txt        
  
    #@2,0,find2,u,sbatch -p short -c 1 -t 50:0
    grep -H Nick $u >>  Nick.txt; grep -H Julia $u >>  Julia.txt

done 
    
#@3,1.2,merge          
cat John.txt Mike.txt Nick.txt Julia.txt > all.txt
