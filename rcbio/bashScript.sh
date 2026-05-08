#!/bin/sh

for i in A B; do            
    
    echo John >> John.txt; echo Mike >>  Mike.txt       
    
    echo Nick >> Nick.txt; echo Julia >>  Julia.txt
    
done
    
cat John.txt Mike.txt Nick.txt Julia.txt > all.txt

