#!/bin/sh

for i in A B; do            
    
    u=university$i.txt   
    
    grep -H John $u >>  John.txt; grep -H Mike $u >>  Mike.txt       
    
    grep -H Nick $u >>  Nick.txt; grep -H Julia $u >>  Julia.txt
done
    
cat John.txt Mike.txt Nick.txt Julia.txt > all.txt

