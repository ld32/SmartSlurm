#!/bin/bash

#set -x 

echo "Start..."

for count in {1..5}; do
    value=""
    if [[ $count == 1 ]]; then 
        for index in {1..5000}; do
            value="John Paul\nMike Smith\nNick Will\nJulia Johnson\nJudy Jones\n$value "
        done
    else 
        for c in {1..$count}; do   
            value="`cat university1.txt`$value "  
        done
    fi         
    echo -e "$value" > university$count.txt
    echo Created university$count.txt
done
echo "Done"


