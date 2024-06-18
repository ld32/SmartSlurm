#!/bin/bash

echo "Begin..."
rm numbers1.txt 2>/dev/null
for index in $(seq 5); do
    value=$(seq -w -s '\n' $index $(($index + 50000)))
    echo -e "$value" >> numbers1.txt
done

for count in {2..5}; do
    rm numbers$count.txt 2>/dev/null
    for c in $(seq $count); do   
        cat numbers1.txt >> numbers$count.txt 
    done    
done

for count in {6..10}; do
    cp numbers2.txt  numbers$count.txt
done

echo "Done"

ls -lh numbers*
