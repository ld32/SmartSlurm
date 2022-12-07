#!/bin/bash

echo "Begin..."
rm bigText1.txt 2>/dev/null
for index in $(seq 5); do
    value=$(seq -w -s '\n' $index $(($index + 50000)))
    echo -e "$value" >> bigText1.txt
done

for count in {2..5}; do
    rm bigText$count.txt 2>/dev/null
    for c in $(seq $count); do   
        cat bigText1.txt >> bigText$count.txt 
    done    
done
echo "Done"

ls -lh bigText*

