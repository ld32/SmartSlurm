#!/bin/bash

echo "Begin..."
rm bigText1.txt 2>/dev/null
for index in $(seq 50000); do
    value="$value John Mike Smith David\n"
    
done
echo -e "$value" >> bigText1.txt

for count in {2..5}; do
    rm bigText$count.txt 2>/dev/null
    for c in $(seq $count); do   
        cat bigText1.txt >> bigText$count.txt 
    done    
done
echo "Done"

ls -lh bigText*

