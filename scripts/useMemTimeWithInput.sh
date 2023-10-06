#!/bin/bash

start=$SECONDS

totalSize=0

echo "Begin allocating memory..."

for input in "$@"; do

    txt=$txt.`cat $input`
    
    fileSize=`du --apparent-size -c -L $input | tail -n 1 | cut -f 1`
    
    totalSize=$(( fileSize + totalSize ))

done 

totalSize=$(( $totalSize/1024 )); # convert to M turned

delay=$((start + totalSize * 60 - SECONDS))  # 60 seconds per

echo "...end allocating memory. Begin sleeping for $delay seconds..."

[ "$delay" -ge 1 ] && sleep $delay

echo "Done"

