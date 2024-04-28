#!/bin/bash

#set -x 

totalSize=0

echo "Start allocating memory..." >&2

for input in "$2 $3 $4"; do
    [ -z "$input" ] && break

    txt=$(cat $input); txt="$txt.$txt.$txt$txt.$txt.$txt"
    
    fileSize=`du --apparent-size -c -L $input | tail -n 1 | cut -f 1`
    
    totalSize=$(( fileSize + totalSize ))

    sleep 10
done          

grep $1 university1.txt

delay=$((totalSize * 3 / 100))  

echo "...end allocating memory. Begin sleeping for $delay seconds..." >&2

sleep $delay

echo "Done" >&2