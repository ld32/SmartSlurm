#!/bin/bash

start=$SECONDS

size=$(($1 * 5))
delay=$(($1 * 30))

echo "Begin allocating memory..."
for index in $(seq $size); do
    value=$(seq -w -s '' $index $(($index + 100000)))
    eval array$index=$value
    echo $value >> bigText.txt
done

delay=$(( start + delay - SECONDS ))

echo "...end allocating memory. Begin sleeping... for $delay seonds"

sleep $delay 

echo "Done"

