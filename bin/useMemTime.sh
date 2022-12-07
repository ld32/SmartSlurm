#!/bin/bash

start=$SECONDS

size=$(($1 * 50))
delay=$(($2 * 30))

echo "Begin allocating memory..."
for index in $(seq $size); do
    value=$(seq -w -s '' $index $(($index + 100000)))
    eval array$index=$value
done
echo "...end allocating memory. Begin sleeping..."

delay=$(( start + delay - SECONDS )) 

[ "$delay" -ge 1 ] && sleep $delay

echo "Done"

