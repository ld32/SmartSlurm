#!/bin/bash

start=$SECONDS

size=$(($1 * 1000))
delay=$(($2 * 60))

echo "begin allocating memory..."
for index in $(seq $size); do
    value=$(seq -w -s '' $index $(($index + 100000)))
    eval array$index=$value
done
echo "...end allocating memory"

while [ $(( SECONDS - start )) -le $delay ]; do 
    echo Sleeping for 5 seconds...
    sleep 5; 
done
