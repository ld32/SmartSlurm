#!/bin/ash

size=$(($1 * 1000))
delay=$2

echo "begin allocating memory..."
for index in $(seq $size); do
    value=$(seq -w -s '' $index $(($index + 100000)))
    eval array$index=$value
done
echo "...end allocating memory"

echo "sleeping for $delay"
sleep $delay
