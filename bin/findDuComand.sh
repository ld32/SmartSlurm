#!/bin/bash

#set -x

nodes=$(sinfo -hNl | awk '{print $1}' | uniq | tr '\n' ' ')

log=$HOME/ld32.$(date '+%Y%m%d%H%M%S').txt
date >> $log
for node in $nodes; do

    echo "Computing node: $node" >> $log

    ssh $node "ps -Af | grep -v '$USER\|root' | grep ' du -'" >> $log 2>&1 || echo ssh failed or no du is found. >> $log

    # break

    #echo checking $node done >> $log
done

date >> $log

#cat $log
