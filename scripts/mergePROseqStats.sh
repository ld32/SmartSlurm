#!/bin/sh

#set -x

echo start
folders=""
for i in `tail -n+2 $1 | cut -d',' -f1,10`; do 
    echo working on $i; 
    if [[ "$folders" == *" ${i#*,} "* ]]; then 
        echo folder already exist
        cat ${i#*,}/logs/${i%,*}.stats.txt >>  ${i#*,}/logs/stats.txt
    else
        folders="$folders ${i#*,} "
        if [[ "$2" == "true" ]]; then 
            echo -e Sample"\t"Trim \#"\t"Trim %"\t"Spike \#"\t"Spike %"\t"Ref \#"\t"Ref %"\t"TSS \#"\t"TSS %"\t"Insert > ${i#*,}/logs/stats.txt
        else 
            echo -e Sample"\t"Trim \#"\t"Trim %"\t"Spike \#"\t"Spike %"\t"Ref \#"\t"Ref %"\t"TSS \#"\t"TSS %"\t"Insert"\t"Spike dedup %"\t"Ref dedup % > ${i#*,}/logs/stats.txt
        fi
        cat ${i#*,}/logs/${i%,*}.stats.txt >>  ${i#*,}/logs/stats.txt

    fi
done
