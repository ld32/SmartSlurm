#!/bin/sh
# rap_shapes fixture -- ONE pipeline exercising the three rAP DAG shapes in a single
# run: a LINEAR chain, a FAN-OUT/FAN-IN pair (nested loop), and an IF/ELSE that submits
# one of two alternative terminal steps (distinct step numbers, branched on a submit-
# time value). Self-contained: every step body uses only shell tools.
#@pipeline rapShapes
manifest=$1
[ -z "$manifest" ] && { echo "usage: $0 <manifest>" >&2; exit 1; }

samples=$(cut -d' ' -f1 "$manifest")
while read sid rest; do echo "$rest" > doc.$sid.txt; done < "$manifest"

for s in $samples; do
    in=doc.$s.txt

    # --- LINEAR chain: normalize -> tokenize ---
    #@1,0,normalize,,in,sbatch -p short -c 1 --mem 100M -t 5:00
    tr 'A-Z' 'a-z' < $in | tr -s ' ' | sed 's/^ *//; s/ *$//' > norm.$s.txt
    #@end

    #@2,1,tokenize,,in,sbatch -p short -c 1 --mem 100M -t 5:00
    tr ' ' '\n' < norm.$s.txt | sed '/^$/d' > tokens.$s.txt
    #@end

    # --- FAN-OUT prep: round-robin tokens into 3 chunk files (explicit names) ---
    #@3,2,split,,in,sbatch -p short -c 1 --mem 100M -t 5:00
    touch chunk.$s.1 chunk.$s.2 chunk.$s.3
    awk -v s="$s" '{ f = "chunk." s "." (((NR-1)%3)+1); print > f }' tokens.$s.txt
    #@end

    # --- FAN-OUT: one countchunk job per chunk (both loop vars in the flag) ---
    for c in 1 2 3; do
        #@4,3,countchunk,,in,sbatch -p short -c 1 --mem 100M -t 5:00
        w=$(wc -w < chunk.$s.$c)
        printf 'chunk %s words=%s\n' "$c" "$w" > cc.$s.$c.txt
        #@end
    done

    # --- FAN-IN: merge waits for ALL 3 chunk jobs ---
    #@5,4,merge,,in,sbatch -p short -c 1 --mem 100M -t 5:00
    tw=$(cat cc.$s.*.txt | grep -o 'words=[0-9]*' | cut -d= -f2 | awk '{ s+=\$1 } END{ print s+0 }')
    printf 'sample=%s total_words=%s\n' "$s" "$tw" > totals.$s.txt
    #@end

    # --- IF/ELSE: branch on a SUBMIT-time value (input word count); distinct steps ---
    words=$(wc -w < $in)
    if [ "$words" -ge 5 ]; then
        #@6,5,detailed,,in,sbatch -p short -c 1 --mem 100M -t 5:00
        printf '%s DETAILED %s\n' "$s" "$(cat totals.$s.txt)" > report.$s.txt
        #@end
    else
        #@7,5,brief,,in,sbatch -p short -c 1 --mem 100M -t 5:00
        printf '%s brief\n' "$s" > report.$s.txt
        #@end
    fi
done
