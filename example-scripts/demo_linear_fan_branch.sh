#!/bin/sh
# =============================================================================
# rAP EXAMPLE — "linear -> fan-out/in -> branch"  (best-practice template)
# =============================================================================
# A single pipeline that walks through the three shapes you'll actually use:
#   * a short LINEAR chain           (steps 1-2)
#   * a FAN-OUT / FAN-IN pair         (steps 3-5: split, process chunks, merge)
#   * an IF/ELSE that submits one of two ALTERNATIVE steps (steps 6 | 7)
#
# RUN IT
#   runAsPipeline --script "demo_linear_fan_branch.sh demo_manifest.txt" \
#                 --sbatch-options "sbatch -p short -c 1 --mem 100M -t 5:00" \
#                 --tmp noTmp
#   # --mode dryrun previews the conversion without submitting.
#
# MANIFEST  (one sample per line:  <id> <phrase words...>)
# =============================================================================
manifest=$1
[ -z "$manifest" ] && { echo "usage: $0 <manifest>"; exit 1; }

# --- prep on the SUBMIT host: one input file per sample ----------------------
samples=$(cut -d' ' -f1 "$manifest")
while read sid rest; do echo "$rest" > doc.$sid.txt; done < "$manifest"

for s in $samples; do
    in=doc.$s.txt      # this sample's input; named in each marker's "inputs" field

    # ---- LINEAR chain -------------------------------------------------------
    #@1,0,normalize,,in,sbatch -p short -c 1 --mem 100M -t 5:00
    tr 'A-Z' 'a-z' < $in | tr -s ' ' | sed 's/^ *//; s/ *$//' > norm.$s.txt
    #@end

    #@2,1,tokenize,,in,sbatch -p short -c 1 --mem 100M -t 5:00
    tr ' ' '\n' < norm.$s.txt | sed '/^$/d' > tokens.$s.txt
    #@end

    # ---- FAN-OUT prep: round-robin tokens into exactly 3 chunk files --------
    # NOTE: use explicit filenames, NOT a for-loop inside the block: a loop var
    #       assigned *inside* a block is not recognized as node-scope.
    #@3,2,split,,in,sbatch -p short -c 1 --mem 100M -t 5:00
    touch chunk.$s.1 chunk.$s.2 chunk.$s.3
    awk -v s="$s" '{ f = "chunk." s "." (((NR-1)%3)+1); print > f }' tokens.$s.txt
    #@end

    # ---- FAN-OUT: the nested  for c  loop submits one job per chunk. --------
    # rAP includes BOTH enclosing loop vars in the job flag (…split.$s.$c), so
    # the three chunk jobs never collide.
    for c in 1 2 3; do
        #@4,3,countchunk,,in,sbatch -p short -c 1 --mem 100M -t 5:00
        w=$(wc -w < chunk.$s.$c)
        printf 'chunk %s words=%s\n' "$c" "$w" > cc.$s.$c.txt
        #@end
    done

    # ---- FAN-IN: step 5 depends on step 4, so it waits for ALL chunk jobs. --
    #@5,4,merge,,in,sbatch -p short -c 1 --mem 100M -t 5:00
    tw=$(cat cc.$s.*.txt | grep -o 'words=[0-9]*' | cut -d= -f2 | awk '{ s+=\$1 } END{ print s+0 }')
    printf 'sample=%s total_words=%s\n' "$s" "$tw" > totals.$s.txt
    #@end

    # ---- IF/ELSE: pick ONE of two alternative terminal steps ----------------
    # IMPORTANT: the condition must be knowable at SUBMIT time. The driver
    # submits every job before any runs, so you CANNOT branch on a job's output
    # (e.g. totals.$s.txt does not exist yet). Branch on the input instead.
    # IMPORTANT: the two branches use DISTINCT step numbers (6 and 7). Reusing
    # the same step number across if/else makes rAP drop the second block.
    words=$(wc -w < $in)
    if [ "$words" -ge 5 ]; then
        #@6,5,detailed,,in,sbatch -p short -c 1 --mem 100M -t 5:00
        printf '%s: DETAILED (%s)\n' "$s" "$(cat totals.$s.txt)" > report.$s.txt
        #@end
    else
        #@7,5,brief,,in,sbatch -p short -c 1 --mem 100M -t 5:00
        printf '%s: brief note, short input\n' "$s" > report.$s.txt
        #@end
    fi
done
