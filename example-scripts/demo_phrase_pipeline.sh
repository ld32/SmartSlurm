#!/bin/sh
# SmartSlurm DEMO pipeline — a gentle "learn the #@ syntax" starting point.
# A minimal LINEAR chain in the recommended explicit-#@end style (scripts-on-disk):
#   normalize -> tokenize -> measure -> report   (per sample, inside the loop)
#   summary                                       (experiment-level fan-in)
#
# Run it (from its own directory):
#   runAsPipeline --script "demo_phrase_pipeline.sh demo_phrase_manifest.txt" \
#                 --sbatch-options "sbatch -p short -c 1 --mem 100M -t 5:00" --tmp noTmp
#   # --mode dryrun previews the conversion without submitting; --lint checks structure only.
#
# manifest: one sample per line ->  <sampleID> <phrase words...>
#@pipeline phraseDemo
manifest=$1
[ -z "$manifest" ] && { echo "usage: $0 <manifest>"; exit 1; }

# prep (runs on submit host): one phrase file per sample
samples=$(cut -d' ' -f1 "$manifest")
while read sid rest; do echo "$rest" > phrase.$sid.txt; done < "$manifest"

for s in $samples; do
    in=phrase.$s.txt

    #@1,0,normalize,,in,sbatch -p short -c 1 --mem 100M -t 5:00
    # lowercase, collapse spaces, trim ends
    tr 'A-Z' 'a-z' < $in | tr -s ' ' | sed 's/^ *//; s/ *$//' > norm.$s.txt
    #@end

    #@2,1,tokenize,,in,sbatch -p short -c 1 --mem 100M -t 5:00
    # one word per line, drop blanks
    tr ' ' '\n' < norm.$s.txt | sed '/^$/d' > tokens.$s.txt
    #@end

    #@3,2,measure,,in,sbatch -p short -c 1 --mem 100M -t 5:00
    # "word <tab> length", longest first (awk + quotes survive scripts-on-disk)
    awk '{ print length(\$1)"\t"\$1 }' tokens.$s.txt | sort -k1,1nr -k2,2 > measured.$s.txt
    #@end

    #@4,3,report,,in,sbatch -p short -c 1 --mem 100M -t 5:00
    # in-block vars ($n,$u,$long) expand on the node; $s is baked per job
    n=$(wc -l < tokens.$s.txt)
    u=$(sort -u tokens.$s.txt | wc -l)
    long=$(head -1 measured.$s.txt | cut -f2)
    printf '%s\twords=%s\tunique=%s\tlongest=%s\n' "$s" "$n" "$u" "$long" > report.$s.txt
    #@end

done

#@5,4,summary,,,sbatch -p short -c 1 --mem 100M -t 5:00
# experiment-level fan-in across all samples
{ printf 'sample\treport\n'; cat report.*.txt; } > summary.tsv
cat summary.tsv
#@end
