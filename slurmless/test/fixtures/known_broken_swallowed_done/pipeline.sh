#!/bin/sh
# The shipped-example foot-gun: the outer loop's `done` directly abuts the step-2
# block with NO blank line, so the blank-line-terminated parser swallows `done`
# into the wrapped command and the outer for-loop is left unterminated.
# Conversion therefore emits a driver that is INVALID bash. The M4 `#@ … #@end`
# closure fixes this. Kept as a known-issue regression anchor.
number=$1
for j in {1..2}; do
    input=numbers$j.txt
    #@1,0,findNumber,,input,sbatch -p short -c 1 --mem 2G -t 50:0
    findNumber.sh $number $input > $number.$j.txt

    #@2,1,mergeNumber,,,sbatch -p short -c 1 --mem 4G -t 50:0
    cat $number.*.txt > all$number.txt
done
