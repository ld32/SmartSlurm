#!/bin/sh
number=$1
for s in {1..2}; do
    input=numbers$s.txt

    #@1,0,analyze,,input,sbatch -p short -c 1 --mem 2G -t 50:0
    # a multi-line body with a blank line, quotes, and a heredoc -- all of which
    # break the legacy flatten path but are fine as a script on disk.
    label="sample $s"
    grep "$number" $input > $number.$s.txt

    cat <<INNER >> $number.$s.log
    processed $label with quotes: 'single' "double"
INNER
    #@end

done

#@2,1,merge,,,sbatch -p short -c 1 --mem 4G -t 50:0
cat $number.*.txt > all$number.txt
#@end
