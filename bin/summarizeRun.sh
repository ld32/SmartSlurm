#!/bin/sh
Usage="Usage: $0 \nThis script will go through job name list in log/allJobs.txt to see if the jobs finish successfully or not."

set -x

IFS=$'\n'

out=`squeue -u $USER -t PD,R --noheader -o "%.18i-%t"`

lines=`tail -n +2 log/allJobs.txt | awk 'NF>2{print $1, $2, $3}'`

current=0; succ=0; fail=0; running=0; pending=0; unknown=0
toSend="Summery for jobs in log/allJobs.txt:"
for line in $lines; do
    if [ ! -z "${line/ /}" ]; then
        id=${line%% *}; name=${line##* }


        if [ -f log/$name.success ]; then
            toSend="$toSend\nDone: $line"
            succ=$((succ + 1))
        elif [ -f log/$name.failed ]; then
            toSend="$toSend\nFailed $line"
            fail=$((fail + 1))
        elif [[ "$out" == *$id-R* ]]; then # && [[ "$id" != "$SLURM_JOBID" ]]; then
            toSend="$toSend\nRunn: $line"
            running=$((running + 1))
        elif [[ "$out" == *$id-P* ]]; then # && [[ "$id" != "$SLURM_JOBID" ]]; then
            toSend="$toSend\nRunn: $line"
            pending=$((pending + 1))
        else
            toSend="$toSend\nUnkn: $line"
            unknown=$((unknown + 1))
        fi
        [ "$id" == "$SLURM_JOBID" ] && current=$((succ + fail + running + pending + unknown))
    fi
done

total=$((succ + fail + running + pending + unknown))
s="$current/$total Succ:$succ/$total Running:$running/$total Pending:$pending/$total Fail:$fail/$total Unknown:$unknown/$total"

echo -e "$s" > log/summary

# if [ $((running + pending)) -le 5 ]; then
#     echo -e "$toSend" | mail -s $s $USER
#     [ "$USER" != ld32 ] && echo -e "$toSend" | mail -s $s ld32
# fi

set +x

cd log

rm barchartMem.png  barchartTime.png 2>/dev/null
echo Category,Used,Wasted,Saved2,default,Saved1 > dataMem.csv

ls *.out | sort -n | xargs -d '\n' grep ^dataToPlot | awk -F, '{printf "%s-%s-%s,%s,%s,%s,%s\n", substr($15,1,index($15,".")-1), substr($2, length($2)-3), substr($10,1,3),  $8 + $17 *2,   $6-$8-$17 *2, $4-$6, $4}' | sed s/-COM//g | sed s/-OO/-/g >> dataMem.csv

# if less than 0, change to zeor
awk -F',' -v OFS=',' '{ for (i=1; i<=NF; i++) if ($i < 0) $i = 0; print }' dataMem.csv > output.csv

awk -F, -v OFS=',' -v max=$(awk -F, 'BEGIN {max=0} {if (NR!=1 && $5>max) max=$5} END {print max}' output.csv) '{if(NR==1) print $0; else {diff=max-$5; print $0 "," diff "," max}}' output.csv > dataMem.csv #> output.csv

# all job memory
gnuplot -e "set key outside; set key reverse; set key invert; set datafile separator ','; set style data histogram; set style histogram rowstacked gap 2; set style fill solid border rgb 'black'; set xtics rotate by -45; set terminal png size 800,600; set output 'barchartMem.png'; set title 'Job vs. Memmory'; set ylabel 'Memory (MegaBytes)'; plot 'dataMem.csv' using 2:xtic(1) title 'Used' lc rgb 'green', '' using 3:xtic(1) title 'Wasted' lc rgb 'red', '' using 4:xtic(1) title 'Saved2' lc rgb 'yellow', '' using 6:xtic(1) title 'Saved1' lc rgb 'pink'"

echo To see the plot:
echo display log/barchartMem.png

echo Category,Used,Wasted,default,Saved > dataTime.csv

ls *.out | sort -n | xargs -d '\n' grep ^dataToPlot | awk -F, '{printf "%s-%s-%s,%s,%s,%s,%s\n", substr($15,1,index($15,".")-1), substr($2, length($2)-3), substr($10,1,3),  $9*$16,   ($7-$9)*$16, ($5-$7)*$16, $5*$16}' | sed s/-COM//g | sed s/-OO/-/g >> dataTime.csv

awk -F',' -v OFS=',' '{ for (i=1; i<=NF; i++) if ($i < 0) $i = 0; print }' dataTime.csv > output.csv #> dataMem.csv

awk -F, -v OFS=',' -v max=$(awk -F, 'BEGIN {max=0} {if (NR!=1 && $5>max) max=$5} END {print max}' output.csv) '{if(NR==1) print $0; else {diff=max-$5; print $0 "," diff "," max}}' output.csv > dataTime.csv

# all job time
gnuplot -e "set key outside; set key reverse; set key invert; set datafile separator ','; set style data histogram; set style histogram rowstacked gap 2; set style fill solid border rgb 'black'; set xtics rotate by -45; set terminal png size 800,600; set output 'barchartTime.png'; set title 'Job vs. Time'; set ylabel 'Time (Mins)'; plot 'dataTime.csv' using 2:xtic(1) title 'Used' lc rgb 'green', '' using 3:xtic(1) title 'Wasted' lc rgb 'red', '' using 4:xtic(1) title 'Saved' lc rgb 'yellow'" #", '' using 6:xtic(1) title 'Saved' lc rgb 'pink'"

echo To see the plot:
echo display log/barchartTime.png

# time vs. memory for current job
gnuplot -e "set key outside; set key reverse; set key invert; set datafile separator ' '; set style data histogram; set style histogram rowstacked gap 2; set style fill solid border rgb 'black'; set xtics rotate by -45; set terminal png size 800,600; set output 'job_$SLURM_JOBID.mem.png'; set title 'Time vs. Mem for job $SLURM_JOBID'; set xlabel 'Time (Mins)'; set ylabel 'Mem (M)'; plot 'job_$SLURM_JOBID.memCPU.txt' using 2:xtic(1) title 'Used' lc rgb 'green', '' using 3:xtic(1) title 'Wasted' lc rgb 'red', '' using 4:xtic(1) title 'Saved' lc rgb 'yellow'"

# time vs. CPU usage for current job
gnuplot -e "set key outside; set key reverse; set key invert; set datafile separator ' '; set style data histogram; set style histogram rowstacked gap 2; set style fill solid border rgb 'black'; set xtics rotate by -45; set terminal png size 800,600; set output 'job_$SLURM_JOBID.cpu.png'; set title 'Time vs. CPU Usage for job $SLURM_JOBID'; set xlabel 'Time (Mins)'; set ylabel 'CPU Usage (%)'; plot 'job_$SLURM_JOBID.memCPU.txt' using 5:xtic(1) title 'Used' lc rgb 'green'"

