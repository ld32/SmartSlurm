#!/bin/sh
Usage="Usage: $0 \nThis script will go through job name list in log/allJobs.txt to see if the jobs finish successfully or not."

#set -x

echo According to log/allJobs.txt:
names=`tail -n -1 log/allJobs.txt | awk '{print $3}' | tr "\n" " "`
for name in $names; do
   [ -f log/$name.success ] && echo Success! $name  || { echo "Failed! $name #############" && [ -f log/$name.err ] && echo -e "Here is the content of log/$name.err :\n" &&  tail -n 10 log/$name.err; }
   #fi
done

if [ -f log/allJobs.txt.firstx ]; then # || { echo Job list file not exist: log/allJobs.txt.first; exit 1; }
   echo According to log/allJobs.txt.first:
   names=`tail -n -1 log/allJobs.txt.first | awk '{print $3}' | tr "\n" " "`
   for name in $names; do
         if [ -f log/$name.start ]; then
            [ -f log/$name.success ] && echo Success! $name  || { echo "Failed! $name #############" && [ -f log/$name.err ] && echo -e "Here is the content of log/$name.err :\n" &&  tail -n 10 log/$name.err; }
         else
            echo Didn not run! $name
         fi
   done
fi

cd log

rm barchartMem.png  barchartTime.png 2>/dev/null
echo Category,Used,Wasted,Saved2,default,Saved1 > dataMem.csv

ls *.out | sort -n | xargs -d '\n' grep ^dataToPlot | awk -F, '{printf "%s-%s-%s,%s,%s,%s,%s\n", substr($15,1,index($15,".")-1), substr($2, length($2)-3), substr($10,1,3),  $8 + $17 *2,   $6-$8-$17 *2, $4-$6, $4}' | sed s/-COM//g | sed s/-OO/-/g >> dataMem.csv

# if less than 0, change to zeor
awk -F',' -v OFS=',' '{ for (i=1; i<=NF; i++) if ($i < 0) $i = 0; print }' dataMem.csv > output.csv #> dataMem.csv

awk -F, -v OFS=',' -v max=$(awk -F, 'BEGIN {max=0} {if (NR!=1 && $5>max) max=$5} END {print max}' output.csv) '{if(NR==1) print $0; else {diff=max-$5; print $0 "," diff "," max}}' output.csv > dataMem.csv #> output.csv

gnuplot -e "set key outside; set key reverse; set key invert; set datafile separator ','; set style data histogram; set style histogram rowstacked gap 2; set style fill solid border rgb 'black'; set xtics rotate by -45; set terminal png size 800,600; set output 'barchartMem.png'; set title 'Job vs. Memmory'; set ylabel 'Memory (MegaBytes)'; plot 'dataMem.csv' using 2:xtic(1) title 'Used' lc rgb 'green', '' using 3:xtic(1) title 'Wasted' lc rgb 'red', '' using 4:xtic(1) title 'Saved2' lc rgb 'yellow', '' using 6:xtic(1) title 'Saved1' lc rgb 'pink'"

echo To see the plot:
echo display log/barchartMem.png

echo Category,Used,Wasted,default,Saved > dataTime.csv

ls *.out | sort -n | xargs -d '\n' grep ^dataToPlot | awk -F, '{printf "%s-%s-%s,%s,%s,%s,%s\n", substr($15,1,index($15,".")-1), substr($2, length($2)-3), substr($10,1,3),  $9*$16,   ($7-$9)*$16, ($5-$7)*$16, $5*$16}' | sed s/-COM//g | sed s/-OO/-/g >> dataTime.csv

awk -F',' -v OFS=',' '{ for (i=1; i<=NF; i++) if ($i < 0) $i = 0; print }' dataTime.csv > output.csv #> dataMem.csv


awk -F, -v OFS=',' -v max=$(awk -F, 'BEGIN {max=0} {if (NR!=1 && $5>max) max=$5} END {print max}' output.csv) '{if(NR==1) print $0; else {diff=max-$5; print $0 "," diff "," max}}' output.csv > dataTime.csv

gnuplot -e "set key outside; set key reverse; set key invert; set datafile separator ','; set style data histogram; set style histogram rowstacked gap 2; set style fill solid border rgb 'black'; set xtics rotate by -45; set terminal png size 800,600; set output 'barchartTime.png'; set title 'Job vs. CPUTime'; set ylabel 'CPUTime (CPUxMins)'; plot 'dataTime.csv' using 2:xtic(1) title 'Used' lc rgb 'green', '' using 3:xtic(1) title 'Wasted' lc rgb 'red', '' using 4:xtic(1) title 'Saved' lc rgb 'yellow', '' using 6:xtic(1) title 'Saved' lc rgb 'pink'"

echo To see the plot:
echo display log/barchartTime.png

gnuplot -e "set key outside; set key reverse; set key invert; set datafile separator ' '; set style data histogram; set style histogram rowstacked gap 2; set style fill solid border rgb 'black'; set xtics rotate by -45; set terminal png size 800,600; set output 'job_$SLURM_JOBID.mem.png'; set title 'Time vs. Mem for job $SLURM_JOBID'; set xlabel 'Time (Mins)'; set ylabel 'Mem (M)'; plot 'job_$SLURM_JOBID.mem.txt' using 2:xtic(1) title 'Used' lc rgb 'green', '' using 3:xtic(1) title 'Wasted' lc rgb 'red', '' using 4:xtic(1) title 'Saved' lc rgb 'yellow'"
# check running jobs
# checkJobsSlurm  allJobs.txt
