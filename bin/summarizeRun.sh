#!/bin/sh

Usage="Usage: $0 [jobRecordFile, optional, for example: log/allJobs.txt.old]\nNote: this script will go through job name list in a file (eg: log/alljobs), to see if the jobs finish successfully or not." 

#set -x 

if [ ! -z "$1" ]; then 
    [ -f "$1" ] || { echo $Usage; exit 1; }
     
    echo According to $1:
    names=`tail -n -1 $1 | awk '{print $3}' | tr "\n" " "`
    for name in $names; do       
       if [ -f log/$name.start ]; then 
           [ -f log/$name.success ] && echo Success! $name  || { echo "Failed! $name #############" && [ -f log/$name.err ] && echo -e "Here is the content of log/$name.err :\n" &&  tail -n 10 log/$name.err; }
       else 
          echo Didn not run! $name 
       fi   
    done
    #exit 
fi       

if [ -f log/allJobs.txt.first ]; then # || { echo Job list file not exist: log/allJobs.txt.first; exit 1; }

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


# echo Category,Value1,Value2 > log/dataMem.csv
#                                           #                   software, status, used, wasted
# grep ^dataToPlot log/*.out | awk -F, '{printf "%s-%s,%s,%s\n", $13,     $10,    $8,   $6-$8}'> log/dataMem.csv

# gnuplot -e "set datafile separator ','; set style data histogram; set style histogram rowstacked gap 2; set style fill solid border rgb 'black'; set xtics rotate by -45; set terminal png size 800,600; set output 'log/barchartMem.png'; plot 'log/dataMem.csv' using 2:xtic(1) title 'Value1' lc rgb 'green', '' using 3:xtic(1) title 'Value2' lc rgb 'red'"



cd log 



rm barchartMem.png  barchartTime.png 
echo Category,Used,Wasted,default,Saved > dataMem.csv
                                                                            
ls *.out | sort -n | xargs -d '\n' grep ^dataToPlot | awk -F, '{printf "%s-%s-%s,%s,%s,%s\n", substr($15,1,index($15,".")-1), substr($2, length($2)-3), substr($10,1,3),  $8 + $17,   $6-$8-$17, $4}' | sed 's/,-.*/,0/' | sed s/-COM//g | sed s/-OO/-/g >> dataMem.csv



#awk -F, -v OFS=',' 'NR==1 {print $0 ",Saved"; next} {if ($4>max) max=$4; print $0} END {print "," max}' input.csv > output.csv


awk -F, -v OFS=',' -v max=$(awk -F, 'BEGIN {max=0} {if (NR!=1 && $4>max) max=$4} END {print max}' dataMem.csv) '{if(NR==1) print $0; else {diff=max-$2-$3; if(diff<0) diff=0; print $0 "," diff}}' dataMem.csv > output.csv



#awk -F, -v OFS=',' '{diff=$2-$1; if(diff<0) diff=0; print $0, diff}' input.csv > output.csv


mv output.csv dataMem.csv

#ls * | sort -n | xargs -d '\n' grep hello

#awk -F, '{sum=$4+$5; substr($3,1,3); printf "%s%s,%s,%s,%s,%s\n", $2, $3, $4, $5, substr($3,1,3), sum}' input.csv > output.csv


gnuplot -e "set datafile separator ','; set style data histogram; set style histogram rowstacked gap 2; set style fill solid border rgb 'black'; set xtics rotate by -45; set terminal png size 800,600; set output 'barchartMem.png'; set title 'Job vs. Memmory'; set ylabel 'Memory (MegaBytes)'; plot 'dataMem.csv' using 2:xtic(1) title 'Used' lc rgb 'green', '' using 3:xtic(1) title 'Wasted' lc rgb 'red', '' using 5:xtic(1) title 'Saved' lc rgb 'yellow'"


# gnuplot -e "set datafile separator ','; set style data histogram; set style histogram rowstacked gap 1; set style fill solid border rgb 'black'; set xtics ('A' 0, 'B' 1, 'C' 2, 'D' 3, 'E' 4) rotate by -45; set terminal png size 800,600; set output 'barchart.png'; set title 'Sample vs. Time'; plot 'data.csv' using 2:xtic(1) title 'Value1' lc rgb 'green', '' using 3:xtic(1) title 'Value2' lc rgb 'red'"


echo To see the plot: 
echo display log/barchartMem.png 

echo Category,Used,Wasted,default,Saved > dataTime.csv

ls *.out | sort -n | xargs -d '\n' grep ^dataToPlot | awk -F, '{printf "%s-%s-%s,%s,%s,%s\n", substr($15,1,index($15,".")-1), substr($2, length($2)-3), substr($10,1,3),  $9,   $7-$9, $5}' | sed s/-COM//g | sed s/-OO/-/g >> dataTime.csv

awk -F, -v OFS=',' '{if (NF==1) {print $0} else {diff=$4-$2-$3; if(diff<0) diff=0; print $0 "," diff}}' dataTime.csv > output.csv

mv output.csv dataTime.csv

gnuplot -e "set datafile separator ','; set style data histogram; set style histogram rowstacked gap 2; set style fill solid border rgb 'black'; set xtics rotate by -45; set terminal png size 800,600; set output 'barchartTime.png'; set title 'Job vs. Time'; set ylabel 'Time (Mins)'; plot 'dataTime.csv' using 2:xtic(1) title 'Used' lc rgb 'green', '' using 3:xtic(1) title 'Wasted' lc rgb 'red', '' using 5:xtic(1) title 'Saved' lc rgb 'yellow'"

echo To see the plot: 
echo display log/barchartTime.png 

#if [ -f /tmp/job_$SLURM_JOBID.mem.txt ]; then
   #gnuplot -e "set datafile separator ' '; set style data histogram; set style histogram rowstacked gap 2; set style fill solid border rgb 'black'; set xtics rotate by -45; set terminal png size 800,600; set output 'log/job_$SLURM_JOBID.mem.png'; set title 'Time vs. Memory'; set ylabel 'Mem (M)'; set xlabel 'Time (m)';  plot '/tmp/job_$SLURM_JOBID.mem.txt' using 2:xtic(1) lc rgb 'green'"

gnuplot -e "set datafile separator ' '; set style data histogram; set style histogram rowstacked gap 2; set style fill solid border rgb 'black'; set xtics rotate by -45; set terminal png size 800,600; set output 'job_$SLURM_JOBID.mem.png'; set title 'Time vs. Mem for job $SLURM_JOBID'; set xlabel 'Time (Mins)'; set ylabel 'Mem (M)'; plot 'job_$SLURM_JOBID.mem.txt' using 2:xtic(1) title 'Used' lc rgb 'green', '' using 3:xtic(1) title 'Wasted' lc rgb 'red', '' using 4:xtic(1) title 'Saved' lc rgb 'yellow'"

#fi

# check running jobs 
# checkJobsSlurm  allJobs.txt
