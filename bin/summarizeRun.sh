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
echo Category,Used,Wasted > dataMem.csv
                                                                            
ls *.out | sort -n | xargs -d '\n' grep ^dataToPlot | awk -F, '{printf "%s-%s-%s,%s,%s\n", substr($15,1,index($15,".")-1), substr($2, length($2)-3), substr($10,1,3),  $8 + $17,   $6-$8-$17}' | sed 's/,-.*/,0/' >> dataMem.csv

#ls * | sort -n | xargs -d '\n' grep hello

#awk -F, '{sum=$4+$5; substr($3,1,3); printf "%s%s,%s,%s,%s,%s\n", $2, $3, $4, $5, substr($3,1,3), sum}' input.csv > output.csv


gnuplot -e "set datafile separator ','; set style data histogram; set style histogram rowstacked gap 2; set style fill solid border rgb 'black'; set xtics rotate by -45; set terminal png size 800,600; set output 'barchartMem.png'; set title 'Job vs. Memmory'; set ylabel 'Memory (MegaBytes)'; plot 'dataMem.csv' using 2:xtic(1) title 'Used' lc rgb 'green', '' using 3:xtic(1) title 'Wasted' lc rgb 'red'"


# gnuplot -e "set datafile separator ','; set style data histogram; set style histogram rowstacked gap 1; set style fill solid border rgb 'black'; set xtics ('A' 0, 'B' 1, 'C' 2, 'D' 3, 'E' 4) rotate by -45; set terminal png size 800,600; set output 'barchart.png'; set title 'Sample vs. Time'; plot 'data.csv' using 2:xtic(1) title 'Value1' lc rgb 'green', '' using 3:xtic(1) title 'Value2' lc rgb 'red'"


echo To see the plot: 
echo display log/barchartMem.png 

echo Category,Used,Wasted > dataTime.csv

ls *.out | sort -n | xargs -d '\n' grep ^dataToPlot | awk -F, '{printf "%s-%s-%s,%s,%s\n", substr($15,1,index($15,".")-1), substr($2, length($2)-3), substr($10,1,3),  $9,   $7-$9}' >> dataTime.csv

gnuplot -e "set datafile separator ','; set style data histogram; set style histogram rowstacked gap 2; set style fill solid border rgb 'black'; set xtics rotate by -45; set terminal png size 800,600; set output 'barchartTime.png'; set title 'Job vs. Time'; set ylabel 'Time (Mins)'; plot 'dataTime.csv' using 2:xtic(1) title 'Used' lc rgb 'green', '' using 3:xtic(1) title 'Wasted' lc rgb 'red'"

echo To see the plot: 
echo display log/barchartTime.png 

# check running jobs 
checkJobsSlurm  allJobs.txt
