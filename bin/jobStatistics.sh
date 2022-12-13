#!/bin/bash

echo Running
echo $0 $@ 

# usage: command <software> <ref> <numberOfJobRecordsToUse>

[[ -z "$1" || -z "$2" ]] && echo Usage: checkJobRecord.sh software reference && exit 1   

# if not rsynced today, sync some jobRecords from other users 
#if [ ! -f ~/smartSlurm/stats/jobRecord.txt ]; then
#    echo File not exist:  ~/smartSlurm/stats/jobRecord.txt, sync some some data ...
#    cat ~/smartSlurm/stats/myJobRecord.txt > ~/smartSlurm/stats/jobRecord.txt
#    if [ ! -f ~/smartSlurm/stats/jobRecord.txt ]; then
#      touch ~/smartSlurm/stats/jobRecord.txt #
#    fi  

#else
#    # file modified in last 2 minutes
    [ ! -z "`find ~/smartSlurm/jobRecord.txt -mmin -2`" ] && echo jobRecord.txt synced within 20 hour. No need to re-sync || cat $HOME/smartSlurm/myJobRecord.txt > ~/smartSlurm/jobRecord.txt  
#fi

OUT="$(mktemp -d)"

#filter by software and reference
grep COMPLETED ~/smartSlurm/jobRecord.txt | awk -F"," -v a=$1 -v b=$2 '{ if($12 == a && $13 == b) {print $2, $7 }}' | sort -r  -k1,1 -k2,2 | sort -u -k1,1 > $OUT/mem.txt
grep COMPLETED ~/smartSlurm/jobRecord.txt | awk -F"," -v a=$1 -v b=$2 '{ if($12 == a && $13 == b) {print $2, $8 }}' | sort -r  -k1,1 -k2,2 | sort -u -k1,1 > $OUT/time.txt

#awk -v a=$1 -v b=$2 '{ if($2 == a && $3 == b) {print $5, $13 }}' ~/smartSlurm/stats/jobRecord.txt ~/smartSlurm/stats/myJobRecord.txt | grep COMPLETED | uniq > $OUT/time.txt

echo "Got mem data from jobRecord.txt (content of mem.txt):"
cat $OUT/mem.txt 

echo "Got time data from jobRecord.txt (content of time.txt):"
cat $OUT/time.txt 

if [[ $(wc -l <$OUT/mem.txt) -le $3 ]]; then
    echo There is less than $3 records. No way to fit a curve. Exiting...
    exit 
fi

cd $OUT

# make plot and calculate statistics
gnuplot -e 'set term pdf; set output "mem.pdf"; set title "Input Size vs. Run Time" font "Helvetica Bold,18"; set xlabel "Input Size(K)"; set ylabel "Run Time(Min)"; f(x)=a*x+b; fit f(x) "mem.txt" u 1:2 via a, b; t(a,b)=sprintf("f(x) = %.2fx + %.2f", a, b); plot "mem.txt" u 1:2,f(x) t t(a,b); print "Finala=", a; print "Finalb=",b; stats "mem.txt" u 1 ' 2>&1 | grep 'Final\| M' | awk 'NF<4{print $1, $2}' |sed 's/:/=/' | sed 's/ //g' > mem.stat.txt

# make plot and calculate statistics
gnuplot -e 'set term pdf; set output "time.pdf"; set title "Input Size vs. Memory Usage" font "Helvetica Bold,18"; set xlabel "Input Size(K)"; set ylabel "Memory(M)"; f(x)=a*x+b; fit f(x) "time.txt" u 1:2 via a, b; t(a,b)=sprintf("f(x) = %.2fx + %.2f", a, b); plot "time.txt" u 1:2,f(x) t t(a,b); print "Finala=", a; print "Finalb=",b; stats "time.txt" u 1 ' 2>&1 | grep 'Final\| M' | awk 'NF<4{print $1, $2}' |sed 's/:/=/' | sed 's/ //g' > time.stat.txt

# after this file is created, don't need to calculate stats anymore 
if [[ $(wc -l <mem.txt) -ge $3 ]]; then 
  mv $OUT/mem.stat.txt ~/smartSlurm/stats/$1.$2.mem.stat.final 
  mv $OUT/time.stat.txt ~/smartSlurm/stats/$1.$2.time.stat.final  
  echo There are more than $3 jobs already run for this software, statics is ready for current job: 
  echo Memeory statisics:
  cat ~/smartSlurm/stats/$1.$2.mem.stat.final
  echo
  echo Time statistics:
  cat ~/smartSlurm/stats/$1.$2.time.stat.final
else
  mv $OUT/mem.stat.txt ~/smartSlurm/stats/$1.$2.mem.stat 
  mv $OUT/time.stat.txt ~/smartSlurm/stats/$1.$2.time.stat
  echo There are less than $3 jobs already run for this software, statics is not ready to use yet: 
  echo 
  echo Memeory statisics:
  echo "inputsize mem(M)"
  cat ~/smartSlurm/stats/$1.$2.mem.stat
  echo Time statistics:
  echo "inputsize time(minute)"
  cat ~/smartSlurm/stats/$1.$2.time.stat
fi  

#mv $OUT/mem.pdf ~/smartSlurm/stats/$1.$2.mem.pdf 
mv $OUT/mem.txt ~/smartSlurm/stats/$1.$2.mem.txt
#mv $OUT/time.pdf ~/smartSlurm/stats/$1.$2.time.pdf 
mv $OUT/time.txt ~/smartSlurm/stats/$1.$2.time.txt

convert $OUT/mem.pdf -background White -flatten ~/smartSlurm/stats/$1.$2.mem.pdf
convert $OUT/time.pdf -background White -flatten ~/smartSlurm/stats/$1.$2.time.pdf
echo
echo You can see the plot using commands:
echo display ~/smartSlurm/stats/$1.$2.mem.pdf
echo display ~/smartSlurm/stats/$1.$2.time.pdf

rm -r $OUT 2>/dev/null

# Finala=0.03
# Finalb=5.0
# Mean=250.0000
# Minimum=200.0000
# Maximum=300.0000
# Median=250.0000
#rm -r $OUT 
echo got files in ~/smartSlurm/stats:  
ls -lrt ~/smartSlurm/stats/


