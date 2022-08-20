#!/bin/sh

Usage="Usage: $0 full_path_to_flag_folder  flag(job name or job id) [remove_success_flag_file] \n  Note: this script will go through job id list file, find the downstream jobs, and kill them if the current job failed or remove success flag file. "

[ -z "$2" ] && { echo -e "job name can not be empty! \n$Usage"; exit; }
               
[ -f $1/allJobs.txt.first ] || { echo Job record file not exist: $1/allJobs.txt.first; exit 1; } 

text=`cat $1/allJobs.txt.first`

job=$2

IFS=$' '; 

# check the third column for the job name, then find the the job id in column 1
id=`echo $text | awk '{if ($3 ~ /'"$job/"') print $1,$3;}'`
#echo job id $id

[ -z "$id" ] && { echo -e "job id for $job can not found! \n$Usage"; exit 1; }

echo 
declare -A seen
path=$1
function findID {
    #echo function start for $1
    ids="$1" 
    [ -z "$ids" ] &&  return

    IFS=$'\n';   
    for k in $ids; do
    	i=${k% *}; j=${k#* };
        #echo working on $k ,   get  i = $i , j = $j 
        if [ ! "${seen[$i]}" ]; then
              seen[$i]=1
              echo rm $path/$j.success; rm $path/$j.* 2>/dev/null 
        else
              echo job killed before
        fi
   
        IFS=$' ';             # check second column for the job id
    	ids=`echo $text | awk '{if ($2 ~ /'"$i/"') print $1,$3;}'`
        #echo working on $i find: $ids | tr '\n' ' '; echo

        findID "$ids"
    done
}

findID "$id"
echo Will re-run the down stream steps even if they are done before.
