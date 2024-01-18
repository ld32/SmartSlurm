#!/bin/bash
set -x 
set -e



#I think you can run tests on /n/groups/htem/tier2/cb3 (but please don't actually archive it); /n/groups/htem/tier2/cb3/sections/170218121547_cb3_0293 is just one section out of ~500, and our more recent datasets are about 1000-4000 sections long. Actually the datasets in /n/groups/htem/tier2 are the ones we prepared to go on tier2 but did not happen because of IT difficulties so they could also make good case studies.


function archiveFiles() {
    files=${1#\\n}                           # remove leading /n
    
    item=${files%%\\n*}-${files##*\\n}      # firstFile-lastFile

    [ -f "$2/$item.tar" ] && echo done earlier && return  
  
    while true; do
         local jobID=0
         if { set -C; 2>/dev/null >$dFolderTmp/lock.0; }; then
            trap "rm -f $dFolderTmp/lock.0" EXIT
            for i in $(seq 1 $core); do 
                if [ ! -f $dFolderTmp/lock.$i ]; then 
                    touch $dFolderTmp/lock.$i
                    jobID=$i
                    break
                fi
            done
            rm -f $dFolderTmp/lock.0
        fi
        echo sleeping
        [ "$jobID" -eq "0" ] && sleep 10 || break
    done
    echo -e "$files" > "$2/$item.list.txt"
    local tmp=`mktemp` &&
    tar --create --preserve-permissions --file "$tmp" -T $2/$item.list.txt && 
    md5sum  $tmp > $tmp.md5sum &&
    mv "$tmp"  "$2/$item.tar" && mv "$tmp.md5sum"  "$2" &&
    rm -f  $dFolderTmp/lock.$jobID &&
    echo $(date) job: $jobID | tee -a  $dFolder/archive.log  &&
    echo cd `pwd` | tee -a $dFolder/archive.log  &&
    echo tar --create --preserve-permissions --file "$2/$item.tar" -T $2/$item.list.txt | tee -a $dFolder/archive.log &
    
}

function archiveFolder() {
    mkdir -p "$2"

    files=""    

    totalSize="0"

    cd $1

    for line in `find -L . -maxdepth 1 -mindepth 1 -type d -printf "%f\n" | sort -n`; do
        (archiveFolder $line "$2/$line")
    done    

    for line in `find -L . -maxdepth 1 -mindepth 1 -type f -printf "%f-%k\n" | sort -n`; do  # %f file name, %k file size

        #echo working on item `pwd`/$line item: ${line%%-*} size: ${line#*-}
        
        item=${line%%-*}
    
        # assume this is no .tar files
        #[[ "$item" == *tar ]] && cp "$item" "$2/$1" &&  continue
        
        size=${line#*-}

        [ "$size" -gt "1048576" ] && cp "$item" "$2" && continue # bigger than 1G, inore it

        files="$files\n$item"

        totalSize=$((totalSize + size))

        if [ "$totalSize" -gt "1048576" ]; then # bigger than 1G 
            (archiveFiles "$files" "$2")
            files=""
            totalSize="0"
        fi
    done
    if [ ! -z "$files" ]; then 
        (archiveFiles "$files" "$2")
    fi 
} 

sFolder="$2"

dFolder="$3"

core=$1  

[ -d "$sFolder" ] || { echo "Usage: $0 <cores> <sourceFolder> [destinationFolder]"; exit 1; }

[ -z "$dFolder" ] && dFolder="$sFolder-tar"

[ -d "$dFolder" ] || mkdir -p $dFolder 

dFolder=`realpath $dFolder`

dFolderTmp=`mktemp -d`

#rm -f $dFolderTmp/lock.* 2>/dev/null

startTime=`date` 

archiveFolder "$sFolder" "$dFolder"

wait

endTime=`date`

echo "Time used: $(($(date -d "$endTime" '+%s') - $(date -d "$startTime" '+%s')))" >> $dFolder/archive.log

