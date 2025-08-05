#!/bin/bash

#set -x 
#set -e

function unArchiveFolder() {
    cd $1
    rm *.list.txt 2>/dev/null
    for item in `find . -maxdepth 1 -mindepth 1`; do 
        [ -d "$item" ] && (unArchiveFolder "$item") && continue
        if ([[ "$item" == *tar ]]); then 
            while true; do
                jobID=0
                for i in $(seq 1 $core); do
                    if `mkdir $dFolder/lock.$i 2>/dev/null`; then
                        jobID=$i
                        break
                    fi
                done
                [ "$jobID" -eq "0" ] && sleep 10 || break
            done 
            md5sum -c ${item%.tar}.md5sum &&
            tar xf "$item" && rm "$item" && rm -r $dFolder/lock.$jobID && echo xvf "$item" | tee -a  $dFolder/unArchive.log || 
            echo Failed checksum for $item >> $dFolder/unArchive.log &
        fi
    done    
}

sFolder="$2"
dFolder="$3"

[ -d "$sFolder" ] || { echo "Usage: $0 <cores> <sourceFolder> <tarFolder>"; exit 1;  }
[ -d "$dFolder" ] || { echo "Usage: $0 <cores> <sourceFolder> <tarFolder>"; exit 1;  }

dFolder=`realpath $dFolder`

cwd=`pwd`
core=$1
echo untar start:  >> $dFolder/unArchive.log
startTime=`date`
unArchiveFolder "$dFolder"
echo diff -r $sFolder $dFolder | tee -a $dFolder/unArchive.log

cd $cwd
diff -r $sFolder $dFolder | tee -a $dFolder/unArchive.log
endTime=`date`
echo "Time used: $(($(date -d "$endTime" '+%s') - $(date -d "$startTime" '+%s')))" | tee -a  $dFolder/unArchive.log
