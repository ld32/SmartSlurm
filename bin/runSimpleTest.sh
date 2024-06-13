#!/bin/bash

set -x 

usage() { echo -e "Usage: \n${0##*/} Rename jobRecord folder and job log folder, and start new run using bashScriptV2.sh"; exit 1; }

[ -z "$1" ] || usage

if [ -z "$smartSlurmJobRecordDir" ]; then
    if [ -f ~/.smartSlurm/config/config.txt ]; then
        source ~/.smartSlurm/config/config.txt
    else
        source $(dirname $0)/../config/config.txt || { echo Config list file not found: config.txt; exit 1; }
    fi
fi

if [ -d $smartSlurmJobRecordDir ]; then 
    echo Do you want to rename job record folder? 
    echo y: Click y to rename
    echo Enter: Skip and do not rename
    read -p "" x </dev/tty
    [[ "$x" == "y" ]] && mv $smartSlurmJobRecordDir $smartSlurmJobRecordDir.$(stat -c '%.19z' $smartSlurmJobRecordDir | tr " " "." | tr ":" "-")
fi 

[ -d $smartSlurmLogDir ] && mv -r $smartSlurmLogDir/* $smartSlurmLogDir.$(stat -c '%.19z' $smartSlurmLogDir | tr " " "." | tr ":" "-")

[ -d $smartSlurmLogDir ] && rm -fr $smartSlurmLogDir/* 

runAsPipeline "$(dirname $0)/../scripts/bashScriptV2x.sh 1234" "sbatch -A rccg -p short -c 1 --mem 2G -t 50:0" noTmp run
