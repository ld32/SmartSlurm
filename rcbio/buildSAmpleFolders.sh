#!/bin/bash 

[ -z "$1" ] && { echo "Usage: buildSampleFolders.sh <sample_sheet.xlsx>"; exit 1; }

[ -f "$1" ] || { echo "Usage: buildSampleFolders.sh <sample_sheet.xlsx>"; exit 1; }

module load conda/miniforge3

conda activate /n/shared_db/misc/rcbio/rcbioEnv

buildSampleFoldersFromSampleSheet.py $1 
