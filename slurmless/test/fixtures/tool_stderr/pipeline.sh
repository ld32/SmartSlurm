#!/bin/sh
number=$1
#@1,0,noisyTool,,,sbatch -p short -c 1 --mem 2G -t 50:0
noisyTool.sh $number

