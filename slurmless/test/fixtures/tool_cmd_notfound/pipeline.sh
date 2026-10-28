#!/bin/sh
number=$1
#@1,0,badcmd,,,sbatch -p short -c 1 --mem 2G -t 50:0
totally_not_a_real_command --do-stuff $number

