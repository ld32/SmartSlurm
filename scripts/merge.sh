#!/bin/sh

set -x 

# deactivate base env
source deactivate

#export PATH=/n/data1/cores/bcbio/eclip/fastq-tools-0.8.3/bin:/n/data1/cores/bcbio/eclip/eCLIP/bin:$PATH
module load gcc/6.2.0 perl/5.24.0 miniconda3/23.1.0 samtools/1.9 R/3.3.3

source activate /n/data1/cores/bcbio/eclip/eclipEnvPython3.9

sed "s#eclipOut#$PWD/eclipOut#g" /n/data1/cores/bcbio/eclip/merge_peaks/wf/merge_peaks_2inputs.yaml > merge_peaks_2inputs.yaml

chmod +x merge_peaks_2inputs.yaml

export PATH=/n/data1/cores/bcbio/eclip/merge_peaks/wf:/n/data1/cores/bcbio/eclip/merge_peaks/bin:/n/data1/cores/bcbio/eclip/merge_peaks/bin/perl:$PATH

#@1,0,merge,,,sbatch -c 12 -p short -t 12:0:0 --mem 40G
rm -r merge_peaks_2inputs || echo;  ./merge_peaks_2inputs.yaml
