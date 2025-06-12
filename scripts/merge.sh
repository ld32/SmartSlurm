#!/bin/sh

set -x 

# deactivate base env
source deactivate

#export PATH=/n/data1/cores/bcbio/eclip/fastq-tools-0.8.3/bin:/n/data1/cores/bcbio/eclip/eCLIP/bin:$PATH
module load gcc/6.2.0 perl/5.24.0 miniconda3/23.1.0 samtools/1.9 R/3.3.3

conda activate /n/data1/cores/bcbio/eclip/eclipEnvPython3.9

sed "s#eclipOut#$PWD/eclipOut#g" /n/data1/cores/bcbio/eclip/merge_peaks/wf/merge_peaks_2inputs.yaml > merge_peaks_2inputs.yaml

chmod +x merge_peaks_2inputs.yaml

#export PATH=/opt/merge_peaks/wf:/opt/merge_peaks/bin:/opt/merge_peaks/bin/perl:$PATH

export PATH=/n/data1/cores/bcbio/eclip/merge_peaks/wf:/n/data1/cores/bcbio/eclip/merge_peaks/bin:/n/data1/cores/bcbio/eclip/merge_peaks/bin/perl:$PATH

#@1,0,merge,,,sbatch -c 20 -p short -t 12:0:0 --mem 10G
env; rm -r merge_peaks_2inputs || true; ./merge_peaks_2inputs.yaml

#which eCLIP_full_IDR_pipeline_2inputs_scatter_singleNode
#rm -r merge_peaks_2inputs || echo; singularity exec --bind /n/data1 /n/app/singularity/containers/ld32/idr_2.0.2.sif "export PATH=/n/data1/cores/bcbio/eclip/merge_peaks/wf:/n/data1/cores/bcbio/eclip/merge_peaks/bin:/n/data1/cores/bcbio/eclip/merge_peaks/bin/perl:$PATH; ./merge_peaks_2inputs.yaml"
