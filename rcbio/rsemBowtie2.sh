#!/bin/sh

#set -x 
usage() { echo -e "Usage :\n${0##*/} [withTPMoutput, optional] <-r species_index (required, such as mm10, hg19 or GRCz10. Let us know if you need other references)>"; exit 1;} 

while getopts ":r:" o; do
    case "${o}" in
         r)
            r=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done

#module load gcc/6.2.0 rsem/1.3.0 bowtie2/2.5.4  samtools/1.3.1 igvtools/2.3.88 R/3.4.1 

module load gcc/14.2.0 bowtie2/2.5.4 samtools/1.21 R/4.4.2
module load conda/miniforge3
conda activate /n/shared_db/misc/rcbio/rcbioEnv

echo Current loaded modules: `module list`

# set up bowtie2 index paths 
path=`which sbatchRun`
source ${path%\/bin\/sbatchRun}/config/config.txt

case "$r" in
    "mm10")index="$rsemBowtie2mm10"
    ;;
  
    "hg19") index="$rsemBowtie2hg19" 
    ;;
    
    "GRCz10")index="$rsemBowtie2GRCz10"
    ;;    

    *)  echo "Index '$r' is not supported. Please email rchelp@hms.harvard.edu for help."; exit
    ;;
esac

#[ -d group2 ] || { echo group2 is not found. You need at least two groups to run this pipeline; exit 1; }

mkdir -p out

# go through sample groups
for group in `ls -v -d group*`; do
    echo working on group:  $group
    
    COUNTER=0 
    for sample in `ls -v -d $group/*  | xargs -n 1 basename`; do 
        echo working on sample: $sample
       
        ls  $group/$sample/*_1.fastq >/dev/null 2>&1 || ls $group/$sample/*_1.fq  >/dev/null 2>&1 || { echo Read file not found for $sample! Please make sure the fastq files are named as xxx_1.fastq, xxx_2.fastq or xxx_1.fq, xxx_2.fq;  exit 1; }
        
        reads1=""; reads2=""        
        for r1 in `ls $group/$sample/*_1.fastq $group/$sample/*_1.fq 2>/dev/null | xargs -n 1 basename`; do 
            readgroup=${r1##*/}
            readgroup=${readgroup%_*}
            r2=${r1%_*}_2${r1##*_1}
            echo working on readgroup: $readgroup
			reads1="$reads1$group/$sample/$r1," 
            
            [[ -f $group/$sample/$r2 ]] && reads2="$reads2$group/$sample/$r2," || { echo -e "\n!!!Warning: read2 file '$r2' not exist, ignore this warning if you are working with single-end data\n"; }
        done
        
        #echo reads1 $reads1
        #echo reads2 $reads2         
        
        [ -z $reads2 ] || reads1="--paired-end $reads1"
                 
        #@1,0,expression,index,sbatch -c 4 -p short -t 8:0:0 --mem 32G 
        rsem-calculate-expression -p 4 -q  --bowtie2 --estimate-rspd --output-genome-bam  ${reads1%,} ${reads2%,} $index out/$group.$sample && mkdir -p /tmp/$USER/$group/$sample && samtools sort -@ 4 -m 1G out/$group.$sample.transcript.bam -T /tmp/$USER/$group/$sample -o out/$group.$sample.transcript.sorted.bam && samtools index out/$group.$sample.transcript.sorted.bam && samtools sort -@ 4 -m 1G out/$group.$sample.genome.bam -T /tmp/$USER/$group/$sample -o out/$group.$sample.genome.sorted.bam && samtools index out/$group.$sample.genome.sorted.bam
 
        
        # bowtie1 
        #rsem-calculate-expression -p 4 -q  --bowtie-chunkmbs 200 --output-genome-bam ${reads1%,} ${reads2%,} $index out/$group.$sample 
               
        #@2,1,wigTranscript
        rsem-bam2wig out/$group.$sample.transcript.sorted.bam out/$group.$sample.transcript.wig out/$group.$sample.transcript && igvtools toTDF out/$group.$sample.transcript.wig out/$group.$sample.transcript.tdf $index.transcripts.genome
        
        #@3,1,wigGenome
        rsem-bam2wig out/$group.$sample.genome.sorted.bam out/$group.$sample.genome.wig out/$group.$sample.genome 
        
        #@4,1,plot
        rsem-plot-model out/$group.$sample out/$group.$sample.plot.pdf
        
        COUNTER=$((COUNTER + 1))          
        generesults="$generesults out/$group.$sample.genes.results"
        isoresults="$isoresults out/$group.$sample.isoforms.results"     
    done 
    [[ "$COUNTER" -eq 0 ]] || counts="$counts,$COUNTER"
done 

#@5,1,iso
rsem-generate-data-matrix $isoresults > iso.count.matrix && { [[ ! "$1" == withTPMoutput ]] || rsem-generate-data-matrix-tpm  out/*isoforms.results > iso.tpm.matrix.txt; } && rsem-run-ebseq --ngvector $index.ngvec iso.count.matrix ${counts#,} all.iso.results && rsem-control-fdr all.iso.results 0.05 final.iso.results.fdr.txt

#@6,1,gene
rsem-generate-data-matrix $generesults > gene.count.matrix && { [[ ! "$1" == withTPMoutput ]] ||rsem-generate-data-matrix-tpm  out/*genes.results > gene.tpm.matrix.txt; } && rsem-run-ebseq gene.count.matrix ${counts#,} all.gene.results && rsem-control-fdr all.gene.results 0.05 final.gene.results.fdr.txt
