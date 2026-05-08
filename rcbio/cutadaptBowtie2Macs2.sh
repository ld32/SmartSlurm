#!/bin/sh

#set -x
usage() { echo -e "Usage :\n${0##*/} [-r species_index (required if no -b. Such as: dm3, dm6, mm10, hg18, hg19 or hg38. Let us know if you need other references)] [-b bowtie2IndexWithPath(required if no -r, don't need this if -r is given)]"; exit 1;}

while getopts ":r:b:" o; do
    case "${o}" in
        r)
            reference=${OPTARG}
            ;;
        b)
            bowtieIndex=${OPTARG}
            ;;
    esac
done

#module load gcc/6.2.0 bowtie2/2.2.9 samtools/0.1.19
SHARED_DATABASES=/n/shared_db/igenome/03032016/

if [ -z "${reference}" ]; then
    if [ ! -z "$bowtieIndex" ]; then
        bowtie2-inspect -n ${bowtieIndex} &> /dev/null || { echo -e "Error: \ngenome bowtie2 index could not be found: $bowtieIndex"; usage; }
        index=$bowtieIndex
    else
        usage
    fi

else
  case "$reference" in
    "mm10")index="$SHARED_DATABASES/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome"
    ;;

    "dm3") index="$SHARED_DATABASES/Drosophila_melanogaster/UCSC/dm3/Sequence/Bowtie2Index/genome"
    ;;

    "dm6") index="$SHARED_DATABASES/Drosophila_melanogaster/UCSC/dm6/Sequence/Bowtie2Index/genome"
    ;;

    "hg18") index="$SHARED_DATABASES/Homo_sapiens/UCSC/hg18/Sequence/Bowtie2Index/genome"
    ;;

    "hg19") index="$SHARED_DATABASES/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome"
    ;;

    "hg38") index="$SHARED_DATABASES/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome"
    ;;

    *)  echo "Index '$r' is not supported. Please email rchelp@hms.harvard.edu for help."; exit
    ;;

  esac
fi

#module load gcc/6.2.0 samtools/1.3.1 bowtie2/2.2.9 fastx/0.0.13 python/2.7.12  macs2/2.1.1.20160309 deeptools/3.0.2


module load conda/miniforge3

conda activate /n/shared_db/misc/rcbio/macsEnv

mkdir -p fq bam macs2 

rm fq/* 2>

for lib in lib/*; do
     #@1,0,split,,sbatch -n 1 -p short -t 2:0:0 --mem 8G
    cutadapt -g ^file:barcodes.fa -o fq/{name}.fq $lib
    break
done

for bqName  in `grep "^>" barcodes.fa | grep "Co$" | sed 's/^>//'`; do
    
    #@2,1,bowtie2,index,sbatch -n 4 -p short -t 12:0:0 --mem 40G
    bowtie2 -p 4 -q --local -x $index -U fq/$bqName.fq -S bam/$bqName.sam 
    
    #@3,2,samtools,,sbatch -n 1 -p short -t 2:0:0 --mem 8G
    samtools view -h -S bam/$bqName.sam -b -o bam/$bqName.bam && samtools sort bam/$bqName.bam > bam/$bqName.sorted.bam && samtools index bam/$bqName.sorted.bam 
    
    #@4,3,bamCoverage,,sbatch -n 4 -p short -t 12:0:0 --mem 40G
    bamCoverage -b bam/$bqName.sorted.bam -o bam/$bqName.sorted.bam.bw --binSize 5 --normalizeUsing RPKM --smoothLength 30 --extendReads 200 -p 4
    
    bqName1=${bqName%Co}Tr
    
    #@5,1,bowtie21,index,sbatch -n 4 -p short -t 12:0:0 --mem 40G
    bowtie2 -p 4 -q --local -x $index -U fq/$bqName1.fq -S bam/$bqName1.sam 
    
    #@6,5,samtools1,,sbatch -n 1 -p short -t 2:0:0 --mem 8G
    samtools view -h -S bam/$bqName1.sam -b -o bam/$bqName1.bam && samtools sort bam/$bqName1.bam > bam/$bqName1.sorted.bam && samtools index bam/$bqName1.sorted.bam  
    
    #@7,6,bowtie21,index,sbatch -n 4 -p short -t 12:0:0 --mem 40G
    bamCoverage -b bam/$bqName1.sorted.bam -o bam/$bqName1.sorted.bam.bw --binSize 5 --normalizeUsing RPKM --smoothLength 30 --extendReads 200 -p 4 
    
    #@8,3.6,macs2,,sbatch -n 4 -p short -t 2:0:0 --mem 40G
    macs2 callpeak -t bam/$bqName1.sorted.bam -c bam/$bqName.sorted.bam -f BAM -g hs --nomodel --broad -p 1e-9 --broad-cutoff 0.05 -n K27 --outdir macs2/${bqName%Co} -p 4
    #break
done 
