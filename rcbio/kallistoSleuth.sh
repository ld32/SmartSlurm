#!/bin/sh
usage() { echo -e "\nUsage : `basename $0` <-i transcriptome index (required: currently we have: GRCh38, GRCm38, BDGP6 and GRCz10)> [ -l fragmentLength, optional. Only need this for single end reads] [-s Estimated standard deviation of fragment length, optional. Only need this for single send reads]\n"; exit 1;} 

while getopts ":i:l:s:" o; do
    case "${o}" in
        i)
            i=${OPTARG}
            ;;

        l)  l=${OPTARG}
	        ;;
	        
	    s)  s=${OPTARG}
	        ;;

        *)
            usage
            ;;

    esac
done

[ -z "$i" ] && usage

[ ! -z "$l" ] && [ -z "$s" ] && usage

[ -z "$l" ] && [ ! -z "$s" ] && usage

#module load gcc/6.2.0 kallisto/0.43.1  R/3.4.1 
module load  gcc/14.2.0 R/4.4.2 
module load conda/miniforge3
conda activate /n/shared_db/misc/rcbio/rcbioEnv

echo Current loaded modules: `module list`

# set up kallisto index paths 
path=`which sbatchRun`
source ${path%\/bin\/sbatchRun}/config/config.txt


SHARED_DATABASES=/n/shared_db

# set the correct index and gene annotation file
if [ -f "${i}" ]; then
   index=$i
else 
  case "$i" in
    "GRCh38")index="$kallisto0_51_1GRCh38"
            sp=hsapiens_gene_ensembl
    ;;
    
    "GRCm38")index="$kallisto0_51_1GRCm38"
            sp=mmusculus_gene_ensembl
    ;;
    
    "BDGP6")index="$kallisto0_51_1BDGP6"
            sp=dmelanogaster_gene_ensembl
    ;;
    
    "GRCz10")index="$kallisto0_51_1GRCz10"
            sp=drerio_gene_ensembl
    ;;
  
  
    *)  echo "Index '$i' is not supported. Please email rchelp@hms.harvard.edu for help."; usage
    ;;
  esac
fi

[ -d group2/ ] || { echo group2 is not found. You need at least two groups to run this pipeline; exit 1; }

echo -e "sample\tgroup\tpath" > sample.lst 

mkdir -p kallistoOut

for group in `ls -v -d group*/|sed 's|[/]||g'`; do
  
    echo working on group:  $group
    for sample in `ls -d $group/*/ | xargs -n 1 basename`; do
         
        echo working on sample: $sample
        fileList=`ls  $group/$sample/*{_1.fastq,_1.fq,_r1.fastq,_r1.fq,_R1.fastq,_R1.fq}* 2>/dev/null`
        echo fileList: $fileList
        
        [ -z "$fileList" ] && { echo Read file not found for $group/$sample! Please make sure the fastq files are named as xxx_1.fastq, xxx_2.fastq or xxx_1.fq, xxx_2.fq or xxx_r1.fastq, xxx_r2.fastq or xxx_R1.fq, xxx_R2.fq or xxx_1.fastq.gz, xxx_2.fastq.gz or xxx_1.fq.gz, xxx_2.fq.gz or xxx_r1.fastq.gz, xxx_r2.fastq.gz or xxx_R1.fq.gz, xxx_R2.fq.gz;  exit 1; }
        
        reads=""
        for r1 in `echo $fileList | xargs -n 1 basename`; do 
            #echo working on read one file: $r1
            
            readgroup=${r1%_*}
            echo working on readgroup: $readgroup
            
            if [  -z "$l" ] ; then
                ext=${r1##*_}; ext=${ext/1/2}
                r2=${readgroup}_$ext
                reads="$reads  $group/$sample/$r1 $group/$sample/$r2"
            else                       
                reads="$reads  $group/$sample/$r1"     
            fi
        done
        
        echo -e "$sample\t$group\tkallistoOut/$group$sample" >> sample.lst 
        
        [ -z "$l" ] ||  reads="--single $reads -l $l -s $s"
        
        #@1,0,kallisto2,,sbatch -p short --mem 32G -c 4 -t 2:0:0 
        rm -r kallistoOut/$group$sample 2>/dev/null ; kallisto quant -o kallistoOut/$group$sample -b 100 -t 4 -i $index $reads
          
    done     
done

#@2,1,sleuth,,sbatch -p short --mem 10G -t 2:0:0
Rscript sleuth.r $sp
