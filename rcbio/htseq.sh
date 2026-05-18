#!/bin/sh

#set -x

usage() { echo -e "Usage :\n${0##*/} <-p path_of_alignment, such as: bowtieOut, bwaOut or starOut, required> [-r species_index (required if no -g. Such as: dm3, dm6, mm10, hg18, hg19 or hg38. Let us know if you need other references)] [-g .gtf file with full path (required if no -r, don't need this if -r is given)] <-s strand(yes,no,reverse)>"; exit 1;}

while getopts ":p:r:s:g:" o; do
    case "${o}" in
        r)
            r=${OPTARG}
            ;;
        g)
            g=${OPTARG}
            ;;
        p)  
            p=${OPTARG}
            ;;
        s)
            s=${OPTARG}
            ;; 
        *)
            usage
            ;;
    esac
done

[ -d group2 ] || { echo group2 is not found. You need at least two groups to run this pipeline; exit 1; }

# rename flag folder from earlier bowtie/star/bwa alignment run 
#[ -f flag/alljobs.jid.first ] && { grep 2.1.mergeHTSeq flag/alljobs.jid.first >/dev/null || { mv -f flag flagEarlierRun; mkdir flag; } ; } 

[ -z "$p" ] && { echo Please provide the alignment folder for option -p such as -p bwaOut, -p bowtieOut or -p starOut;  usage; }

[ -d $p ] || { echo Folder not exist: $p; usage; }

[ -z "$s" ] && { echo Please provide the strand for option -s;  usage; }

#module load gcc/6.2.0  python/2.7.12  htseq/0.9.1 samtools/0.1.19
module load conda/miniforge3
conda activate /n/shared_db/misc/rcbio/htseqEnv 

echo Current loaded modules: `module list`

# set up gtf file paths 
path=`which sbatchRun`
source ${path%\/bin\/sbatchRun}/config/config.txt

if [ -z "${r}" ]; then
    if [ ! -z "$g" ]; then 
        gtf=$g
    else
        usage
    fi
    
else 
  case "$r" in
     "dm3")
       gtf="$gtfdm3"
    ;;
    
    "dm6")
       gtf="$gtfdm6"
    ;;
   
    "mm10")
       gtf="$gtfmm10"
    ;;
        
    "hg18")
       gtf="$gtfhg18"
    ;;
    
    "hg19") 
       gtf="$gtfhg19"
    ;;
    
    "hg38")
       gtf="$gtfhg38"
    ;;
     
    *)  echo "Index '$r' is not supported. Please email rchelp@hms.harvard.edu for help."; exit
    ;;
    
  esac
fi

[ -f $gtf ] || { echo .gtf file not exist; exit 1; }

pwdhere=`pwd`

for group in `ls -v -d group*/|sed 's|[/]||g'`; do
    echo working on group:  $group
    cd $pwdhere/$group
    for sample in `ls -d */|sed 's|[/]||g'`; do
        echo working on sample: $sample
        i=$pwdhere/$p/$group$sample/accepted_hits.bam
        [ -f $i ] || { echo bam file not find: $i; exit 1; }
        
        #@1,0,htseq,,sbatch -p short -t 12:0:0 --mem 8G
        samtools view -bf 1 $i > $i.pe.bam && samtools sort -n $i.pe.bam -o $i.sorted.bam && htseq-count -s $s -f bam $i.sorted.bam $gtf > $i.read.count.txt


    done
done

module load gcc/14.2.0 R/4.4.2

cd $pwdhere

#@2,1,mergeHTSeq,,sbatch -p short -t 2:0:0 --mem 4G
mergeHTSeqCount.r $p


