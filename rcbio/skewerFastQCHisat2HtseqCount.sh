#!/bin/sh

usage() { echo -e "\nUsage : $0 <-r species_index -a <adapter.fa> (required. Such as: dm6, GRCz10, GRCz11, mm10 or hg19. Let us know if you need other references)> <-s strand(yes,no,reverse)>"; exit 1;} 

while getopts ":a:r:s:" o; do
    case "${o}" in
        r)
            r=${OPTARG}
            ;;
        a)
            a=${OPTARG}
            ;;    
        s)  
            s=${OPTARG} 
            ;;   
            
        *)
            usage
            ;;
    esac
done

[ -z $r ] && { echo Genome reference is needed; usage; }
[ -f $a ] || { echo Adapter sequence file is missing or does not exist $a; usage; }
[ -z "$s" ] && { echo Please provide the strand for option -s;  usage; }

#module purge 
#module load gcc/6.2.0 skewer/0.2.2 fastqc/0.11.5 hisat2/2.1.0 samtools/1.3.1   python/2.7.12 htseq/0.9.1 R/3.4.1

module load conda/miniforge3
conda activate /n/shared_db/misc/rcbio/rcbioEnv

echo Current loaded modules: `module list`

# set up paths 
path=`which sbatchRun`
source ${path%\/bin\/sbatchRun}/config/config.txt

case "$r" in
    "dm6" ) index="$hisat2dm6Index"
            splice="--known-splicesite-infile $hisat2dm6Splice"
            gtf="$hisat2dm6GTF"
    ;;	
    
    "mm10") index="$hisat2mm10Index"
            splice="--known-splicesite-infile $hisat2mm10Splice"
            gtf="$hisat2mm10GTF"
    ;;
    "hg19") index="$hisat2hg19Index"
            splice="--known-splicesite-infile $hisat2hg19Splice"
            gtf="$hisat2hg19GTF"
    ;;
    
    "GRCz10") index="$hisat2GRCz10Index"
            splice="--known-splicesite-infile $hisat2GRCz10Splice"
            gtf="$hisat2GRCz10GTF"
    ;;
   
    "GRCz11") index="$hisat2GRCz11Index"
            splice="--known-splicesite-infile $hisat2GRCz11Splice"
            gtf="$hisat2GRCz11GTF"
    ;; 
    
   *)  echo "Index '$r' is not supported. Please email rchelp@hms.harvard.edu for help."; usage
    ;;
esac
 
#[ -d group2 ] || { echo group2 is not found. You need at least two groups to run this pipeline; exit 1; }

pwd=`pwd`

for group in `ls -v -d group*/|sed 's|[/]||g'`; do
    echo working on group:  $group
    #loopStart,sample
    for sample in `ls -d $group/*/ | xargs -n 1 basename`; do
         
        echo working on sample: $sample, read1 files: 
        ls  $group/$sample/*_1.fastq* 2>/dev/null || ls $group/$sample/*_1.fq*  2>/dev/null || { echo Read file not found for $sample! Please make sure the fastq files are named as xxx_1.fastq, xxx_2.fastq or xxx_1.fq, xxx_2.fq;  exit 1; }
        read1=""
        read2=""
        
        for r1 in `ls $group/$sample/*_1.fastq* $group/$sample/*_1.fq* 2>/dev/null | xargs -n 1 basename`; do 
            #echo r1 is $r1
            readgroup=${r1%_*}
            echo working on readgroup: $readgroup
            r2=${r1%_*}_2${r1##*_1}
            #echo R2 is $r2
            
			[[ -f $group/$sample/$r2 ]] && r2="$group/$sample/$r2" || { echo -e "\n\n!!!Warning: read2 file '$r2' not exist, ignore this warning if you are working with single-end data\n\n"; r2=""; }
            out=skewer/$group.$sample.$readgroup
            r1=$group/$sample/$r1
           
            #@1,0,skewer,,sbatch -p short -t 2:0:0 -c 4 --mem 20G
            mkdir -p $out; skewer -t 4 -x $a -m pe $r1 $r2 -o $out/SKEWER
            
            read1="$read1,$pwd/$out/SKEWER-trimmed-pair1.fastq"; [ -z $r2 ] || read2="$read2,$pwd/$out/SKEWER-trimmed-pair2.fastq"
            
            out=fastqc/$group.$sample.$readgroup; 
            
            #@2,0,fastqc,,sbatch -p short -t 2:0:0 -c 1 --mem 8G
            mkdir -p $out; fastqc -o $out $r1 $r2  
        done
        
        [ -z "$read2" ] && reads="-U ${read1#,}" || reads="-1 ${read1#,} -2 ${read2#,}"
        
        #echo reads: $reads  
        
        out=hisat/$group.$sample; mkdir -p $out htseq 
        
        #@3,1,hisatCount,,sbatch -c 4 -p short -t 12:0:0 --mem 40G
        cd $out; rm aligns.sorted.bam.tmp* 2>/dev/null; hisat2 $index $reads --phred33 --mm -p 4 --dta  $splice | samtools view -Suh -f 1 - |  samtools sort - -n -o aligns.sorted.bam && htseq-count -s $s -f bam aligns.sorted.bam $gtf > $pwd/htseq/$group.$sample.read.count.txt; cd -
                 

        #cd $out; hisat2 $index $reads --phred33 --mm -p 4 --dta  $splice 2>hisat2.log | samtools view -Suh -f 1 - |  samtools sort - -n -o aligns.sorted.bam; htseq-count -s no -f bam aligns.sorted.bam $gtf > $pwd/htseq/$group.$sample.read.count.txt; cd -            
        
        # old bowtie align command       
        #mkdir -p bowtieOut/$group$sample;  bowtie2  $p  -p $n -x $index -1 ${read1#,} -2 ${read2#,} | samtools view -bS - > bowtieOut/$group$sample/accepted_hits.bam
        
        # old htseq command 
        #samtools view -bf 1 $i > $i.pe.bam && samtools sort -n $i.pe.bam $i.sorted && htseq-count -s $s -f bam $i.sorted.bam $gtf > $i.read.count.txt"
        
        #exit  
    done 
done

#@4,3,mergeCount
merge.count.r htseq
