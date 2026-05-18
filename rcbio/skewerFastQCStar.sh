#!/bin/sh

usage() { echo -e "\nUsage : ${0##*/} <-l readLength, required> -a <adapter.fa (required)>  [-r speciesIndex (required if no -b. Such as: mm10,hg18 or hg19. Let us know if you need other references)] [-b starIndexWithPath(required if no -r, don't need this if -r is given)] [-f genomeFastaFile(required if no -r, don't need this if -r is given)] [-g gtfFileWithPath (optional, don't need this if -r is given)] "; exit 1;} 

while getopts ":r:b:f:g:l:a:" o; do
    case "${o}" in
        r)
            r=${OPTARG}
            ;;
        b)
            b=${OPTARG}
            ;;
        f)
            f=${OPTARG}
            ;;
        g)
            g=${OPTARG}
            ;;
        l)
            l=${OPTARG}
            ;;
        a)
            a=${OPTARG}
            ;; 
       *)
            usage
            ;;
    esac
done

[ -f $a ] || { echo Adapter sequence file is missing or not exist $a; usage; }

[ -z "$l" ] && echo Please give read length && usage 

#module load  gcc/6.2.0 star/2.5.4a  samtools/1.3.1 picard/2.8.0  skewer/0.2.2 fastqc/0.11.5

module load conda/miniforge3
conda activate /n/shared_db/misc/rcbio/rcbioEnv
module load  gatk/4.6.1.0

echo Current loaded modules: `module list`

# set up paths 
path=`which sbatchRun`
source ${path%\/bin\/sbatchRun}/config/config.txt

# set the correct index and gene annotation file
if [ -z "${r}" ]; then
    if [ ! -z "$b" ]; then 
        [ -f $b/genomeParameters.txt ] || { echo -e "Error: \ngenome star index could not be found: $b"; usage; } 
        index=$b
    else
        usage
    fi
    [ -f ${f} ] || { echo -e "Error: \n1genome fasta file could not be found: ${f}"; usage;}  
    
    if [ ! -z "${g}" ]; then
       [ -f "${g}" ] || { echo -e "Error: \ngtf file could not be found: $g"; usage; } 
       gtf="-G $g"
    fi
    
else 
  case "$r" in
    "mm10")index="starmm10index"
       gtf="-G $starmm10gtf"
    ;;

    "hg18") index="$starhg18index"
        gtf="-G $starhg18gtf"
    ;;

    "hg19") index="$starhg19index"
        gtf="-G $starhg19gtf"
    ;;
    
    *)  echo "Index '$r' is not supported. Please email rchelp@hms.harvard.edu for help."; exit
    ;;
  esac
  fa=$index/genome.fa
  
  [ -f $fa ] || { echo genome fasta file not exist: $fa; usage; } 
fi

[ -d group2 ] || { echo group2 is not found. You need at least two groups to run this pipeline; exit 1; }

pwdhere=`pwd`

#loopStart,group
for group in `ls -v -d group*`; do
    echo working on group:  $group
   
    #loopStart,sample
    for sample in `ls -d $group/*  | xargs -n 1 basename`; do 
        echo working on sample: $sample
        inputsams=""
        ls  $group/$sample/*_1.fastq* >/dev/null 2>&1 || ls $group/$sample/*_1.fq*  >/dev/null 2>&1 || { echo Read file not found for $sample! Please make sure the fastq files are named as xxx_1.fastq, xxx_2.fastq or xxx_1.fq, xxx_2.fq;  exit 1; }
          
        #loopStart,readgroup
        for r1 in `ls $group/$sample/*_1.fastq* $group/$sample/*_1.fq* 2>/dev/null | xargs -n 1 basename`; do 
            readgroup=${r1##*/}
            readgroup=${readgroup%_*}
            r2=${r1%_*}_2${r1##*_1}
            echo working on readgroup: $readgroup
            r1=$group/$sample/$r1 
            
            mode=pe 
            [[ -f $group/$sample/$r2 ]] && r2=$group/$sample/$r2 || { echo -e "\n\n!!!Warning: read2 file '$r2' not exist, ignore this warning if you are working with single-end data\n\n"; r2=""; mode=any;  }
            
            out=skewer/$group.$sample.$readgroup
            
            unzipcmd=""
            [[ $r1 == *gz ]] && unzipcmd="zcat $r1  > $out/unzip/r1.fq; r1=$out/unzip/r1.fq;"  
            [[ $r1 == *bz2 ]] && unzipcmd="bzcat $r1  > $out/unzip/r1.fq; r2=$out/unzip/r2.fq;"
            
            unzipcmd1=""
            [[ $r2 == *gz ]] && unzipcmd1="zcat $r2  > $out/unzip/r2.fq; r1=$out/unzip/r1.fq;"  
            [[ $r2 == *bz2 ]] && unzipcmd1="bzcat $r2  > $out/unzip/r2.fq; r2=$out/unzip/r2.fq;"
           
            #@1,0,skewer,,sbatch -p short -t 2:0:0 -c 4 --mem 20G
            mkdir -p $out/unzip; $unzipcmd $unzipcmd1 skewer -t 4 -x $a -m $mode $r1 $r2 -o $out/SKEWER
            
            out=fastqc/$group.$sample.$readgroup; 
            
            #@2,0,fastqc,,sbatch -p short -t 2:0:0 -c 1 --mem 8G
            mkdir -p $out; fastqc -o $out $r1 $r2  
            
            [ -z $r2 ] && r1="$pwdhere/skewer/$group.$sample.$readgroup/SKEWER-trimmed.fastq" || { r1="$pwdhere/skewer/$group.$sample.$readgroup/SKEWER-trimmed-pair1.fastq";  r2="$pwdhere/skewer/$group.$sample.$readgroup/SKEWER-trimmed-pair2.fastq";}
            
            zipcmd=""
            #[[ $r1 == *gz ]] && zipcmd="--readFilesCommand zcat"  
            #[[ $r1 == *bz2 ]] && zipcmd="--readFilesCommand bzcat"    
            
            mkdir -p starOut/$group$sample$readgroup-star.p1
            mkdir -p starOut/$group$sample$readgroup-star.ref
            mkdir -p starOut/$group$sample$readgroup-star.p2
            
            # first pass 
            #@3,1,star,index.gtf,sbatch -p short -c 4 -t 0-8:0 --mem 40G 
            rm -r $pwdhere/starOut/$group$sample$readgroup-star.p1/* $pwdhere/starOut/$group$sample$readgroup-star.ref/* $pwdhere/starOut/$group$sample$readgroup-star.p2/*; cd $pwdhere/starOut/$group$sample$readgroup-star.p1 && STAR --genomeDir $index --readFilesIn $r1 $r2 --runThreadN 4 $zipcmd  && cd $pwdhere/starOut/$group$sample$readgroup-star.ref && STAR --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles $fa --sjdbFileChrStartEnd $pwdhere/starOut/$group$sample$readgroup-star.p1/SJ.out.tab --sjdbOverhang $(($l - 1)) --runThreadN 4  &&  cd $pwdhere/starOut/$group$sample$readgroup-star.p2 && STAR --genomeDir $pwdhere/starOut/$group$sample$readgroup-star.ref --readFilesIn $r1 $r2 --runThreadN 4  $zipcmd --outSAMstrandField intronMotif    
            
            # # the awk command is to correct a issue: https://groups.google.com/forum/#!topic/rna-star/Ta1Z2u4bPfc
            # second pass
            ##@3,2,star3:star3 command 
            # && cd starOut/$group$sample$readgroup-star.p2 && STAR --genomeDir ../../starOut/$group$sample$readgroup-star.ref --readFilesIn $r1 $r2 --runThreadN $n  $zipcmd --outSAMstrandField intronMotif --genomeLoad LoadAndRemove && runAwk && cd $pwdhere          
                   
            # add readgroup
            ##@4,3,star4:star4 command 
            #java -Xmx4g -Djava.io.tmpdir=/tmp -jar $PICARD/picard.jar AddOrReplaceReadGroups I=$sample/star.p2.$readgroup/Aligned.out.sam O=$sample/star.p2.$readgroup/out.bam SO=coordinate  RGID=$sample.$readgroup RGLB=rglb  RGPU=rgpu RGPL=illumina RGSM=$sample
            
            inputsams="$inputsams INPUT=starOut/$group$sample$readgroup-star.p2/Aligned.out.sam"
            
        #loopEnd 	
        done
        
        
        cd $pwdhere           
        mkdir starOut/$group$sample
        
        #@4,3,merge,,sbatch -p short -c 1 -t 0-4:0 --mem 10G  
        gatk MergeSamFiles VALIDATION_STRINGENCY=SILENT OUTPUT=starOut/$group$sample/accepted_hits.bam  $inputsams SORT_ORDER=coordinate && samtools index starOut/$group$sample/accepted_hits.bam
        
        #break
    done 
    #break
done


