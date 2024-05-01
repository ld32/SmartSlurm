#!/bin/bash
# RNA-seq job Submission pipeline

#set -x
#set -e

if [[ "$1" == "-h" ]] | [[ "$1" == "--help" ]] | [[ -z "$1" ]]; then
    echo "Usage: $(basename "$0") -h or --help

    Job submission including manifest:
    $(basename "$0") manifest.txt reference celltype numbermismatches spike trimpos1 trimposN STARversion yes|no numeratorCond denominatorCond

    Job submission per sample:
    $(basename "$0") read1.fastq read2.fastq samplename reference celltype numbermismatches spike trimpos1 trimposN STARversion yes|no"

    exit
fi

export PATH=/n/data1/cores/ntc/scripts/AdelmanLab/RNAseqMappingSmartSlurm:$PATH

settngsandpaths=`which settings_and_paths.sh`
echo settings file $settngsandpaths
scriptsPath=$(dirname $settngsandpaths)
echo script path pulled from settings file $scriptsPath

gtfPath=`awk '{if($1=="gtfdir"){print $2}}' $settngsandpaths`
echo "gtf path pulled from settings file ${gtfPath}"
spikeinbase=`awk '{if($1=="spikein"){print $2}}' $settngsandpaths`
indexdir=`awk '{if($1=="indexdir"){print $2}}' $settngsandpaths`
rRNA=`awk '{if($1=="rRNA"){print $2}}' ${settngsandpaths}`

module purge
module load gcc/6.2.0 python/2.7.12 cutadapt/1.14 bowtie/1.2.2 samtools/1.9 bedtools/2.27.1 R/3.6.1 ucsc-tools/363 fastqc/0.11.9

output=STAR_mapping_stats.txt
if [ ! -f $output ]; then
    header="Samplename""\t""Number Total raw reads""\t""Number Spike Reads uniquely mapped to genome""\t""% Spike Reads uniquely mapped""\t""Number rRNA Ref Reads uniquely mapped to genome""\t""% rRNA Ref Reads uniquely mapped to genome""\t""Number Ref Reads uniquely mapped to genome""\t""% Ref Reads uniquely mapped to genome""\t""Number Reads Failed to Align against Ref genome""\t""% Perc Reads Failed Align to Ref genome""\t""Number Reads Multiply Aligned to Ref genome""\t""% Perc Reads multiply mapped to Ref genome""\t""Number Duplicated Reads mapped to Ref genome""\t""% Duplicated Reads mapped to Ref genome""\t""Total Number of Splice Junctions""\t""% Reads with Splice Junctions""\t""# Splice Junctions with uniq mapped reads""\t""Number of Mito Reads""\t""% of Mito Reads"
    echo -e $header >> $output
fi
 

if [[ "$1" != *".txt" ]]; then
    echo -e "$1 $2 $3 empty\n" > tmpSampleSheet.txt
    set -- "${@::1}" 'tmpSampleSheet.txt' "${@:4}"
else
    echo "manifest file detected. run dos2unix on manifest"
    dos2unix $1
    endline=`tail -c 1 "$1"`
    if [ -z "$(tail -c 1 "$1")" ]; then
        echo "Newline already at end of file"
    else
        echo "Adding newline at end of file"
        echo -e "\n" >> $1
    fi
fi

mkdir -p logs

echo "Input parameters to the pipeline:
manifest = $1
Genome ref = $2
cell type = $3
number mismatches = $4
genome Spike-In = $5
min trim = $6
max trim = $7
STAR version = ${8}
Cutadapt used? = ${9}
Numerator conditions = ${10}
Denominator conditions = ${11} " > logs/pipeline_parameters.log

cat logs/pipeline_parameters.log

manifest=$1

if [ ! -f $manifest ]; then
    echo "Usage: $(basename "$0") -h or --help

    Job submission including manifest:
    $(basename "$0") manifest.txt reference celltype numbermismatches spike trimpos1 trimposN STARversion yes|no numeratorCond denominatorCond
    
    Job submission per sample:
    $(basename "$0") read1.fastq read2.fastq samplename reference celltype numbermismatches spike trimpos1 trimposN STARversion yes|no"

    exit
fi
genomeRef=$2
celltype=$3
numbermismatches=$4
genomeSpikeIn=$5
mintrim=$6
maxtrim=$7
STARversion=$8
useCutadapt=${9}
numerator_cond=${10}
denominator_cond=${11}
gtf=${gtfPath}/${genomeRef}/${genomeRef}.basic.gtf

#Check what STAR version is being used and load that module
if [ -z "$STARversion" ]; then
    STARversion="STAR273a"
    echo "no STAR version specified, setting to version to 2.7.3a"
    module load star/2.7.3a
elif [ ${STARversion} == "STAR273a" ]
then
    echo "STAR version 2.7.3a detected"
    module load star/2.7.3a
elif [ ${STARversion} == "STAR270f" ]
then
    echo "STAR version 2.7.0f detected"
    module load star/2.7.0f
else
    echo "You must enter either STAR273a or STAR270f for STAR version for this to work"
    exit
fi

{
while read -r f1 f2 f3 f4
do
    [[ -n "$f3" ]] || continue #only submit jobs with sample name

    ls ${f1/,/ } || { echo Read1 file ${f1/,/ } does not exist; exit 1; }
    ls ${f2/,/ } || { echo Read2 file ${f2/,/ } does not exist; exit 1; }

    #Check if one or multiple fastqs submitted (per R1 or R2)
    numFastqR1=`echo $f1 | awk 'BEGIN {FS=","} {print NF}'`
    numFastqR2=`echo $f2 | awk 'BEGIN {FS=","} {print NF}'`

    #if # of fastqs per read group are uneven, give an error message and terminate the pipeline
    if [ $numFastqR1 -ne $numFastqR2 ]; then
        echo "unequal numbers of fastqs, this is likely a serious problem so terminating the pipeline"
        exit
    elif [ "$numFastqR1" -eq 0 ]; then
        echo "No fastqs detected - terminating the pipeline"
        exit
    fi

done
} < "$manifest"

wkdir=`pwd`; mkdir -p fastqc fastq bowtierRNA raw_coverage_tracks SpikeIn mapping DEseqOutput logs
echo "bam" > DEseq_bams.txt
echo "condition" > DEseq_conditions.txt
echo "label" > DEseq_labels.txt

#loopStart:flag1
{
while read -r f1 f2 f3 f4; do

    [[ -n "$f3" ]] || continue #only submit jobs with sample name
    cd $wkdir/

    read1Files=$f1
    flag1="${read1Files%%,*}" # user the first file name as job flag
    read2Files=$f2
    samplename=$f3

    trimandfilterpath=$scriptsPath

    #rRNAbowtieout=${samplename}_bowtierRNA
    SpikeIn=${spikeinbase}/${genomeSpikeIn}/${genomeSpikeIn}
    #SpikeInbowtieout=${samplename}_bowtieSpikeIn
    Reference=${indexdir}/${STARversion}/${genomeRef}_${maxtrim}
    chrsizes=$scriptsPath/chromSizes/${genomeRef}.chrom.sizes
    #echo $read1Files $fileR2 $samplename

    IFS=',' read -r -a reads1 <<< $read1Files; IFS=',' read -r -a reads2 <<< $read2Files
    bams=""

    for i in "${!reads1[@]}"; do
        cd $wkdir
        r1=${reads1[i]}; r2=${reads2[i]}

        #@1,0,fastqc,,r1,sbatch -c 2 -p short -t 1:0:0 --mem 5G
        fastqc -t 2 $r1 $r2 -o fastqc/;

        #r1fastq=${r1##*/}; r2fastq=${r2##*/};
        #r1fastq=${r1fastq%.gz}; r2fastq=${r2fastq%.gz}
        r1fastq=${r1##*/}; r2fastq=${r2##*/};
        [[ $r1 == *.gz ]] && { r1fastq=${r1fastq%.gz}; r2fastq=${r2fastq%.gz}; }
        [[ $r1 == *.bz2 ]] && { r1fastq=${r1fastq%.bz2}; r2fastq=${r2fastq%.bz2}; }

        #@2,0,split,,r1,sbatch -c 2 -p short -t 1:0:0 --mem 5G
        set -x; zcat $r1 > fastq/$r1fastq || bzcat $r1 > fastq/$r1fastq; \
        zcat $r2 > fastq/$r2fastq || bzcat $r2 > fastq/$r2fastq; \
        rowCout=$(wc -l < fastq/$r1fastq);  rowCout=$((rowCout/4/2)); rowCout=$((rowCout*4)); \
        split -a 1 -l $rowCout fastq/$r1fastq fastq/$r1fastq; split -a 1 -l $rowCout fastq/$r2fastq fastq/$r2fastq; \
        [ `wc -l < fastq/${r1fastq}a` == `wc -l < fastq/${r2fastq}a` ] && [ `wc -l < fastq/${r1fastq}b` == `wc -l < fastq/${r2fastq}b` ]

        bams1=""
        for index in a b; do
            cd $wkdir
            fileR1=$r1fastq$index; fileR2=$r2fastq$index;
            in=$wkdir/fastq/$fileR1

            trimandfilterR1=$fileR1.1.trim_${mintrim}_${maxtrim}.minQS_20.fastq
            trimandfilterR2=$fileR1.2.trim_${mintrim}_${maxtrim}.minQS_20.fastq

            #@3,2,trimFilter,,in,sbatch -c 1 -p short -t 1:0:0 --mem 5G
            perl ${trimandfilterpath}/trim_and_filter_PE.pl -1 fastq/$fileR1 -2 fastq/$fileR2 -a ${mintrim} -b ${maxtrim} -c ${mintrim} -d ${maxtrim} -m 20 -q sanger -o fastq/$fileR1; ls fastq/$trimandfilterR1 fastq/$trimandfilterR2

            rRNAbowtieout=${fileR1}_bowtierRNA
            SpikeInbowtieout=${fileR1}_bowtieSpikeIn
            
            #exit 

            if [[ $useCutadapt == "yes" ]]; then
                cutadaptR1=$fileR1.trim.paired.fastq
                cutadaptR2=$fileR2.trim.paired.fastq

                #@4,3,cutadapt,,in,sbatch -c 1 -p short -t 1:0:0 --mem 5G
                cutadapt -f fastq --match-read-wildcards -m 20 -q 10 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o fastq/${cutadaptR1} -p fastq/${cutadaptR2} fastq/${trimandfilterR1} fastq/${trimandfilterR2} 2>&1 | tee fastq/${fileR1}_cutadaptLog.out; \
                rm fastq/${trimandfilterR1} fastq/${trimandfilterR2}

            else
                cutadaptR1=${trimandfilterR1}
                cutadaptR2=${trimandfilterR2}
            fi

            cd bowtierRNA

            rnaIndex=${rRNA}/rRNA${genomeRef}

            #@5,3.4,bowtie,rnaIndex,in,sbatch -c 5 -p short -t 1:0:0 --mem 40G
            bowtie -v2 -X1000 -p 5 --best -3 1 ${rRNA}/rRNA${genomeRef} -1 ../fastq/${cutadaptR1} -2 ../fastq/${cutadaptR2} ${rRNAbowtieout} 2>&1 | tee ${fileR1}_bowtierRNALog.out

            cd ../SpikeIn

            #@6,3.4,bowtie,SpikeIn,in,sbatch -c 5 -p short -t 1:0:0 --mem 40G
            bowtie -n2 -l 40 -X1000 -p 5 --best -3 1 --un "${fileR1}_NOT_${genomeSpikeIn}" ${SpikeIn} -1 ../fastq/${cutadaptR1} -2 ../fastq/${cutadaptR2} ${SpikeInbowtieout} 2>&1 | tee ${fileR1}_bowtieSpikeInLog.out; \
             echo $(grep "at least one" ${fileR1}_bowtieSpikeInLog.out | cut -d" " -f9) > ${fileR1}_spikeccount

            #@7,6,bzip2,,in,sbatch -c 1 -p short -t 1:0:0 --mem 4G
            rm ../fastq/${cutadaptR1}.bz2 ../fastq/${cutadaptR2}.bz2 2>/dev/null || true; \
            bzip2 ../fastq/${cutadaptR1}; bzip2 ../fastq/${cutadaptR2}

            cd ../mapping

            #@8,7,mapping,Reference,in,sbatch -c 5 -p short -t 12:0:0 --mem 45G
            STAR --runThreadN 5 --runMode alignReads --genomeDir ${Reference} --readFilesIn ../SpikeIn/${fileR1%.fastq*}_1.fastq${fileR1#*.fastq}_NOT_${genomeSpikeIn} ../SpikeIn/${fileR1%.fastq*}_2.fastq${fileR1#*.fastq}_NOT_${genomeSpikeIn} --quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 42949672960 --outFileNamePrefix $fileR1 --outMultimapperOrder Random --outSAMattrIHstart 0 --outFilterType BySJout --outFilterMismatchNmax ${numbermismatches} --alignSJoverhangMin 8 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outWigType bedGraph --outWigNorm None --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0; \
            rm ../SpikeIn/${fileR1%.fastq*}_1.fastq${fileR1#*.fastq}*.bz2 ../SpikeIn/${fileR1%.fastq*}_2.fastq${fileR1#*.fastq}*.bz2 2>/dev/null || true; \
            bzip2 ../SpikeIn/${fileR1%.fastq*}_1.fastq${fileR1#*.fastq}_NOT_${genomeSpikeIn}; bzip2 ../SpikeIn/${fileR1%.fastq*}_2.fastq${fileR1#*.fastq}_NOT_${genomeSpikeIn}; \
            samtools index -@ 5 ${fileR1}Aligned.sortedByCoord.out.bam; \
            mkdir -p Dedup; cd Dedup; \
            STAR --runThreadN 5 --limitBAMsortRAM 42949672960 --runMode inputAlignmentsFromBAM --bamRemoveDuplicatesType UniqueIdenticalNotMulti --inputBAMfile ../${fileR1}Aligned.sortedByCoord.out.bam --outFileNamePrefix $fileR1.flagdedup --outSAMtype BAM SortedByCoordinate; \
            samtools view -b -F 0x400 $fileR1.flagdedupProcessed.out.bam > $fileR1.dedup.Processed.out.bam; \
            samtools index -@ 5 $fileR1.dedup.Processed.out.bam

            bams1="$bams1 Dedup/$fileR1.dedup.Processed.out.bam"

        done

        cd $wkdir/mapping
        #@9,8,merge,,,sbatch -c 1 -p short -t 2:0:0 --mem 10G
        samtools merge -f $r1fastq.dedup.Processed.out.bam $bams1
        bams="$bams $r1fastq.dedup.Processed.out.bam"
    done
    #@10,9,merge1,,,sbatch -c 5 -p short -t 12:0:0 --mem 45G
    samtools merge -f $samplename.dedup.Processed.out.bam $bams; \
    samtools index -@ 5 $samplename.dedup.Processed.out.bam; \
    samtools view -c $samplename.dedup.Processed.out.bam chrM > ../$samplename.chrM.counts; \
    mkdir -p bedGraphsDedup; cd bedGraphsDedup; \
    STAR --runThreadN 5 --runMode inputAlignmentsFromBAM --inputBAMfile ../$samplename.dedup.Processed.out.bam --outFileNamePrefix $samplename.dedup --outWigType bedGraph --outWigNorm None --outWigStrand Stranded; \
    mv $samplename.dedupSignal.UniqueMultiple.str1.out.bg ../../raw_coverage_tracks/${samplename}_rawcounts_dedup_coverage_minus.bedGraph; \
    mv $samplename.dedupSignal.UniqueMultiple.str2.out.bg ../../raw_coverage_tracks/${samplename}_rawcounts_dedup_coverage_plus.bedGraph;

    cd $wkdir/raw_coverage_tracks

    #Raw counts tracks (convert into bigwig and then create 25bp binned data)
    #@11,10,sortAndGraph,,,sbatch -c 1 -p short -t 2:0:0 --mem 8G
    sort -k 1,1 -k 2,2n ${samplename}_rawcounts_dedup_coverage_minus.bedGraph > ${samplename}_sorted_dedup_coverage_minus.bedGraph; \
    bedGraphToBigWig ${samplename}_sorted_dedup_coverage_minus.bedGraph ${chrsizes} ${samplename}_rawcounts_dedup_coverage_minus.bw; \
    sort -k 1,1 -k 2,2n ${samplename}_rawcounts_dedup_coverage_plus.bedGraph > ${samplename}_sorted_dedup_coverage_plus.bedGraph; \
    bedGraphToBigWig ${samplename}_sorted_dedup_coverage_plus.bedGraph $chrsizes ${samplename}_rawcounts_dedup_coverage_plus.bw; \
    sh $scriptsPath/RNA_stats.sh $samplename

    [[ $f4 == empty ]] && exit

    cd $wkdir

    echo "mapping/$samplename.dedup.Processed.out.bam" >> DEseq_bams.txt
    echo $f4 >> DEseq_conditions.txt
    echo $samplename >> DEseq_labels.txt
done
} < "$manifest"

paste DEseq_bams.txt DEseq_conditions.txt | paste - DEseq_labels.txt > DEseq_metadata.txt

##@12,11,deseq,,,sbatch -p short -t 30 --mem 8G -c 5
#[ -d DEseqOutput ] && rm -r DEseqOutput; Rscript $scriptsPath/deseq_analysis.r DEseq_metadata.txt ${gtf} TRUE 2 DEseq 2 0.001 '${numerator_cond}' '${denominator_cond}' 2>&1 | tee rscript.out; \
#$scriptsPath/STAR_stats_join_size_factors.sh DEseqOutput/DEseq_DEseq_size_factors.txt
