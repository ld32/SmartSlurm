#!/bin/sh

#set -x

usage()
{
    echo "usage: $(basename "$0") [-h | --help]
            [-m | --metadata METADATA_CSV]
            [-p | --partition SBATCH_PARTITION]
            [-t | --jobTime SBATCH_TIME]
            [-c | --cpuMem SBATCH_MEM_PER_CPU]
            [-n | --noDedup TURN_OFF_UMI_DEDUP]

required arguments:
    -m, --metadata        CSV conforming to format described in README.md

optional arguments:
    -p, --partition        sbatch partition [default = medium]
    -t, --jobTime        total job time allotted [default = 0-20]
    -c, --cpuMem        memory per CPU allotted (total mem = c x 6) [default = 16G]
    -n, --noDedup        turn off UMI deduplication [default: UMI deduplicaton on]
    "
}

repoPath=/n/data1/cores/ntc #$(dirname $0)
if [[ $repoPath =~ /n/data1/cores/ntc ]]; then
   groups_path=/n/data1/cores/ntc
   scriptsPath=/n/data1/cores/ntc/scripts
   bwtPath=/n/data1/cores/ntc/scripts/bowtie2-indexes
   pythonPath=/n/data1/cores/ntc/ntc-conda
elif [[ $repoPath =~ /n/data1/hms/bcmp/adelman ]]; then
   groups_path=/n/data1/hms/bcmp/adelman
   scriptsPath=/n/data1/hms/bcmp/adelman/Scripts/ntc
   bwtPath=/n/data1/hms/bcmp/adelman/Scripts/bowtie2-indexes
   pythonPath=/n/data1/hms/bcmp/adelman/ntc-conda
else
   echo exit 1
fi

# groups_path=/home/ld32/scratch3/data/ntc
# scriptsPath=/home/ld32/scratch3/data/ntc/scripts
# bwtPath=/home/ld32/scratch3/data/ntc/scripts/bowtie2-indexes
# pythonPath=/home/ld32/scratch3/data/ntc/ntc-conda

no_dedup=""
while [ "$1" != "" ]; do
    case $1 in
        -m | --metadata    ) shift
                            metadata=$1
                            ;;
        -p | --partition ) shift
                            partition=$1
                            ;;
        -t | --jobTime ) shift
                            job_time=$1
                            ;;
        -c | --cpuMem ) shift
                            job_mem=$1
                            ;;
        -n | --noDedup )    no_dedup=true
                            ;;
        -h | --help )  usage
                       exit
                            ;;
        * )                    usage
                            exit 1
    esac
    shift
done

if [ -z "$metadata" ]; then
    exit 1
else
    grep -l '^M$' $metadata && dos2unix "$metadata"
fi

# load modules
module purge
module load conda2/4.2.13
module load gcc/6.2.0 python/2.7.12
unset PYTHONPATH
module load cutadapt/1.14
module load bowtie2/2.3.4.3
module load samtools/1.3.1
module load R/4.0.1
export R_LIBS="${groups_path}/R-4.0.1/library/"
module load picard/2.8.0
module load subread/1.6.2
# module load fastqc/0.11.9
# module load fastq_screen/0.11.2
# module load multiqc/1.5
# module load bowtie/1.2.2

# remove emtpy rows
sed -i '/^\s*$/d' "$metadata"

# read metadata and check to see if there is any error
{
read
while IFS=$',' read -r -a m
do
    sampleName=${m[0]}
    name=${m[1]}
    group=${m[2]}
    genomeRef=${m[3]}
    genomeSpike=${m[4]}
    owner=${m[5]}
    cellType=${m[6]}
    analysisType=${m[7]}
    tssList=${m[8]}
    outDir=${m[9]}
    fileR1=${m[10]}
    fileR2=${m[11]}

    [ -f "$fileR1" ] || { echo Read1 file $fileR1 does not exist; exit 1; }
    [ -f "$fileR2" ] || { echo Read2 file $fileR2 does not exist; exit 1; }
    [ -f "$tssList" ] || { echo tss file $tssList does not exist; exit 1; }

    mkdir -p $outDir/logs || { echo outDir $outDir could not be made; exit 1; }
    refIndex=$bwtPath/${genomeRef}_$genomeSpike/${genomeRef}_$genomeSpike
    bowtie2-inspect -n $refIndex >/dev/null || { echo bowtie2 index $refIndex does not exist or not correct; exit 1; }

    [[ "$analysisType" == "Start-seq" ]] || [[ "$analysisType" == "PRO-seq" ]] || { echo "ERROR: 'analysisType' must be one of: Start-seq, PRO-seq"; exit 1; }

done
} < "$metadata"

#loopStart:fileR1

# read metadata and submit sample-level pipeline jobs
{
read
while IFS=$',' read -r -a m
do
    sampleName=${m[0]}
    name=${m[1]}
    group=${m[2]}
    genomeRef=${m[3]}
    genomeSpike=${m[4]}
    owner=${m[5]}
    cellType=${m[6]}
    analysisType=${m[7]}
    tssList=${m[8]}
    outDir=${m[9]}
    outDir=${outDir%/} # trim ending /
    fileR1=${m[10]}
    fileR2=${m[11]}

    declare -A ends=(["r1"]="" ["r2"]="")

    if [[ "$analysisType" == "Start-seq" ]]; then
        ends["r1"]+="5pr"
        ends["r2"]+="3pr"
        endsR1="5pr"
        endsR2="3pr"
    elif [[ "$analysisType" == "PRO-seq" ]]; then
        refPrefix=${sampleName}_$genomeRef
        spikePrefix=${sampleName}_$genomeSpike
        endsR1="3pr"
        endsR2="5pr"
    else
        echo "ERROR: 'analysisType' must be one of: Start-seq, PRO-seq"
        exit 1
    fi

    # set species-specific chromosomes for filtering ref from spike reads
    declare -A chrs=(["hg38"]=`seq 1 22 | sed 's/^/chr/'`
        ["hg19"]=`seq 1 22 | sed 's/^/chr/'`
        ["mm39"]=`seq 1 19 | sed 's/^/chr/'`
        ["mm10"]=`seq 1 19 | sed 's/^/chr/'`
        ["mm9"]=`seq 1 19 | sed 's/^/chr/'`
        ["dm6"]="2L 2R 3L 3R 4 X Y"
        ["dm3"]="2L 2R 3L 3R 4 X"
        ["pJBmix1"]="pJB"
        ["XL9_2"]="chr1L chr2L chr1S chr2S chr5L chr6L chr3L chr4L chr5S chr6S chr7L chr4S chr3S chr8L chr9_10L chr9_10S chr8S chr7S"
        ["bosTau9"]=`seq 1 29 | sed 's/^/chr/'`
        ["hg38Akata"]=`seq 1 22 | sed 's/^/chr/'`
        ["hg38hsv17"]=`seq 1 22 | sed 's/^/chr/'`
        ["canFam4"]=`seq 1 38 | sed 's/^/chr/'`)

    chrs["hg38"]+=" chrX chrY"
    chrs["hg19"]+=" chrX chrY"
    chrs["mm39"]+=" chrX chrY"
    chrs["mm10"]+=" chrX chrY"
    chrs["mm9"]+=" chrX chrY"
    chrs["bosTau9"]+=" chrX"
    chrs["hg38Akata"]+=" chrX chrY LN824208.1"
    chrs["hg38hsv17"]+=" chrX chrY NC_001806.2"
    chrs["canFam4"]+=" chrX"   
    
    mkdir -p $outDir/fastq

    # remove adapters
    # (UMI-specific) trim UMI length off 3' ends to ensure no UMI read-through left over, (Both) trim 1bp off 3' ends so R2 always upstream of reverse complement of R1
    if [[ "$no_dedup" == "true" ]]; then

        mapInputR1=$outDir/fastq/$sampleName.1.noadap.fastq
        mapInputR2=$outDir/fastq/$sampleName.2.noadap.fastq

        #@1,0,cutadaptSeqtk,,fileR1.fileR2,sbatch -c 1 -p short -t 12:0:0 --mem 8G
        rm $mapInputR1 $mapInputR2 &>/dev/null || : ; cutadapt -f fastq -O 1 --match-read-wildcards -m 20 -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG -A GATCGTCGGACTGTAGAACTCTGAAC -o $outDir/fastq/$sampleName.1.trim.paired.fastq -p $outDir/fastq/$sampleName.2.trim.paired.fastq $fileR1 $fileR2 2>&1 | tee $outDir/logs/${sampleName}_cutadaptLog.out; \
        $scriptsPath/seqtk/seqtk trimfq -e 1 $outDir/fastq/$sampleName.1.trim.paired.fastq > $mapInputR1; cp $outDir/fastq/$sampleName.2.trim.paired.fastq $mapInputR2
    else

        mapInputR1=$outDir/fastq/$sampleName.1.noumi.fastq
        mapInputR2=$outDir/fastq/$sampleName.2.noumi.fastq

        # extract dual UMIs
        #@2,0,umiCutAdaptSeqtk,,fileR1.fileR2,sbatch -c 1 -p short -t 12:0:0 --mem 8G
        rm $mapInputR1 $mapInputR2 2>/dev/null || : ; source activate $pythonPath; \
        umi_tools extract -L $outDir/logs/${sampleName}_umi_extract.log -p NNNNNN --bc-pattern2=NNNNNN -I $fileR1 -S $outDir/fastq/${sampleName}_R1_UMI_extract.fastq.gz --read2-in=$fileR2 --read2-out=$outDir/fastq/${sampleName}_R2_UMI_extract.fastq.gz; \
        conda deactivate; \
        cutadapt -f fastq -O 1 --match-read-wildcards -m 26 -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG -A GATCGTCGGACTGTAGAACTCTGAAC -o $outDir/fastq/$sampleName.1.trim.paired.fastq -p $outDir/fastq/$sampleName.2.trim.paired.fastq $outDir/fastq/${sampleName}_R1_UMI_extract.fastq.gz $outDir/fastq/${sampleName}_R2_UMI_extract.fastq.gz 2>&1 | tee $outDir/logs/${sampleName}_cutadaptLog.out; \
        $scriptsPath/seqtk/seqtk trimfq -e 7 $outDir/fastq/$sampleName.1.trim.paired.fastq > $mapInputR1; \
        $scriptsPath/seqtk/seqtk trimfq -e 6 $outDir/fastq/$sampleName.2.trim.paired.fastq > $mapInputR2
    fi

#exit
    mkdir -p $outDir/mapping
    bowtieOut=$outDir/mapping/${refPrefix}_$genomeSpike.bam
    refIndex=$bwtPath/${genomeRef}_$genomeSpike/${genomeRef}_$genomeSpike
    #@3,1.2,bowtie2,refIndex,mapInputR1.mapInputR2,sbatch -c 6 -p short -t 12:0:0 --mem 12G
    rm $bowtieOut* &>/dev/null || : ; \
    bowtie2 -p 6 -x $refIndex -1 $mapInputR1 -2 $mapInputR2 2> $outDir/logs/${sampleName}_bowtie2.log | samtools view -u -f 0x2 | samtools sort -@ 6 -o $bowtieOut && cat $outDir/logs/${sampleName}_bowtie2.log || { cat $outDir/logs/${sampleName}_bowtie2.log; exit 1; }
#exit
    if [[ "$genomeRef" == "fm_library" ]]; then
    
        chrsr=`echo ${chrs[$genomeSpike]}`
        # filter ref from spike reads
        #@4,3,samtools4,,bowtieOut,sbatch -c 1 -p short -t 12:0:0 --mem 8G
        rm $outDir/mapping/$spikePrefix.bam.tmp* $outDir/mapping/${refPrefix}.bam.tmp* &>/dev/null || : ; \
        samtools index $bowtieOut; \
        samtools view -u $bowtieOut $chrsr  | samtools sort -@ 6 -o $outDir/mapping/$spikePrefix.bam; \
        samtools index $outDir/mapping/$spikePrefix.bam; \
        samtools view -u -L $scriptsPath/chromSizes/fm_library.bed $bowtieOut | samtools sort -@ 6 -o $outDir/mapping/${refPrefix}.bam; \
        samtools index $outDir/mapping/${refPrefix}.bam

    else
        chrsr=`echo ${chrs[$genomeSpike]}`
        chrsr1=`echo ${chrs[$genomeRef]}`
        # filter ref from spike reads
        #@5,3,samtools5,,bowtieOut,sbatch -c 1 -p short -t 12:0:0 --mem 8G
        rm $outDir/mapping/$spikePrefix.bam* $outDir/mapping/${refPrefix}.bam* &>/dev/null || echo; \
        samtools index $bowtieOut; \
        samtools view -u $bowtieOut $chrsr | samtools sort -@ 6 -o $outDir/mapping/$spikePrefix.bam; \
        samtools index $outDir/mapping/$spikePrefix.bam; \
        samtools view -u $bowtieOut $chrsr1 | samtools sort -@ 6 -o $outDir/mapping/${refPrefix}.bam; \
        samtools index $outDir/mapping/${refPrefix}.bam
    fi

    # UMI dedup
    if [[ "$no_dedup" != "true" ]]; then
        #@6,4.5,umiDedup6,,bowtieOut,sbatch -c 1 -p short -t 12:0:0 --mem 8G
        source activate $pythonPath; \
        rm $outDir/mapping/${spikePrefix}_dedup.bam* $outDir/mapping/${refPrefix}_dedup.bam* &>/dev/null || : ; \
        umi_tools dedup -L $outDir/logs/${spikePrefix}_dedup.log --paired -I $outDir/mapping/$spikePrefix.bam -S $outDir/mapping/${spikePrefix}_dedup.bam; \
        umi_tools dedup -L $outDir/logs/${refPrefix}_dedup.log --paired -I $outDir/mapping/${refPrefix}.bam -S $outDir/mapping/${refPrefix}_dedup.bam; \
        conda deactivate

        refPrefix=${refPrefix}_dedup
        spikePrefix=${spikePrefix}_dedup
    fi

    # calculate insert size distributions
    #@7,4.5.6,picard7.NoCheckpoint,,,sbatch -c 1 -p short -t 12:0:0 --mem 8G
    java -jar $PICARD/picard-2.8.0.jar CollectInsertSizeMetrics I=$outDir/mapping/$spikePrefix.bam H=$outDir/logs/${spikePrefix}_picardhistogram.pdf O=$outDir/logs/${spikePrefix}_picardoutput.txt 2>&1 | tee $outDir/logs/${spikePrefix}_picardLog.txt; \
    java -jar $PICARD/picard-2.8.0.jar CollectInsertSizeMetrics I=$outDir/mapping/${refPrefix}.bam H=$outDir/logs/${refPrefix}_picardhistogram.pdf O=$outDir/logs/${refPrefix}_picardoutput.txt 2>&1 | tee $outDir/logs/${refPrefix}_picardLog.txt

    # separate ref 3/5pr reads
    #@8,4.5.6,samtools8.NoCheckpoint,,bowtieOut,sbatch -c 1 -p short -t 12:0:0 --mem 8G
    rm $outDir/mapping/${refPrefix}_${endsR1}.sam* $outDir/mapping/${refPrefix}_${endsR2}.sam* &>/dev/null || : ; \
    samtools view -u -f 0x40 $outDir/mapping/${refPrefix}.bam | samtools sort -n -@ 6 -O SAM -o $outDir/mapping/${refPrefix}_${endsR1}.sam; \
    samtools view -u -f 0x80 $outDir/mapping/${refPrefix}.bam | samtools sort -n -@ 6 -O SAM -o $outDir/mapping/${refPrefix}_${endsR2}.sam

    # bedgraphs
    mkdir -p $outDir/bedGraphs

    #@9,8,bedGraphs9,,bowtieOut,sbatch -c 1 -p short -t 12:0:0 --mem 140G
    perl $scriptsPath/AdelmanLab/NIH_scripts/bowtie2stdbedgraph/bowtie2stdBedGraph.pl -x -a -o b -b 10 -D -l $scriptsPath/chromSizes/$genomeRef.chrom.sizes $outDir/mapping/${refPrefix}_3pr.sam $outDir/bedGraphs/${refPrefix}_3pr; \
    perl $scriptsPath/AdelmanLab/NIH_scripts/bowtie2stdbedgraph/bowtie2stdBedGraph.pl -a -o b -b 10 -D -l $scriptsPath/chromSizes/$genomeRef.chrom.sizes $outDir/mapping/${refPrefix}_5pr.sam $outDir/bedGraphs/${refPrefix}_5pr

    # matrix
    mkdir -p $outDir/matrix
    tssPrefix=$( echo $tssList | awk '{split(a[split($0,a,"/")],b,".txt"); print b[1]}' )

    #@10,9,heatmap10,,,sbatch -c 1 -p short -t 1:0:0 --mem 4G
    rm $outDir/logs/$sampleName.stats.txt &>/dev/null || : ; \
    $scriptsPath/AdelmanLab/NIH_scripts/make_heatmap/make_heatmap -t 6 -l s -s s --nohead -p $outDir/bedGraphs/${refPrefix}_3pr_forward.bedGraph -m $outDir/bedGraphs/${refPrefix}_3pr_reverse.bedGraph -- $tssList $outDir/matrix/${cellType}_${tssPrefix}_${refPrefix}_3pr_25mer_+-2kb.txt -2000 25 160;\
    $scriptsPath/AdelmanLab/NIH_scripts/make_heatmap/make_heatmap -t 6 -l s -s s --nohead -p $outDir/bedGraphs/${refPrefix}_5pr_forward.bedGraph -m $outDir/bedGraphs/${refPrefix}_5pr_reverse.bedGraph -- $tssList $outDir/matrix/${cellType}_${tssPrefix}_${refPrefix}_5pr_5mer_+-500.txt -500 5 200; \
    sh /n/data1/cores/ntc/pipeline2_smartSlurm/stats.sh $sampleName $genomeSpike $genomeRef $outDir $no_dedup $scriptsPath

    #exit
done
} < "$metadata"

# project-level stats
#@11,10,multistats,,,sbatch -p short -t 30 --mem=8G
/n/data1/cores/ntc/pipeline2_smartSlurm/mergeStats.sh $metadata $no_dedup; \
Rscript $scriptsPath/NascentTranscriptionCore/pipeline2/multiStats.R $metadata
