#!/bin/sh

# deactivate base env
source deactivate

#install:
# most of software are available as modules
module load gcc/6.2.0 bedtools/2.27.1 fastqc/0.11.9 R/3.3.3 samtools/1.9 star/2.7.0a ucsc-tools/363 perl/5.24.0

# perl pacakges are manually installed using cpan
# perl=5.10.1
#     Statistics::Basic 1.6611
#     Statistics::Distributions 1.02
#     Statistics::R 0.34
eval `perl -Mlocal::lib=/n/data1/cores/bcbio/eclip/perl5`

# These softare are installed using conda. But python2 does not work with umi_tools, so we switch to python 3.7
# python=2.7.16
#     pysam=0.15.4
#     numpy=1.16.5
#     seaborn
# umi_tools=1.0.0
export PATH=/n/data1/cores/bcbio/eclip/fastq-tools-0.8.3/bin:/n/data1/cores/bcbio/eclip/eCLIP/bin:$PATH
module load conda/miniforge3/24.11.3-0
conda activate /n/data1/cores/bcbio/eclip/eclipEnv

# these two are manually installed following the instruction on software page

# The following was copied from: https://github.com/YeoLab/eCLIP/blob/master/documentation/eCLIP_analysisSOP_v2.2.1.docx
# Supplementary Protocol 2: eCLIP-seq Processing Pipeline

# Requires:

# 	•	Sequencing data (in fastq format).
# 	•	Genome STAR index directory (fasta files may be downloaded from UCSC; hg19)
# 	•	Repeat element STAR index directory (fasta files may be downloaded from RepBase)
# 	•	FASTA file containing barcodes for demultiplexing reads (described below)
# 	•	chrom.sizes file (tabbed file containing chromosome name and length, can be downloaded from UCSC; hg19)

# Installation (Programs Used & Version Information):
# This pipeline can be run with the following software installed on your environment. We recommend installation through the Anaconda package manager. Complete documentation can be found here: https://github.com/yeolab/eclip

# Yeo Lab Custom Script Versions:
# makebigwigfiles.py: https://github.com/YeoLab/makebigwigfiles/
# clipper: https://github.com/YeoLab/clipper
# clip_analysis: https://github.com/YeoLab/clipper/releases/tag/1.1
# demux.py: https://github.com/YeoLab/eclipdemux
# Input normalization and IDR workflow (described below): https://github.com/YeoLab/merge_peaks/releases/tag/0.0.6

# Other programs used:
# FastQC: v. 0.10.1
# Cutadapt: v. 1.14
# STAR: v. STAR_2.5.2b
# Samtools: v. 1.6
# bedToBigBed: v. 2.7
# Bedtools: v. 2.27.1
# R: v. 3.3.2
# fastq-sort: v. 0.8.0
# umi_tools: v. 1.0.0

# Python 2.7.17 :: Anaconda 2.1.0 (64-bit)
# Pysam 0.12.0.1
# Bx 0.5.0
# HTSeq 0.11.2
# Numpy 1.10.2
# Pandas 0.23.1
# Pybedtools 0.8.1
# Sklearn 0.17.1
# Scipy 0.17.1
# Matplotlib 2.2.2
# Gffutils 0.8.7.1
# Statsmodels 0.6.1
# perl=5.10.1 (changes to sorting in 5.22 may work but cause slightly different peak output)
# Statistics::Basic 1.6611
# Statistics::Distributions 1.02
# Statistics::R 0.34

# Outline of workflow:

# 	•	(paired-end only) Demultiplexes paired-end reads using inline barcodes and extracts unique molecular identifiers (UMI) with eclipdemux. (single-end) extract unique molecular identifiers with umi_tools
# 	•	Trims adapters & adapter-dimers with cutadapt
# 	•	Maps to repeat elements with STAR and filter
# 	•	Maps filtered reads to genome with STAR
# 	•	Removes PCR-duplicates with umi_tools (single-end) or with a custom python script (barcodecollapsepe.py)
# 	•	(paired-end only) Merges multiple inline barcodes and filters R1 (uses only R2 for peak calling)
# 	•	Calls enriched peak regions (peak clusters) with CLIPPER
# 	•	Uses size-matched input sample to normalize and calculate fold-change enrichment within enriched peak regions with custom perl scripts (overlap_peakfi_with_bam_PE.pl, peakscompress.pl)
# Usage:

# Script Details
# Our processing pipeline is performed using Common Workflow Language (CWL) definition files that take sequencing reads and returns per-replicate peaks that can be used as inputs into an IDR-based workflow to merge replicates. Complete documentation can be found here: https://github.com/yeolab/eclip

# Human Readable Description of Steps
# Note: For paired-end data, until the merging step each script is run twice, once for each barcode used
export PS4='Line $LINENO: '

#set -x

#set -e

outDir="`pwd`/eclipOut"
#loopStart:dir
for dir1 in `ls -v -d smartSlurmInputs/*/`; do
    echo working on dir1:  $dir1
    dir1=${dir1%/}
    dir=`basename $dir1`

    for dir2 in `ls -dr $dir1/input/  $dir1/ip/ | xargs -n 1 basename`; do
        echo working on dir2: $dir2
        pwd
        ls  $dir1/$dir2/*_1.fastq* 2>/dev/null || ls $dir1/$dir2/*_1.fq* 2>/dev/null || { echo Read file not found for $dir2! Please make sure the fastq files are named as xxx_1.fastq, xxx_2.fastq or xxx_1.fq, xxx_2.fq;  exit 1; }

        bamList=""
        for r1 in `ls $dir1/$dir2/*_1.fastq* $dir1/$dir2/*_1.fq* 2>/dev/null | xargs -n 1 basename`; do
            #echo r1 is $r1
            readgroup=${r1%_*}
            barcodeID=${readgroup#*.}
            barcodeA=${barcodeID%-*}
            barcodeB=${barcodeID#*-}
            barcodeOption="--expectedbarcodeida $barcodeA --expectedbarcodeidb $barcodeB"
            barcodeIDs="$barcodeA $barcodeB"
            [[ "$barcodeA" == "$barcodeB" ]] && barcodeOption="--expectedbarcodeida $barcodeA --expectedbarcodeidb C08" && barcodeIDs="$barcodeA"

            echo working on readgroup: $readgroup
            r2=${r1%_*}_2${r1##*_1}

            mkdir -p $outDir/$dir/$dir2/$readgroup
            cd $outDir/$dir/$dir2/$readgroup

            # Identify unique molecular identifiers (UMIs) (SE): Use umi_tools to extract unique molecular barcodes.

            # umi_tools extract \
            # --random-seed 1 \
            # --bc-pattern NNNNNNNNNN \
            # --log EXAMPLE_SE.rep1_clip.---.--.metrics \
            # --stdin $r1 \
            # --stdout EXAMPLE_SE.rep1.umi.r1.fq

            #Demultiplexing inline barcodes and identify UMIs (PE): Perform demultiplexing using a supplied barcodes FASTA file. Depending on protocol, --length may be longer than 5.
            input=$outDir/../$dir1/$dir2/$r1

            #@1,0,demux,,input,sbatch -c 1 -p short -t 2:0:0 --mem 8G
            demux \
            --metrics EXAMPLE_PE.rep1_clip.---.--.metrics \
            $barcodeOption \
            --fastq_1 $outDir/../$dir1/$dir2/$r1 \
            --fastq_2 $outDir/../$dir1/$dir2/$r2 \
            --newname rep2_clip \
            --dataset EXAMPLE_PE \
            --barcodesfile /n/data1/cores/bcbio/eclip/eCLIP/example/inputs/yeolabbarcodes_20170101.fasta \
            --length 5
#exit
            # 	•	Use the FASTA ID for --expectedbarcodeida and --expectedbarcodeidb corresponding to each expected barcode of your IP dir2. Use "NIL" for barcode-less size-matched input
            # 	•	--length describes the length of the UMI to save.
            # 	•	--metrics, --newname, --dataset are all labels that can be used to customize output file names
            # 	•	--barcodesfile contents are described as follows:

            # Inline barcode description:
            # Each inline barcode is ligated to the 5’ end of Read1 and its id and sequence are listed below:
            # A01    ATTGCTTAGATCGGAAGAGCGTCGTGT
            # B06    ACAAGCCAGATCGGAAGAGCGTCGTGT
            # C01    AACTTGTAGATCGGAAGAGCGTCGTGT
            # D8f    AGGACCAAGATCGGAAGAGCGTCGTGT
            # A03    ANNNNGGTCATAGATCGGAAGAGCGTCGTGT
            # G07    ANNNNACAGGAAGATCGGAAGAGCGTCGTGT
            # A04    ANNNNAAGCTGAGATCGGAAGAGCGTCGTGT
            # F05    ANNNNGTATCCAGATCGGAAGAGCGTCGTGT
            # RiL19/NIL   AGATCGGAAGAGCGTCGTGT
            # (see eCLIP protocol document for full description of these oligos)

            # We have observed occasional double ligation events on the 5’ end of Read1, and we have found that to fix this requires we run cutadapt twice.  Additionally, because two adapters are used for each library (to ensure proper balancing on the Illumina sequencer), two separate barcodes may be ligated to the same Read1 5’ end (often with 5’ truncations).  To fix this we split the barcodes up into 15bp chunks so that cutadapt is able to deconvolute barcode adapters properly (as by default it will not find adapters missing the first N bases of the adapter sequence)

            # The barcodes file is made by appending one of the barcodes below (these are the same barcode sequences used to demultiplex). See example file (yeolabbarcodes_20170101.fasta):
            # AAGCAAT A01
            # GGCTTGT B06
            # ACAAGTT C01
            # TGGTCCT D8f
            # ATGACCNNNNT  A03
            # TCCTGTNNNNT  G07
            # CAGCTTNNNNT  A04
            # GGATACNNNNT  F05

            # To the 5’ adapter
            # CTTCCGATCT




            for barcodeID in $barcodeIDs; do

                # Cutadapt round 1: Takes output from demultiplexed files.  Run to trim off both 5’ and 3’ adapters on both reads

                #Fastqc round 1: Run and examined by eye to make sure libraries look alright

                # Cutadapt round 2: Takes output from cutadapt round 1. Run to trim off the 3’ adapters on read 2, to control for double ligation events.

                #Fastqc round 2: Takes output from STAR rmRep.  Runs a second round of fastqc to verify that after read grooming the data still is usable.
                #Fastq-sort: Takes cutadapt round 2 output and sorts to reduce randomness

                #@2,1,cutadapt,,,sbatch -c 1 -p short -t 2:0:0 --mem 8G
                sh -c "module load python/2.7.12 cutadapt/1.14; cutadapt \
                -f fastq \
                --match-read-wildcards \
                --times 1 \
                -e 0.1 \
                -O 1 \
                --quality-cutoff 6 \
                -m 18 \
                -a NNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
                -g CTTCCGATCTACAAGTT \
                -g CTTCCGATCTTGGTCCT \
                -A AACTTGTAGATCGGA \
                -A AGGACCAAGATCGGA \
                -A ACTTGTAGATCGGAA \
                -A GGACCAAGATCGGAA \
                -A CTTGTAGATCGGAAG \
                -A GACCAAGATCGGAAG \
                -A TTGTAGATCGGAAGA \
                -A ACCAAGATCGGAAGA \
                -A TGTAGATCGGAAGAG \
                -A CCAAGATCGGAAGAG \
                -A GTAGATCGGAAGAGC \
                -A CAAGATCGGAAGAGC \
                -A TAGATCGGAAGAGCG \
                -A AAGATCGGAAGAGCG \
                -A AGATCGGAAGAGCGT \
                -A GATCGGAAGAGCGTC \
                -A ATCGGAAGAGCGTCG \
                -A TCGGAAGAGCGTCGT \
                -A CGGAAGAGCGTCGTG \
                -A GGAAGAGCGTCGTGT \
                -o EXAMPLE_PE.rep2_clip.$barcodeID.r1.fqTr.fq \
                -p EXAMPLE_PE.rep2_clip.$barcodeID.r2.fqTr.fq \
                EXAMPLE_PE.rep2_clip.$barcodeID.r1.fq.gz \
                EXAMPLE_PE.rep2_clip.$barcodeID.r2.fq.gz";  \
                fastqc -t 2 --extract -k 7 EXAMPLE_PE.rep2_clip.$barcodeID.r1.fqTr.fq -o .; \
                fastqc -t 2 --extract -k 7 EXAMPLE_PE.rep2_clip.$barcodeID.r2.fqTr.fq -o . ; \
                sh -c "module load python/2.7.12 cutadapt/1.14; cutadapt \
                -f fastq \
                --match-read-wildcards \
                --times 1 \
                -e 0.1 \
                -O 5 \
                --quality-cutoff 6 \
                -m 18 \
                -A AACTTGTAGATCGGA \
                -A AGGACCAAGATCGGA \
                -A ACTTGTAGATCGGAA \
                -A GGACCAAGATCGGAA \
                -A CTTGTAGATCGGAAG \
                -A GACCAAGATCGGAAG \
                -A TTGTAGATCGGAAGA \
                -A ACCAAGATCGGAAGA \
                -A TGTAGATCGGAAGAG \
                -A CCAAGATCGGAAGAG \
                -A GTAGATCGGAAGAGC \
                -A CAAGATCGGAAGAGC \
                -A TAGATCGGAAGAGCG \
                -A AAGATCGGAAGAGCG \
                -A AGATCGGAAGAGCGT \
                -A GATCGGAAGAGCGTC \
                -A ATCGGAAGAGCGTCG \
                -A TCGGAAGAGCGTCGT \
                -A CGGAAGAGCGTCGTG \
                -A GGAAGAGCGTCGTGT \
                -o EXAMPLE_PE.rep2_clip.$barcodeID.r1.fqTrTr.fq \
                -p EXAMPLE_PE.rep2_clip.$barcodeID.r2.fqTrTr.fq \
                EXAMPLE_PE.rep2_clip.$barcodeID.r1.fqTr.fq \
                EXAMPLE_PE.rep2_clip.$barcodeID.r2.fqTr.fq; "; \
                fastqc -t 2 --extract -k 7 EXAMPLE_PE.rep2_clip.$barcodeID.r1.fqTrTr.fq -o .; \
                fastqc -t 2 --extract -k 7 EXAMPLE_PE.rep2_clip.$barcodeID.r2.fqTrTr.fq -o . ; \
                fastq-sort --id EXAMPLE_PE.rep2_clip.$barcodeID.r1.fqTrTr.fq > EXAMPLE_PE.rep2_clip.$barcodeID.r1.fqTrTr.sorted.fq; \
                fastq-sort --id EXAMPLE_PE.rep2_clip.$barcodeID.r2.fqTrTr.fq > EXAMPLE_PE.rep2_clip.$barcodeID.r2.fqTrTr.sorted.fq

                #STAR rmRep: Takes sorted output from cutadapt round 2.  Maps to human specific version of RepBase used to remove repetitive elements, helps control for spurious artifacts from rRNA (& other) repetitive reads.

                #Fastq-sort: Takes unmapped output from STAR rmRep and sorts it to account for issues with STAR not outputting first and second mate pairs in order

                rIndex=$outDir/../hg113seqs_repbase_starindex
                input=EXAMPLE_PE.rep2_clip.$barcodeID.r1.fqTrTr.sorted.fq
                #@3,2,star,rIndex,input,sbatch -c 2 -p short -t 2:0:0 --mem 20G
                STAR \
                --runMode alignReads \
                --runThreadN 2 \
                --genomeDir /n/data1/cores/bcbio/eclip/eCLIP/example/inputs/hg113seqs_repbase_starindex \
                --genomeLoad NoSharedMemory \
                --alignEndsType EndToEnd \
                --outSAMunmapped Within \
                --outFilterMultimapNmax 30 \
                --outFilterMultimapScoreRange 1 \
                --outFileNamePrefix EXAMPLE_PE.rep2_clip.$barcodeID.r1.fqTrTr.sorted.STAR \
                --outSAMtype BAM Unsorted \
                --outFilterType BySJout \
                --outBAMcompression 10 \
                --outReadsUnmapped Fastx \
                --outFilterScoreMin 10 \
                --outSAMattrRGline ID:foo \
                --outSAMattributes All \
                --outSAMmode Full \
                --outStd Log \
                --readFilesIn EXAMPLE_PE.rep2_clip.$barcodeID.r1.fqTrTr.sorted.fq EXAMPLE_PE.rep2_clip.$barcodeID.r2.fqTrTr.sorted.fq; \
                mv EXAMPLE_PE.rep2_clip.$barcodeID.r1.fqTrTr.sorted.STARAligned.out.bam EXAMPLE_PE.rep2_clip.$barcodeID.r1.fq.repeat-mapped.bam; \
                mv EXAMPLE_PE.rep2_clip.$barcodeID.r1.fqTrTr.sorted.STARUnmapped.out.mate1 EXAMPLE_PE.rep2_clip.$barcodeID.r1.fq.repeat-unmapped.fq; \
                mv EXAMPLE_PE.rep2_clip.$barcodeID.r1.fqTrTr.sorted.STARUnmapped.out.mate2 EXAMPLE_PE.rep2_clip.$barcodeID.r2.fq.repeat-unmapped.fq; \
                cat EXAMPLE_PE.rep2_clip.$barcodeID.r1.fqTrTr.sorted.STARLog.out; \
                fastq-sort --id EXAMPLE_PE.rep2_clip.$barcodeID.r1.fq.repeat-unmapped.fq > EXAMPLE_PE.rep2_clip.$barcodeID.r1.fq.repeat-unmapped.sorted.fq; \
                fastq-sort --id EXAMPLE_PE.rep2_clip.$barcodeID.r2.fq.repeat-unmapped.fq > EXAMPLE_PE.rep2_clip.$barcodeID.r2.fq.repeat-unmapped.sorted.fq

                #STAR genome mapping: Takes output from STAR rmRep.  Maps unique reads to the human genome
                #Name sort BAM: sort output from STAR by name to ensure read pairs are adjacent.
                #Barcode_collapse_pe (PE): takes output from STAR genome mapping.  Custom random-mer-aware script for PCR duplicate removal.
                #Position sort BAM:  Takes output from barcode collapse PE (or from SE namesort bam).  Sorts resulting bam file for use downstream.

                gIndex=hg19chr19kbp550_starindex
                input=EXAMPLE_PE.rep2_clip.$barcodeID.r1.fq.repeat-unmapped.sorted.fq
                #@4,3,star1,gIndex,,sbatch -c 2 -p short -t 2:0:0 --mem 20G
                STAR \
                --runMode alignReads \
                --runThreadN 2 \
                --genomeDir /n/data1/cores/bcbio/eclip/eCLIP/example/inputs/hg19chr19kbp550_starindex \
                --genomeLoad NoSharedMemory \
                --readFilesIn \
                EXAMPLE_PE.rep2_clip.$barcodeID.r1.fq.repeat-unmapped.sorted.fq \
                EXAMPLE_PE.rep2_clip.$barcodeID.r2.fq.repeat-unmapped.sorted.fq \
                --outSAMunmapped Within \
                --outFilterMultimapNmax 1 \
                --outFilterMultimapScoreRange 1 \
                --outFileNamePrefix EXAMPLE_PE.rep2_clip.$barcodeID.r1.fq.repeat-unmapped.sorted.STAR \
                --outSAMattributes All \
                --outSAMtype BAM Unsorted \
                --outFilterType BySJout \
                --outReadsUnmapped Fastx \
                --outFilterScoreMin 10 \
                --outSAMattrRGline ID:foo \
                --outStd Log \
                --alignEndsType EndToEnd \
                --outBAMcompression 10 \
                --outSAMmode Full; \
                mv EXAMPLE_PE.rep2_clip.$barcodeID.r1.fq.repeat-unmapped.sorted.STARAligned.out.bam EXAMPLE_PE.rep2_clip.$barcodeID.r1.fq.genome-mapped.bam; \
                cat EXAMPLE_PE.rep2_clip.$barcodeID.r1.fq.repeat-unmapped.sorted.STARLog.out; \
                samtools sort -n -o EXAMPLE_PE.rep2_clip.$barcodeID.r1.fq.genome-mappedSo.bam EXAMPLE_PE.rep2_clip.$barcodeID.r1.fq.genome-mapped.bam; \
                barcodecollapsepe.py \
                -o EXAMPLE_PE.rep2_clip.$barcodeID.r1.fq.genome-mappedSo.rmDup.bam \
                -m EXAMPLE_PE.rep2_clip.$barcodeID.r1.fq.genome-mappedSo.rmDup.metrics \
                -b EXAMPLE_PE.rep2_clip.$barcodeID.r1.fq.genome-mappedSo.bam 2>&1 > $barcodeID.barcodecollapsepe.log; \
                samtools sort -o EXAMPLE_PE.rep2_clip.$barcodeID.r1.fq.genome-mappedSo.rmDupSo.bam EXAMPLE_PE.rep2_clip.$barcodeID.r1.fq.genome-mappedSo.rmDup.bam; \
                samtools index EXAMPLE_PE.rep2_clip.$barcodeID.r1.fq.genome-mappedSo.rmDupSo.bam

                #[ -f $outDir/$dir/$dir2/$readgroup/EXAMPLE_PE.rep2_clip.$barcodeID.r1.fq.genome-mappedSo.rmDupSo.bam ] &&

                bamList="$bamList $outDir/$dir/$dir2/$readgroup/EXAMPLE_PE.rep2_clip.$barcodeID.r1.fq.genome-mappedSo.rmDupSo.bam"
                #Barcode_collapse_se (SE): takes output from STAR genome mapping. Use umi_tools dedup to identify the extracted random-mer from the previous step and perform PCR duplicate removal.

                #umi_tools dedup \
                #--random-seed 1 \
                #---I EXAMPLE_SE.rep1_clip.umi.r1.fq.genome-mappedSoSo.bam \
                #--method unique \
                #--output-stats EXAMPLE_SE.rep1_clip.umi.r1.fq.genome-mappedSoSo.txt \
                #-S EXAMPLE_SE.rep1_clip.umi.r1.fq.genome-mappedSoSo.rmDup.bam
                #break
            done
            cd -
        done

        cd $outDir/$dir/$dir2/

        #Samtools  merge (PE only): Takes inputs from multiple final bam files.  Merges the two technical replicates for further downstream analysis.

        #samtools  merge EXAMPLE_PE.rep2_clip.$barcodeID.r1.fq.genome-mappedSo.rmDupSo.merged.bam EXAMPLE_PE.rep2_clip.$barcodeID.r1.fq.genome-mappedSo.rmDupSo.bam EXAMPLE_PE.rep2_clip.D8f.r1.fq.genome-mappedSo.rmDupSo.bam; \

        #Samtools view (PE only):  Takes output from samtools merge.  Only outputs the second read in each pair for use with a single stranded peak caller.  This is the final bam file to perform analysis on.

        #Make normalized read density bigwig files: Takes input from samtools view.  Makes bw files to be uploaded to the genome browser or for other visualization. Use --direction f for SE clip as reads are not reversed.

        #@5,4,samtools,,,sbatch -c 1 -p short -t 2:0:0 --mem 8G
        samtools  merge -f EXAMPLE_PE.rep2_clip.r1.fq.genome-mappedSo.rmDupSo.merged.bam $bamList; \
        samtools  index  EXAMPLE_PE.rep2_clip.r1.fq.genome-mappedSo.rmDupSo.merged.bam; \
        samtools view -f 128 -b -o EXAMPLE_PE.rep2_clip.r1.fq.genome-mappedSo.rmDupSo.merged.r2.bam EXAMPLE_PE.rep2_clip.r1.fq.genome-mappedSo.rmDupSo.merged.bam; \
        makebigwigfiles \
        --bw_pos EXAMPLE_PE.rep2_clip.r1.fq.genome-mappedSo.rmDupSo.merged.r2.norm.pos.bw \
        --bw_neg EXAMPLE_PE.rep2_clip.r1.fq.genome-mappedSo.rmDupSo.merged.r2.norm.neg.bw \
        --bam EXAMPLE_PE.rep2_clip.r1.fq.genome-mappedSo.rmDupSo.merged.r2.bam \
        --genome /n/data1/cores/bcbio/eclip/eCLIP/example/inputs/hg19chr19kbp550_starindex/chrNameLength.txt \
        --direction r


        #@6,5,samtools1,,,sbatch -c 1 -p short -t 2:0:0 --mem 8G
        samtools view -cF 4 EXAMPLE_PE.rep2_clip.r1.fq.genome-mappedSo.rmDupSo.merged.r2.bam > mapped_readnum.txt

        cd -
    done

    #Input normalization: Compares the number of reads within the IP dir2 to the number of reads within the size-matched INPUT dir2 across Clipper-called peak clusters. This step is performed both within this pipeline as well as within the merge_peaks pipeline using the same perl scripts.

    #Clipper:  Takes results from samtools view.  Calls peaks on those files.

    #@7,6,cliper,,,sbatch -c 12 -p short -t 12:0:0 --mem 20G
    sh -c "source deactivate; conda activate /n/data1/cores/bcbio/eclip/clipperEnv; \
    clipper --processors=12 --quiet --species hg19 \
    --bam $outDir/$dir/ip/EXAMPLE_PE.rep2_clip.r1.fq.genome-mappedSo.rmDupSo.merged.r2.bam \
    --save-pickle \
    --outfile $outDir/$dir/ip/EXAMPLE_PE.rep2_clip.r1.fq.genome-mappedSo.rmDupSo.merged.r2.peakClusters.bed"

    #EXAMPLE_PE.rep2_clip.r1.fq.genome-mappedSo.rmDupSo.merged.r2.bam \
    #EXAMPLE_PE.rep2_input.NIL.r1.fq.genome-mappedSo.rmDupSo.r2.bam \
    ##@8,7,overlap,,,sbatch -c 1 -p short -t 2:0:0 --mem 8G
    #overlap_peakfi_with_bam_PE.pl $outDir/$dir/ip/EXAMPLE_PE.rep2_clip.r1.fq.genome-mappedSo.rmDupSo.merged.r2.bam \
    #$outDir/$dir/input/EXAMPLE_PE.rep2_clip.r1.fq.genome-mappedSo.rmDupSo.merged.r2.bam \
    #$outDir/$dir/ip/EXAMPLE_PE.rep2_clip.r1.fq.genome-mappedSo.rmDupSo.merged.r2.peakClusters.bed \
    #$outDir/$dir/ip/mapped_readnum.txt $outDir/$dir/input/mapped_readnum.txt \
    #$outDir/$dir/EXAMPLE_PE.rep2_clip.r1.fq.genome-mappedSo.rmDupSo.merged.r2.peakClusters.normed.bed; \
    #rm $outDir/$dir/EXAMPLE_PE.rep2_clip.r1.fq.genome-mappedSo.rmDupSo.merged.r2.peakClusters.normed.compressed.bed.full; \
    #compress_l2foldenrpeakfi_for_replicate_overlapping_bedformat_outputfull.pl \
    #$outDir/$dir/EXAMPLE_PE.rep2_clip.r1.fq.genome-mappedSo.rmDupSo.merged.r2.peakClusters.normed.bed.full \
    #$outDir/$dir/EXAMPLE_PE.rep2_clip.r1.fq.genome-mappedSo.rmDupSo.merged.r2.peakClusters.normed.compressed.bed  \
    #$outDir/$dir/EXAMPLE_PE.rep2_clip.r1.fq.genome-mappedSo.rmDupSo.merged.r2.peakClusters.normed.compressed.bed.full; \
    #make_informationcontent_from_peaks.pl EXAMPLE_PE.rep2_clip.r1.fq.genome-mappedSo.rmDupSo.merged.r2.peakClusters.normed.compressed.bed.full $outDir/$dir/ip/mapped_readnum.txt $outDir/$dir/input/mapped_readnum.txt  entropy.full entropy.excessreads

    #note: .compressed.be is a bed6 formatted bed file with the following columns: chromosome, start, end,
    # -log10 pvalueð Þ, log2 fold enrichmentð Þ, strand. For each peak, we recommend instituting a filter for both
    # -log10 pvalueð Þ and log2 fold enrichmentð Þ to be 3 or greater (pvalue d 0.001 and fold enrichment e 8).
    # Note that this step will also produce a .full file that will contain additional information corresponding
    # to each peak. The columns are described from left to right for each CLIPper peak cluster:
    #  Chromosome
    #  Start
    #  End
    #  Peak name
    #  Number of IP-aligned reads in peak
    #  Number of SMInput-aligned reads in peak
    #  Fisher-exact or chi-square P value
    #  Chi-value (if chi-square test used) or (F)isher exact test or depleted (DEPL) if IP-aligned reads SMInput-aligned reads
    #  Test applied: (C)hi-square test or (F)isher exact test
    #  enriched if IP-aligned reads > SMInput-aligned reads; depleted otherwise
    #  -pvalue
    #  fold enrichment

done


