#!/bin/bash

if [ "$1" == "-h" ]; then
  echo "Usage:
bash `basename $0` samplename genomeSpike genomeRef outputDir noDedup

   -samplename            samplename (e.g. Pol2_ControlA)
   -genomeSpike           bowtie2 index (e.g. dm6)
   -genomeRef             bowtie2 index (e.g. hg38, mm10)
   -outputDir             Output directory (e.g. scratch folder)
   -noDedup				  'true' if samples were not UMI deduplicated"
  exit 0
fi

samplename=$1
genomeSpike=$2
genomeRef=$3
outputDir=$4
noDedup=$5

repoPath="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
if [[ $repoPath =~ /n/data1/cores/ntc ]]; then
	scriptsPath=/n/data1/cores/ntc/scripts
elif [[ $repoPath =~ /n/data1/hms/bcmp/adelman ]]; then
	scriptsPath=/n/data1/hms/bcmp/adelman/Scripts/ntc
else
	exit 1
fi

if [[ "$noDedup" != "true" ]]; then
	dedup_naming=_dedup
else
	dedup_naming=
fi

# trimmed reads
total=$(tr -s ' ' <${outputDir}/logs/${samplename}_cutadaptLog.out | grep "Total read pairs" | cut -d' ' -f5 | sed 's/,//g')
trim=$(tr -s ' ' <${outputDir}/logs/${samplename}_cutadaptLog.out | grep "Pairs written" | cut -d' ' -f5 | sed 's/,//g')
trim_perc=$(printf %.2f $(echo "($trim/$total) * 100" | bc -l))

# spike reads
spike=$(samtools view -F 4 -c ${outputDir}/mapping/${samplename}_${genomeSpike}.bam)
spike=$(printf %.0f $(echo "$spike/2" | bc -l))
if [[ "$noDedup" != "true" ]]; then
	spike_dedup=$(samtools view -F 4 -c ${outputDir}/mapping/${samplename}_${genomeSpike}_dedup.bam)
	spike_dedup=$(printf %.0f $(echo "$spike_dedup/2" | bc -l))
	dedup_perc_spike=$(printf %.2f $(echo "(1-($spike_dedup/$spike)) * 100" | bc -l))
	spike_perc=$(printf %.2f $(echo "($spike_dedup/$trim) * 100" | bc -l))
else
	spike_perc=$(printf %.2f $(echo "($spike/$trim) * 100" | bc -l))
fi

# ref reads
ref=$(samtools view -F 4 -c ${outputDir}/mapping/${samplename}_${genomeRef}.bam)
ref=$(printf %.0f $(echo "$ref/2" | bc -l))
if [[ "$noDedup" != "true" ]]; then
	ref_dedup=$(samtools view -F 4 -c ${outputDir}/mapping/${samplename}_${genomeRef}_dedup.bam)
	ref_dedup=$(printf %.0f $(echo "$ref_dedup/2" | bc -l))
	dedup_perc_ref=$(printf %.2f $(echo "(1-($ref_dedup/$ref)) * 100" | bc -l))
	ref_perc=$(printf %.2f $(echo "($ref_dedup/$trim) * 100" | bc -l))
else
	ref_perc=$(printf %.2f $(echo "($ref/$trim) * 100" | bc -l))
fi

tss_prox_list=${scriptsPath}/tss/${genomeRef}/${genomeRef}.basic_+-150_uniq.txt
if [[ "$genomeRef" == "mm9" ]] || [[ "$genomeRef" == "dm3" ]]; then
	tss_prox_list=${tss_prox_list/.basic/}
elif [[ "$genomeRef" == "fm_library" ]]; then
	tss_prox_list=${tss_prox_list/.basic_+-150/_-125+60}
fi

${scriptsPath}/AdelmanLab/NIH_scripts/make_heatmap/make_heatmap -t 6 -b v -l s -s s --nohead -p ${outputDir}/bedGraphs/${samplename}_${genomeRef}${dedup_naming}_3pr_forward.bedGraph -m ${outputDir}/bedGraphs/${samplename}_${genomeRef}${dedup_naming}_3pr_reverse.bedGraph $tss_prox_list ${outputDir}/logs/${samplename}_${genomeRef}${dedup_naming}_3pr_+-150.txt 1
tss=$(awk -F'\t' '{total+=$7}END{print total}' ${outputDir}/logs/${samplename}_${genomeRef}${dedup_naming}_3pr_+-150.txt)

# insert size
insert=$(printf %.0f $(grep -A 1 "MEAN_INSERT_SIZE" ${outputDir}/logs/${samplename}_${genomeRef}${dedup_naming}_picardoutput.txt | tail -n 1 | cut -f5))

if [[ "$noDedup" != "true" ]]; then
	tss_perc=$(printf %.2f $(echo "($tss/$ref_dedup) * 100" | bc -l))
	echo -e $samplename"\t"$trim"\t"$trim_perc"\t"$spike_dedup"\t"$spike_perc"\t"$ref_dedup"\t"$ref_perc"\t"$tss"\t"$tss_perc"\t"$insert"\t"$dedup_perc_spike"\t"$dedup_perc_ref > ${outputDir}/logs/$samplename.stats.txt
else
	tss_perc=$(printf %.2f $(echo "($tss/$ref) * 100" | bc -l))
	echo -e $samplename"\t"$trim"\t"$trim_perc"\t"$spike"\t"$spike_perc"\t"$ref"\t"$ref_perc"\t"$tss"\t"$tss_perc"\t"$insert > ${outputDir}/logs/$samplename.stats.txt
fi
