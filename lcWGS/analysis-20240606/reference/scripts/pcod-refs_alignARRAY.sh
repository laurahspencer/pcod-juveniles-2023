#!/bin/bash

#SBATCH --job-name=align
#SBATCH --cpus-per-task=10
#SBATCH --output=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/job_outfiles/pcod-refs_alignment_%A-%a.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --time=7-00:00:00
#SBATCH --array=1-655%48

module unload aligners/bwa/0.7.17 bio/samtools/1.11 bio/bamtools/2.5.1 bio/picard/2.23.9 bio/bamutil/1.0.5
module load aligners/bwa/0.7.17 bio/samtools/1.11 bio/bamtools/2.5.1 bio/picard/2.23.9 bio/bamutil/1.0.5

JOBS_FILE=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/scripts/pcod-refs_alignARRAY_input.txt
IDS=$(cat ${JOBS_FILE})

for sample_line in ${IDS}
do
	job_index=$(echo ${sample_line} | awk -F ":" '{print $1}')
	fq_r1=$(echo ${sample_line} | awk -F ":" '{print $2}')
	fq_r2=$(echo ${sample_line} | awk -F ":" '{print $3}')
	if [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then
		break
	fi
done

sample_id=$(echo $fq_r1 | sed 's!^.*/!!')
sample_id=${sample_id%%_*}

bwa mem -M -t 10 /home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/bwa/GCF_031168955.1_ASM3116895v1_genomic ${fq_r1} ${fq_r2} 2> /home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/bwa/pcod-refs_${sample_id}_bwa-mem.out > /home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/bamtools/pcod-refs_${sample_id}.sam

samtools view -bS -F 4 /home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/bamtools/pcod-refs_${sample_id}.sam > /home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/bamtools/pcod-refs_${sample_id}.bam
rm /home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/bamtools/pcod-refs_${sample_id}.sam

samtools view -h /home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/bamtools/pcod-refs_${sample_id}.bam | samtools view -buS - | samtools sort -o /home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/bamtools/pcod-refs_${sample_id}_sorted.bam
rm /home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/bamtools/pcod-refs_${sample_id}.bam

java -jar $PICARD MarkDuplicates I=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/bamtools/pcod-refs_${sample_id}_sorted.bam O=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/bamtools/pcod-refs_${sample_id}_sorted_dedup.bam M=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/bamtools/pcod-refs_${sample_id}_dups.log VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true
rm /home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/bamtools/pcod-refs_${sample_id}_sorted.bam

bam clipOverlap --in /home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/bamtools/pcod-refs_${sample_id}_sorted_dedup.bam --out /home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/bamtools/pcod-refs_${sample_id}_sorted_dedup_clipped.bam --stats
rm /home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/bamtools/pcod-refs_${sample_id}_sorted_dedup.bam

samtools depth -aa /home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/bamtools/pcod-refs_${sample_id}_sorted_dedup_clipped.bam | cut -f 3 | gzip > /home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/bamtools/pcod-refs_${sample_id}.depth.gz

samtools index /home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/bamtools/pcod-refs_${sample_id}_sorted_dedup_clipped.bam
