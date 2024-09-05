#!/bin/bash

#SBATCH --job-name=trim
#SBATCH --cpus-per-task=4
#SBATCH --output=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/job_outfiles/pcod-exp_trimming_%A-%a.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --time=0-12:00:00
#SBATCH --array=1-157%48

module unload bio/trimmomatic/0.39 bio/fastp/0.23.2
module load bio/trimmomatic/0.39 bio/fastp/0.23.2
JOBS_FILE=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/scripts/pcod-exp_trimARRAY_input.txt
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

java -jar ${TRIMMOMATIC} PE -threads 4 -phred33 ${fq_r1} ${fq_r2} /home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/trimmed/${sample_id}_trimmed_R1_paired.fq.gz /home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/trimmed/${sample_id}_trimmed_R1_unpaired.fq.gz /home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/trimmed/${sample_id}_trimmed_R2_paired.fq.gz /home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/trimmed/${sample_id}_trimmed_R2_unpaired.fq.gz ILLUMINACLIP:/home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/adapters.txt:2:30:10:1:true MINLEN:40
fastp --trim_poly_g -L -A --cut_right -i /home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/trimmed/${sample_id}_trimmed_R1_paired.fq.gz -o /home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/trimmed/${sample_id}_trimmed_clipped_R1_paired.fq.gz -I /home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/trimmed/${sample_id}_trimmed_R2_paired.fq.gz -O /home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/trimmed/${sample_id}_trimmed_clipped_R2_paired.fq.gz -h /home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/trimmed/${sample_id}_trimmed_clipped_paired_report.html