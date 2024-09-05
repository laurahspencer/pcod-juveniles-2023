#!/bin/bash

#SBATCH --job-name=fqc_array_pcod-exp
#SBATCH --cpus-per-task=1
#SBATCH --output=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/job_outfiles/pcod-exp-raw_fastqc_%A-%a.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --time=0-03:00:00
#SBATCH --array=1-314%24

module unload bio/fastqc/0.11.9
module load bio/fastqc/0.11.9

JOBS_FILE=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/scripts/pcod-exp-raw_fqcARRAY_input.txt
IDS=$(cat ${JOBS_FILE})

for sample_line in ${IDS}
do
	job_index=$(echo ${sample_line} | awk -F ":" '{print $1}')
	fq=$(echo ${sample_line} | awk -F ":" '{print $2}')
	if [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then
		break
	fi
done

fastqc ${fq} -o /home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/fastqc/raw/