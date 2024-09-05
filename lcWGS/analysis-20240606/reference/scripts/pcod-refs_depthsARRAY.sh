#!/bin/bash

#SBATCH --job-name=depth
#SBATCH --cpus-per-task=5
#SBATCH --output=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/job_outfiles/pcod-refsdepths_%A-%a.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --time=0-12:00:00
#SBATCH --array=1-655%24

JOBS_FILE=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/scripts/pcod-refs_depthsARRAY_input.txt
IDS=$(cat ${JOBS_FILE})

for sample_line in ${IDS}
do
	job_index=$(echo ${sample_line} | awk -F ":" '{print $1}')
	depth_file=$(echo ${sample_line} | awk -F ":" '{print $2}')
	if [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then
		break
	fi
done

touch /home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/bamtools/pcod-refs_depths.csv
/home/lspencer/lcWGS-pipeline/mean_cov_ind.py -i ${depth_file} -o /home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/bamtools/pcod-refs_depths.csv
