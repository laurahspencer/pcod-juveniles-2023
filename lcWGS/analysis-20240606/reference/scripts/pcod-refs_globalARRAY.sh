#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=7-00:00:00
#SBATCH --job-name=global_pcod-refs
#SBATCH --output=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/job_outfiles/pcod-refs_global_%A-%a.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --array=1-48%24

module unload bio/angsd/0.933
module load bio/angsd/0.933

JOBS_FILE=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/scripts/pcod-refs_angsdARRAY_input.txt
IDS=$(cat ${JOBS_FILE})

for sample_line in ${IDS}
do
	job_index=$(echo ${sample_line} | awk -F ":" '{print $1}')
	contig=$(echo ${sample_line} | awk -F ":" '{print $2}')
	if [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then
		break
	fi
done

angsd -b /home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/pcod-refs_filtered_bamslist.txt -ref /home/lspencer/references/pcod-ncbi/GCF_031168955.1_ASM3116895v1_genomic.fa -r ${contig}: -out /home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/gls/pcod-refs_${contig}_global -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -minMapQ 15 -minQ 15 -doCounts 1 -setminDepth 633 -setmaxDepth 12660.0 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doDepth 1 -dumpCounts 3 -only_proper_pairs 1
