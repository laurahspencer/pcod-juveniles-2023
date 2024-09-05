#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=0-20:00:00
#SBATCH --job-name=pcod-refs_wgp-admix
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --output=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/job_outfiles/pcod-refs_wholegenome_wgassign_admix_%A-%a.out

#SBATCH --array=1-10%12

module unload bio/ngsadmix
module load bio/ngsadmix

for k_val in {1..10}
do
	if [[ ${SLURM_ARRAY_TASK_ID} == ${k_val} ]]; then
		break
	fi
done

NGSadmix -likes /home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/gls_wgassign/pcod-refs_wholegenome_wgassign_filtered.beagle.gz -K ${k_val} -outfiles /home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/admixture/pcod-refs_wholegenome_wgassign_filtered_k${k_val}-0 -P 10 -minMaf 0
NGSadmix -likes /home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/gls_wgassign/pcod-refs_wholegenome_wgassign_filtered.beagle.gz -K ${k_val} -outfiles /home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/admixture/pcod-refs_wholegenome_wgassign_filtered_k${k_val}-1 -P 10 -minMaf 0
NGSadmix -likes /home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/gls_wgassign/pcod-refs_wholegenome_wgassign_filtered.beagle.gz -K ${k_val} -outfiles /home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/admixture/pcod-refs_wholegenome_wgassign_filtered_k${k_val}-2 -P 10 -minMaf 0
