#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=0-05:00:00
#SBATCH --job-name=pcod-refs_wgp-pca
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --output=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/job_outfiles/pcod-refs_wholegenome_polymorphic_%A.out

module unload bio/pcangsd/0.99
module load bio/pcangsd/0.99
source /opt/bioinformatics/venv/pcangsd-0.99/bin/activate

pcangsd.py -threads 10 -beagle /home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/gls_wgassign/pcod-refs_wholegenome_wgassign_filtered.beagle.gz -o /home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/pca/pcod-refs_wholegenome-wgassign -sites_save -pcadapt
