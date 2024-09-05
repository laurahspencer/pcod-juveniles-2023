#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=0-05:00:00
#SBATCH --job-name=pcod-all_pca
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --output=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/wgsassign/err-out/pcod-refs-and-exp_wholegenome_pca_%A.out

module unload bio/pcangsd/0.99
module load bio/pcangsd/0.99
source /opt/bioinformatics/venv/pcangsd-0.99/bin/activate

pcangsd.py -threads 10 -beagle /home/lspencer/pcod-lcwgs-2023/analysis-20240606/wgsassign/refs-and-exp.beagle.gz -o /home/lspencer/pcod-lcwgs-2023/analysis-20240606/wgsassign/pca/pcod-refs-and-exp_wholegenome-pca -sites_save -pcadapt
