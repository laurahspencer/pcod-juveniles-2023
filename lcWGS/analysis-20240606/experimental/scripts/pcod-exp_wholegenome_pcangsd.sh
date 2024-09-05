#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=0-05:00:00
#SBATCH --job-name=pcod-exp_wgp-pca
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --output=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/job_outfiles/pcod-exp_wholegenome_polymorphic_%A.out

module unload bio/pcangsd/0.99
module load bio/pcangsd/0.99
source /opt/bioinformatics/venv/pcangsd-0.99/bin/activate

pcangsd.py -threads 10 -beagle /home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/gls_wgassign/pcod-exp_wholegenome_wgassign.beagle.gz -o /home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/pca/pcod-exp_wholegenome-wgassign -sites_save -pcadapt
