#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --job-name=multiQC
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --output=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/job_outfiles/pcod-exp-raw_multiQC.out

source /home/lspencer/venv/bin/activate
multiqc /home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/fastqc/raw/
