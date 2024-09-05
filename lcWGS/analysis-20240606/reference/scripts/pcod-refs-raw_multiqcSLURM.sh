#!/bin/bash

#SBATCH -p himem
#SBATCH --job-name=multiQC
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --output=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/job_outfiles/pcod-refs-raw_multiQC.out

source /home/ltimm/bin/hydraQC/bin/activate
multiqc /home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/fastqc/raw/
