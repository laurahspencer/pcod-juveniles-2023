#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --job-name=bwa_index_GCF_031168955.1_ASM3116895v1_genomic
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --output=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/job_outfiles/bwa-index_GCF_031168955.1_ASM3116895v1_genomic.out

module unload aligners/bwa/0.7.17
module load aligners/bwa/0.7.17

bwa index -p /home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/bwa/GCF_031168955.1_ASM3116895v1_genomic /home/lspencer/references/pcod-ncbi/GCF_031168955.1_ASM3116895v1_genomic.fa
