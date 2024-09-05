#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --job-name=fai_GCF_031168955.1_ASM3116895v1_genomic
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --output=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/job_outfiles/fai_GCF_031168955.1_ASM3116895v1_genomic.out

module unload bio/samtools/1.11
module load bio/samtools/1.11

samtools faidx /home/lspencer/references/pcod-ncbi/GCF_031168955.1_ASM3116895v1_genomic.fa
