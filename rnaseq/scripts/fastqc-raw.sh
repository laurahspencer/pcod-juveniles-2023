#!/bin/bash

#SBATCH --job-name=pcod_fastqc-raw
#SBATCH --output=/home/lspencer/pcod-juv-temp/raw-data/sbatch-fastqc-raw.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 10

# This script is for running FastQC and MultiQC on raw (but concatenated)
# RNASeq data from Pacific cod juvenile tissues

# Load modules and virtual environments
module load bio/fastqc
source /home/lspencer/venv/bin/activate

# run fastqc on each raw read file for data from lane "5010"
cd /home/lspencer/pcod-juv-temp/raw-data/

fastqc \
--threads 8 \
*.fastq.gz \
--outdir /home/lspencer/pcod-juv-temp/raw-data/

# Run multiqc to summarize fastqc reports
multiqc \
/home/lspencer/pcod-juv-temp/raw-data/ \
--outdir /home/lspencer/pcod-juv-temp/raw-data/
