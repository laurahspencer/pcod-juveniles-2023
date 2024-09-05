#!/bin/bash

#SBATCH --job-name=dexseq
#SBATCH --output=/home/lspencer/pcod-juv-temp/dexseq/dexseq-log.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -p himem
#SBATCH -c 24
#SBATCH -t 21-0:0:0

module load R
cd /home/lspencer/pcod-juv-temp/dexseq/

R CMD BATCH /home/lspencer/pcod-juv-temp/dexseq/DEXSeq.R /home/lspencer/pcod-juv-temp/dexseq/dexseq-output.txt
