#!/bin/bash

#SBATCH --job-name=blast
#SBATCH --output=/home/lspencer/references/pcod-ncbi/blast-cds.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 20
#SBATCH -t 21-0:0:0

# This script is for blasting the Pacific cod coding sequences against
# the Uniqprot/Swissprot database

module load bio/blast/2.11.0+
module load bio/bedtools/2.29.2

BASE=/home/lspencer/references/pcod-ncbi/

# Blast cds for genes against uniprot/swissprot
blastx \
-query ${BASE}/cds_from_genomic.fna \
-db /home/lspencer/references/blast/uniprot_sprot_20220111_protein \
-out ${BASE}/cds_from_genomic_blastx.tab \
-evalue 1E-5 \
-num_threads 20 \
-outfmt 6
