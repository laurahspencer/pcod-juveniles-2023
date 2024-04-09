#!/bin/bash

#SBATCH --job-name=fastqc
#SBATCH --output=/home/lspencer/pcod-juv-temp/trimmed-data/sbatch-fastqc.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 20
#SBATCH -t 5-0:0:0

# This script is for running fastqc and multqc on trimmed data

module load bio/fastqc
source /home/lspencer/venv/bin/activate

IN=/home/lspencer/pcod-juv-temp/raw-data
OUT=/home/lspencer/pcod-juv-temp/trimmed-data

SAMPLES=$(ls ${IN}/*_R1_001.fastq.gz | \
awk -F "/" '{print $NF}' | \
awk -F "." '{print $1}' | \
sed -e 's/_R1_001//')

for sample in ${SAMPLES}
do
# Run fastqc on trimmed data files
    fastqc \
        --threads 10 \
        -o ${OUT} \
        ${OUT}/${sample}.trimmed.R1.fastq.gz \
        ${OUT}/${sample}.trimmed.R2.fastq.gz \
        &> ${OUT}/fastqc.${sample}.log
done

# Run multiqc to summarize fastqc reports
multiqc \
${OUT} \
--outdir ${OUT}
