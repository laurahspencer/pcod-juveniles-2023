#!/bin/bash

#SBATCH --job-name=pcod_trim
#SBATCH --output=/home/lspencer/pcod-juv-temp/trimmed-data/sbatch-trim.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 20
#SBATCH -t 5-0:0:0

# This script is for trimming raw (but concatenated) RNA-Seq data and
# filtering for quality and length.
# It has been slightly adapted from code written by Giles Goetz

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
    # Trimming the Illumina adapters
    # Quality-trim  5’ end with cutoff=15 & 3’ end with cutoff=10
    # Trimming out leftover N's
    # Filtering out sequences shorter then 50bp
    cutadapt \
        -o ${OUT}/${sample}.trimmed.R1.fastq.gz \
        -p ${OUT}/${sample}.trimmed.R2.fastq.gz \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -q 15,10 \
        -m 50 \
        --trim-n \
        --cores=20 \
        ${IN}/${sample}_R1_001.fastq.gz \
        ${IN}/${sample}_R2_001.fastq.gz \
        &> ${OUT}/cutadapt.${sample}.log

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
