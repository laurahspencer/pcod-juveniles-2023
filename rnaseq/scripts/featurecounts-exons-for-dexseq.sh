#!/bin/bash

#SBATCH --job-name=featurecounts
#SBATCH --output=/home/lspencer/pcod-juv-temp/aligned/featurecounts-exon-for-dexseq.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 20
#SBATCH -t 5-0:0:0

module load bio/subread/2.0.3
module load bio/samtools/1.11

IN=/home/lspencer/pcod-juv-temp/aligned
GTF=/home/lspencer/references/pcod-ncbi

# Exon Counts for DEXSeq
# Summarize paired-end reads and count fragments (instead of reads) which have both ends successfully aligned, i.e. don't count singletons, and don't count chimeras
# NOTE: I'm only counting reads that overlap with exons, and also summarizing them at the exon level (by specifying -f i'm summarizing at the "feature" level, i.e. exon)
featureCounts \
-p --countReadPairs \
--primary \
-B \
-T 20 \
-C \
-O \
-t exon \
-f \
-F gtf \
-a ${GTF}/GCF_031168955.1_ASM3116895v1_genomic_flat.gtf \
-o ${IN}/featurecounts_exon_dexseq \
${IN}/*.Aligned.sortedByCoord.out.bam

