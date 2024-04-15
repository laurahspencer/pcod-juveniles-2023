#!/bin/bash

#SBATCH --job-name=featurecounts
#SBATCH --output=/home/lspencer/pcod-juv-temp/aligned/featurecounts.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 20
#SBATCH -t 5-0:0:0

module load bio/subread/2.0.3
module load bio/samtools/1.11

IN=/home/lspencer/pcod-juv-temp/aligned
GTF=/home/lspencer/references/pcod-ncbi

# Gene Counts
# Summarize paired-end reads and count fragments (instead of reads) which have both ends successfully aligned, i.e. don't count singletons, and don't count chimeras
# NOTE - I'm only counting reads that overlap with exons, then summarizing them at the gene level (here we use the db_xref, which is the ncbi gene ID (ie Entrez Gene ID), as they are unique gene identifiers)

featureCounts \
-p --countReadPairs \
--primary \
-B \
-T 20 \
-C \
-O \
-t exon \
-g db_xref \
-a ${GTF}/GCF_031168955.1_ASM3116895v1_genomic_edited.gtf \
-o ${IN}/featurecounts_gene \
${IN}/*.Aligned.sortedByCoord.out.bam

# Exon Counts
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
-g ncbi_id_exon \
-a ${GTF}/GCF_031168955.1_ASM3116895v1_genomic_edited.gtf \
-o ${IN}/featurecounts_exon \
${IN}/*.Aligned.sortedByCoord.out.bam
