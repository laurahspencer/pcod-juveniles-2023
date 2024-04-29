#!/bin/bash

#SBATCH --job-name=zp3
#SBATCH --output=/home/lspencer/pcod-juv-temp/aligned/zp3-region.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 20
#SBATCH -t 5-0:0:0

module load bio/samtools/1.11

cd /home/lspencer/pcod-juv-temp/aligned/

# # index all .bam files
# SAMPLES=$(ls *.Aligned.sortedByCoord.out.bam)
# for sample in ${SAMPLES}
# do
# samtools index ${sample}
# done

#Pull data for the ZP3 gene region, sort and index resulting .bams.
SAMPLES=$(ls *.Aligned.sortedByCoord.out.bam | sed -e 's/.Aligned.sortedByCoord.out.bam//')
for sample in ${SAMPLES}
do
samtools view -b -h ${sample}.Aligned.sortedByCoord.out.bam "NC_082390.1:1867790-1875256" > ${sample}.zp3.bam
samtools sort ${sample}.zp3.bam > ${sample}.zp3.sorted.bam
samtools index ${sample}.zp3.sorted.bam
done

#Merge the ZP3 alignment data for all fish into one file
samtools merge merged.zp3.bam *.zp3.sorted.bam
samtools sort merged.zp3.bam > merged.zp3.sorted.bam
samtools index merged.zp3.sorted.bam
