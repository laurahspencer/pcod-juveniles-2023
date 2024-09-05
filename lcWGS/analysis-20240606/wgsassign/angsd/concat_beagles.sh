#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --job-name=concat-beagles
#SBATCH --mail-type=ALL
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --output=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/wgsassign/err-out/concat-beagles.txt

base=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/wgsassign/angsd

zcat ${base}/pcod-_Chr1_wgassign.beagle.gz | head -n 1 > ${base}/pcod_wholegenome_wgassign.beagle; for i in ${base}/pcod-_Chr1_wgassign.beagle.gz ${base}/pcod-_Chr2_wgassign.beagle.gz ${base}/pcod-_Chr3_wgassign.beagle.gz ${base}/pcod-_Chr4_wgassign.beagle.gz ${base}/pcod-_Chr5_wgassign.beagle.gz ${base}/pcod-_Chr6_wgassign.beagle.gz ${base}/pcod-_Chr7_wgassign.beagle.gz ${base}/pcod-_Chr8_wgassign.beagle.gz ${base}/pcod-_Chr9_wgassign.beagle.gz ${base}/pcod-_Chr10_wgassign.beagle.gz ${base}/pcod-_Chr11_wgassign.beagle.gz ${base}/pcod-_Chr12_wgassign.beagle.gz ${base}/pcod-_Chr13_wgassign.beagle.gz ${base}/pcod-_Chr14_wgassign.beagle.gz ${base}/pcod-_Chr15_wgassign.beagle.gz ${base}/pcod-_Chr16_wgassign.beagle.gz ${base}/pcod-_Chr17_wgassign.beagle.gz ${base}/pcod-_Chr18_wgassign.beagle.gz ${base}/pcod-_Chr19_wgassign.beagle.gz ${base}/pcod-_Chr20_wgassign.beagle.gz ${base}/pcod-_Chr21_wgassign.beagle.gz ${base}/pcod-_Chr22_wgassign.beagle.gz ${base}/pcod-_Chr23_wgassign.beagle.gz ${base}/pcod-_Chr24_wgassign.beagle.gz; do zcat $i | tail -n +2 -q >> ${base}/pcod_wholegenome_wgassign.beagle; done
gzip ${base}/pcod_wholegenome_wgassign.beagle

