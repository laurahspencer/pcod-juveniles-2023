#!/bin/bash

#SBATCH --job-name=wgsassign
#SBATCH --output=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/wgsassign/err-out/pop-assign.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 20
#SBATCH -t 2-0:0:0

source ~/.bashrc
mamba activate WGSassign

base=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/wgsassign
best_snps_beagle=${base}/snp-testing/testing-top-5000.beagle.gz
best_snps_afs=${base}/snp-testing/testing-LOOs/top-5000.pop_af.npy
exp_beagle_all=${base}/join-beagles-temp/rehead_beagle2.gz  #this was a temporary file from my join-beagle.sh script that I kept! It has sample IDs in header row.  
exp_beagle_best=${base}/assignment/pcod-exp_top-5000.beagle
outname=${base}/assignment/pcod-experimental-assign_top-5k-snps

# Create sites file containing the best SNPs to use for population assignment
zcat $best_snps_beagle | cut -f1 > ${base}/assignment/assign.sites

# combine header row with filtered sites
zcat ${exp_beagle_all} | head -n 1 > ${exp_beagle_best}
awk 'NR==FNR{c[$1]++;next};c[$1]' ${base}/assignment/assign.sites <(zcat ${exp_beagle_all}) >> ${exp_beagle_best}
gzip ${exp_beagle_best}

# Generate file listing samples in order they appear in beagle file using sample IDs in header row
zcat ${exp_beagle_best}.gz | cut --complement -f1-3 | head -n 1 | tr '\t' '\n' | \
sed -e 's/_AA//g' -e 's/_AB//g' -e 's/_BB//g' | uniq > ${base}/assignment/pcod-exp-sample-order.txt

# Run assignment to identify source population of experimental fish 
WGSassign --beagle ${exp_beagle_best}.gz --pop_af_file ${best_snps_afs} --get_pop_like --out ${outname} --threads 20
