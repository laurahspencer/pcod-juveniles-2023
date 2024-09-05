#!/bin/bash

#SBATCH --job-name=snp-testing
#SBATCH --output=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/wgsassign/err-out/snp-testing-loos-region-one-beagle.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 20
#SBATCH -t 2-0:0:0

source ~/.bashrc
mamba activate WGSassign

base=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/wgsassign/snp-testing/by-marine-region
mkdir ${base}/using-one-beagle/testing-LOOs
ids=${base}/test-IDs2.txt
numbers=(10 50 100 500 1000 5000 10000 25000 50000)

for n in "${numbers[@]}"
do

	# ex. amre.testing.ind85.ds_2x.sites-filter.top_100000_each.beagle.gz
	input_beagle=${base}/using-one-beagle/testing-top-${n}.beagle.gz
	outname=${base}/using-one-beagle/testing-LOOs/top-${n}

	# Get likelihoods for leave-one-out assignment within known reference populations
	# Output = 1) reference.popAF.npy, 2) reference.pop_like_LOO.txt
	WGSassign --beagle ${input_beagle} --pop_af_IDs ${ids} \
	--get_reference_af --loo --out ${outname} --threads 20

done

# Separately run this on all sites (not part of for loop)
WGSassign --beagle ${base}/using-one-beagle/testing-all-sites.beagle.gz --pop_af_IDs ${ids} \
--get_reference_af --loo --out ${base}/using-one-beagle/testing-LOOs/all-sites --threads 20
