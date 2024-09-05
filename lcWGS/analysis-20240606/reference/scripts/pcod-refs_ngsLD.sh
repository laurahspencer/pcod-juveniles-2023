#!/bin/bash
#SBATCH --job-name=ngsLD-mem
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --output=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/job_outfiles/ngsLD.txt
#SBATCH -t 60:00:00
#SBATCH -p himem
#SBATCH -c 24

module load bio/ngsld

pos=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/gls_wgassign/ngsLD/pcod-refs_wholegenome_wgassign_4ngsLD.sites.gz
geno=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/gls_wgassign/ngsLD/pcod-refs_wholegenome_wgassign_4ngsLD.beagle.gz

outdir=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/gls_wgassign/ngsLD
outname=pcod-refs_wholegenome_LD

ngsLD --geno ${geno} --pos ${pos} --probs --n_ind 633 --n_sites 431691 --max_kb_dist 50 --n_threads 24 --out ${outdir}/${outname}
