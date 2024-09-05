#!/bin/bash
#SBATCH --job-name=prune-LD
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --output=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/job_outfiles/LD-prune.txt
#SBATCH --time=7:00:00

module load bio/prune_graph

BASEDIR=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference
infile=$BASEDIR/gls_wgassign/ngsLD/pcod-refs_wholegenome_LD
outfile=$BASEDIR/gls_wgassign/ngsLD/pcod-refs_wholegenome_unlinked

prune_graph --in ${infile} --weight-field column_7 --weight-filter "column_7 >= 0.5" --out ${outfile} --verbose
