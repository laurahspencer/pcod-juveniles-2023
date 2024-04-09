#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --partition=himem
#SBATCH --time=7-00:00:00
#SBATCH --job-name=pcodGMAC2_wgp-pca-admix
#SBATCH --output=/home/sschaal/pcod/20230621/job_outfiles/pcodGMAC2_wholegenome_polymorphic_%A.out

module unload bio/pcangsd/0.99
module load bio/pcangsd/0.99
source /opt/bioinformatics/venv/pcangsd-0.99/bin/activate

pcangsd.py -threads 10 -beagle /home/sschaal/pcod/20230621/gls/pcodGMAC2_wholegenome_polymorphic.beagle.gz -o /home/sschaal/pcod/20230621/pca/pcodGMAC2_wholegenome-polymorphic -sites_save