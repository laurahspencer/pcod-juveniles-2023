#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=7-00:00:00
#SBATCH --job-name=angsd_4_wgassign
#SBATCH --output=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/wgsassign/err-out/angsd-refs-exp_%A-%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --array=1-24%24

## ANGSD Settings
# -setminDepth 1000 is based on 50% of the individuals and the average sequencing coverage of 3X (i.e. 1185 ~ 3x * 0.5*790)
# -setmaxDepth 5000 is based on 2 times all the individuals and the average depth, to avoid sites overly sequenced (i.e. 4740 ~ 2 * 790 * 3x)
# -minInd 395 sets calling only variants with ~50% of individuals (n=790)
# -minMaqQ was 15, now 30
# -minQ was 15, now 33
# -GL was 1, now 2
# -sites restricts variant calling to sites we previously identified in both reference AND experimental fish when angsd was run separately on each data set 


module unload bio/angsd/0.933
module load bio/angsd/0.933

base=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/wgsassign

#JOBS_FILE=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/scripts/pcod-exp_angsdARRAY_input.txt
#mapfile -t IDS <${JOBS_FILE}

JOBS_FILE=${base}/angsd/chroms.txt
IDS=$(cat ${JOBS_FILE})

for sample_line in ${IDS}
do
        job_index=$(echo ${sample_line} | awk -F ":" '{print $1}')
        chrom=$(echo ${sample_line} | awk -F ":" '{print $2}')
        if [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then
                break
        fi
done

angsd -b ${base}/angsd/refs-exp_filtered-bamslist.txt \
-ref /home/lspencer/references/pcod-ncbi/GCF_031168955.1_ASM3116895v1_genomic.fa \
-r ${chrom}: \
-out ${base}/angsd/pcod-_Chr${job_index}_wgassign \
-nThreads 10 \
-sites ${base}/angsd/refs-exp.sites \
-uniqueOnly 1 \
-remove_bads 1 \
-trim 0 \
-C 50 \
-minMapQ 30 \
-minQ 33 \
-doCounts 1 \
-setminDepth 1000 \
-setmaxDepth 5000 \
-GL 2 \
-doGlf 2 \
-doMaf 1 \
-doMajorMinor 1 \
-doDepth 1 \
-dumpCounts 3 \
-only_proper_pairs 1 \
-minInd 395 \
-minmaf 0.05 \
-SNP_pval 1e-10 \
-baq 1
