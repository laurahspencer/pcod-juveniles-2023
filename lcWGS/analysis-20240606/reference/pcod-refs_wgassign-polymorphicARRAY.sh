#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=7-00:00:00
#SBATCH --job-name=angsd_wgassign
#SBATCH --output=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/job_outfiles/pcod-refs_wgassign-polymorphic_%A-%a.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --array=1-24%24

## ANGSD Settings
# -setminDepth 950 is based on 50% of the individuals and the average sequencing coverage of 1.6X (i.e. 950 ~ 3x * 0.5*633)
# -setmaxDepth 9000 is based on 2 times all the individuals and the average depth, to avoid sites overly sequenced (i.e. 4000 ~ 2 * 633 * 3x)
# -minInd 315 sets calling only variants with ~50% of individuals
# -minMaqQ was 15, now 30
# -minQ was 15, now 33
# -GL was 1, now 2

module unload bio/angsd/0.933
module load bio/angsd/0.933

JOBS_FILE=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/scripts/pcod-refs_angsdARRAY_input.txt
mapfile -t IDS <${JOBS_FILE}

# This loop to accommodate the spaces in the genome fasta headers, but ANGSD does not, so don't actually use it
#for sample_line in "${IDS[@]}"
#do
#        job_index=$(echo "$sample_line" | awk -F ":" '{print $1}')
#        contig=$(echo "$sample_line" | awk -F ":" '{print $2}')
#        if [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then
#                break
#        fi
#done

JOBS_FILE=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/scripts/pcod-refs_angsdARRAY_input.txt
IDS=$(cat ${JOBS_FILE})

for sample_line in ${IDS}
do
        job_index=$(echo ${sample_line} | awk -F ":" '{print $1}')
        contig=$(echo ${sample_line} | awk -F ":" '{print $2}')
        if [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then
                break
        fi
done

angsd -b /home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/pcod-refs_filtered_bamslist.txt \
-ref /home/lspencer/references/pcod-ncbi/GCF_031168955.1_ASM3116895v1_genomic.fa \
-r ${contig}: \
-out /home/lspencer/pcod-lcwgs-2023/analysis-20240606/reference/gls_wgassign/pcod-refs_Chr${job_index}_wgassign \
-nThreads 10 \
-uniqueOnly 1 \
-remove_bads 1 \
-trim 0 \
-C 50 \
-minMapQ 30 \
-minQ 33 \
-doCounts 1 \
-setminDepth 950 \
-setmaxDepth 4000 \
-GL 2 \
-doGlf 2 \
-doMaf 1 \
-doMajorMinor 1 \
-doDepth 1 \
-dumpCounts 3 \
-only_proper_pairs 1 \
-minInd 315 \
-minmaf 0.05 \
-SNP_pval 1e-10 \
-baq 1
