#!/bin/bash

#SBATCH --job-name=wgsassign
#SBATCH --output=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/wgsassign/err-out/region-assign-using-one-beagle-top25k.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 20
#SBATCH -t 2-0:0:0

source ~/.bashrc
mamba activate WGSassign

# Define input files
base="/home/lspencer/pcod-lcwgs-2023/analysis-20240606/wgsassign"
#input_beagle="${base}/refs-exp-merged.beagle.gz" # combined beagle that I created previously
input_beagle="${base}/angsd/pcod_wholegenome_wgsassign_rehead.beagle.gz" #beagle created with all ref and exp samples fed into angsd at same time

# Create file that lists all experimental IDs; they are pulled from the combined beagle file and listed in order that they appear in that file 
zcat ${input_beagle} | head -n 1 | cut --complement -f1-3 | tr '\t' '\n' | grep "GM" | sed -e 's/_AA//g' -e 's/_AB//g' -e 's/_BB//g' | uniq > ${base}/assignment/exp-IDs.txt

id_file="${base}/assignment/exp-IDs.txt"
exp_beagle_all="${base}/assignment/exp-all-sites_using-one-beag.beagle.gz"
best_snps_beagle="${base}/snp-testing/by-marine-region/using-one-beagle/testing-top-25000.beagle.gz"
best_snps_afs="${base}/snp-testing/by-marine-region/using-one-beagle/testing-LOOs/top-25000.pop_af.npy"
exp_beagle_best="${base}/assignment/pcod-exp_top-25k-region_using-one-beag.beagle"
outname="${base}/assignment/pcod-experimental-assign-region_top-25k-snps_using-one-beag"

# First I need to subset our combined beagle file to get only our experimental samples. This is to make sure I'm using the exact same sites and GLs are correctly assigned to major/minor alleles

# Extract IDs from the id_file into a regex pattern (exact match)
ids=$(awk '{print "^" $1 "$"}' $id_file | paste -sd '|' -)

# Process the input file
zcat $input_beagle | awk -v ids="$ids" '
BEGIN { FS=OFS="\t"; }
NR==1 {
    # Print the first three columns (marker, allele1, allele2)
    header_line = $1 OFS $2 OFS $3;
    count = 0;
    for (i=4; i<=NF; i++) {
        # Remove suffixes and check if the column matches any ID
        col_name = gensub(/_AA$|_AB$|_BB$/, "", "g", $i)
        if (col_name ~ ids) {
            header_line = header_line OFS $i
            cols_to_print[count++] = i
        }
    }
    print header_line
}
NR>1 {
    line = $1 OFS $2 OFS $3;
    for (i=0; i<count; i++) {
        line = line OFS $cols_to_print[i]
    }
    print line
}
' | gzip > $exp_beagle_all


# Filter newly created experimental beagle file for sites we identified as predicting assignment 

# Create sites file containing the best SNPs to use for population assignment
zcat $best_snps_beagle | cut -f1 > ${base}/assignment/assign-region-25k_using-one-beag.sites

# Filter beagle with top N best sites
# combine header row with filtered sites
zcat ${exp_beagle_all} | head -n 1 > ${exp_beagle_best}
awk 'NR==FNR{c[$1]++;next};c[$1]' ${base}/assignment/assign-region-25k_using-one-beag.sites <(zcat ${exp_beagle_all}) >> ${exp_beagle_best}
gzip ${exp_beagle_best}

# Generate file listing samples in order they appear in beagle file using sample IDs in header row
zcat ${exp_beagle_best}.gz | cut --complement -f1-3 | head -n 1 | tr '\t' '\n' | \
sed -e 's/_AA//g' -e 's/_AB//g' -e 's/_BB//g' | uniq > ${base}/assignment/pcod-exp-sample-order-region_using-one-beag.txt

# Run assignment to identify source population of experimental fish 
WGSassign --beagle ${exp_beagle_best}.gz --pop_af_file ${best_snps_afs} --get_pop_like --out ${outname} --threads 20
