#!/bin/bash

#SBATCH --job-name=subset-beagle
#SBATCH --output=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/wgsassign/err-out/subset-using-one-beagle-region.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 8
#SBATCH -t 2-0:0:0

# Input files
base="/home/lspencer/pcod-lcwgs-2023/analysis-20240606/wgsassign"
#input_beagle="${base}/refs-exp-merged.beagle.gz" # combined beagle that I created previously 
input_beagle="${base}/angsd/pcod_wholegenome_wgsassign_rehead.beagle.gz" # one beagle I made from all ref & exp fish at same time
#cat ${base}/snp-testing/by-marine-region/test_bams-list2.txt | cut -f2,3 > ${base}/snp-testing/by-marine-region/test-IDs2.txt
id_file="${base}/snp-testing/by-marine-region/test-IDs2.txt"
output_beagle="${base}/snp-testing/by-marine-region/using-one-beagle/testing-all-sites.beagle.gz"

### First I need to subset my beagle file for only my test individuals 

# Extract IDs from the id_file into a regex pattern (exact match)
ids=$(awk '{print "^" $1 "$"}' $id_file | paste -sd '|' -)

# Select first three columns from beagle and all columns containing desired sample IDs
zcat $input_beagle | awk -v ids="$ids" '
BEGIN { FS=OFS="\t"; }
NR==1 {
    # Print the first three columns (marker, allele1, allele2)
    header_line = $1 OFS $2 OFS $3;
    for (i=4; i<=NF; i++) {
        # Remove suffixes and check if the column matches any ID
        col_name = gensub(/_AA$|_AB$|_BB$/, "", "g", $i)
        if (col_name ~ ids) {
            header_line = header_line OFS $i
            cols_to_print[i] = 1
        }
    }
    print header_line
}
NR>1 {
    line = $1 OFS $2 OFS $3;
    for (i=4; i<=NF; i++) {
        if (i in cols_to_print) {
            line = line OFS $i
        }
    }
    print line
}
' | gzip > $output_beagle


#### Now I need to filter my test-samples.beagle.gz and create 9 separate beagle files, one each containing the n sites/snps that
# I identified as predicting population:

numbers=(10 50 100 500 1000 5000 10000 25000 50000)

for n in "${numbers[@]}"
do

  # name of file to save
  sites=training.top_${n}_sites
  
  # filter sites 
  awk 'NR==FNR{c[$1]++;next};c[$1]' ${base}/snp-training/by-marine-region/${sites} <(zcat ${output_beagle}) | gzip > ${base}/snp-testing/by-marine-region/using-one-beagle/testing-top-${n}.beagle.gz
 
done

for file in ${base}/snp-testing/by-marine-region/using-one-beagle/testing-top-*.gz; do echo "$file has $(zcat "$file" | wc -l) markers"; done
