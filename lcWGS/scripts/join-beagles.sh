#!/bin/bash

#SBATCH --job-name=join-beagles
#SBATCH --output=/share/afsc/temporary/join-beagles-testing/join-beagles-testing-output.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL

## Join beagles script
# Laura Spencer, laura.spencer@noaa.gov or lhs3@uw.edu
# This script joins two beagle.gz files at markers that are common in the two files, and swaps maj/min homozygous allele columns where they don't match.

## ---  STEP 0: SET DIRECTORIES AND INPUT BEAGLE VARIABLES
base=/share/afsc/temporary/join-beagles-testing
mkdir ${base}/join-beagles-temp
temp=${base}/join-beagles-temp

# beagle1 information.  beagle1 should be derived from a set of individuals with more confidence in major/minor allele designations, e.g. reference population
beagle1=${base}/beagle1.gz
bamlist1=${base}/bamlist1.txt #this file is what you feed into angsd; it contains just a list of paths to bam files
prefix1=ABLG # This code assumes bamlists have full paths and the sample IDs are in the format "/path/to/bams/string_PREFIX####_string.bam", and extracts "PREFIX####" to label beagle columns

# beagle 2 information.
beagle2=${base}/beagle2.gz
bamlist2=${base}/bamlist2.txt
prefix2=GM

# output and other settings
merged_beagle_name=testing-merged
remove_temp_files="No" # if = "Yes" temporary files will be removed. Set to "No" (or anything else) to not remove them
remove_rehead_beagles="No" #if "Yes" the intermediate beagle files that contain all original GLs but include sample IDs will be removed. Set to "No" if you don't want this file.


## --- STEP 1: RENAME COLUMNS IN BEAGLE FILES BASED ON ORDER IN BAMLIST.TXT FILES
# This is done in a for loop over both beagle files.
# NOTE: this creates temporary output beagle files with renamed columns that COULD be copied to a new location if desired!

# Define arrays for the input and output file sets
bamlists=("${bamlist1}" "${bamlist2}")
beagle_names=("${temp}/beagle1.names.txt" "${temp}/beagle2.names.txt")
prefixes=("${prefix1}" "${prefix2}")
beagles=("${beagle1}" "${beagle2}")
output_files=("${temp}/rehead_beagle1.gz" "${temp}/rehead_beagle2.gz")

# Iterate over the indices of the arrays
for i in ${!bamlists[@]}; do
    # Define variables for the current iteration
    bamlist=${bamlists[$i]}
    beagle_names_file=${beagle_names[$i]}
    prefix=${prefixes[$i]}
    beagle=${beagles[$i]}
    output_file=${output_files[$i]}

    # Create a temp file that will store new column names
    rm -f $beagle_names_file # remove in case it already exists
    cat >> $beagle_names_file << END
marker
allele1
allele2
END

    # Add sample IDs repeated 3x with the addition of "_AA", "_AB", "_BB" (major homo., hetero., minor homo.)
    sed -E "s/.*(_${prefix}[0-9]+)_.*$/\1/" $bamlist | \
    sed 's/_//g' | \
    awk -v prefix="$prefix" 'BEGIN {suffixes[1] = "_AA"; suffixes[2] = "_AB"; suffixes[3] = "_BB"} \
    {for (i = 1; i <= 3; i++) print $0 suffixes[i]}' >> $beagle_names_file

    # Read in new column names
    header=$(paste -sd'\t' $beagle_names_file)

    # Create new beagle file with renamed columns
    { echo -e "$header"; zcat $beagle | tail -n +2; } | gzip > $output_file
done

## --- STEP 2: FIND COMMON MARKERS IN THE TWO BEAGLES
zcat $beagle1 | cut -f1 | sort | uniq > ${temp}/col1_beagle1.txt
zcat $beagle2 | cut -f1 | sort | uniq > ${temp}/col1_beagle2.txt
comm -12 ${temp}/col1_beagle1.txt ${temp}/col1_beagle2.txt > ${temp}/common_markers.txt

echo "Number of common markers:"; cat ${temp}/common_markers.txt | wc -l

#clean up
[ "$remove_temp_files" = "Yes" ] && { rm ${temp}/col1_beagle1.txt ${temp}/col1_beagle2.txt; }

## --- STEP 3: FILTER BEAGLES FOR COMMON MARKERS THEN SORT BY MARKER NAME; RENAME ALLELE1 AND ALLELE2 COLUMNS TO INCLUDE "_1" AND "_2" FOR BEAGLE1 AND BEAGLE2, RESPECTIVELY

# BEAGLE1
zcat ${temp}/rehead_beagle1.gz | awk 'NR==FNR {a[$1]; next} $1 in a' ${temp}/common_markers.txt - | sort -k1,1 | sed -e 's/allele1/allele1_1/g' -e 's/allele2/allele2_1/g' > ${temp}/sorted_beagle1_filtered.beagle

# BEAGLE2
zcat ${temp}/rehead_beagle2.gz | awk 'NR==FNR {a[$1]; next} $1 in a' ${temp}/common_markers.txt - | sort -k1,1 | sed -e 's/allele1/allele1_2/g' -e 's/allele2/allele2_2/g'> ${temp}/sorted_beagle2_filtered.beagle


# Double check that markers in sorted filtered beagle files match (columns 1 are identical)

# Extract column 1 from each file and sort them (since comm requires sorted input)
cut -f1 "${temp}/sorted_beagle1_filtered.beagle" > "${temp}/beagle1_col1.txt"
cut -f1 "${temp}/sorted_beagle2_filtered.beagle" > "${temp}/beagle2_col1.txt"

# Compare the two sorted files
if comm -3 "${temp}/beagle1_col1.txt" "${temp}/beagle2_col1.txt" | grep -q '.'; then
    echo "Mismatch found between the two files."
    exit 1
else
    echo "All markers in sorted filtered beagle files match."
fi

# cleanup
rm ${temp}/beagle1_col1.txt ${temp}/beagle2_col1.txt

## --- STEP 4: LEFT-JOIN BEAGLE2 TO BEAGLE1 AND FILTER FOR ONLY MARKERS THAT HAVE THE SAME NAME, MAJOR ALLELE, AND MINOR ALLELE

# Perform the Left Join on Column 1 from each beagle file
join --check-order --header -t $'\t' -1 1 -2 1 -o auto ${temp}/sorted_beagle1_filtered.beagle ${temp}/sorted_beagle2_filtered.beagle > ${temp}/matched_rows.beagle

# Create filtered beagle file with only rows where major and minor alleles match. This should contain MOST of our markers.

rm ${temp}/exact_matches.beagle # remove file if already exists, don't want to add to it!

# Create the new file with header info
head -n 1 ${temp}/matched_rows.beagle > ${temp}/exact_matches.beagle

# Include rows where allele1 and allele2 are exact matches
awk -F'\t' 'BEGIN {col1_index = 0; col2_index = 0; col3_index = 0; col4_index = 0} \
NR == 1 {for (i = 1; i <= NF; i++) {
  if ($i == "allele1_1") col1_index = i; \
  if ($i == "allele1_2") col2_index = i; \
  if ($i == "allele2_1") col3_index = i; \
  if ($i == "allele2_2") col4_index = i
}; next} \
$col1_index == $col2_index && $col3_index == $col4_index {print}' \
${temp}/matched_rows.beagle >> ${temp}/exact_matches.beagle

echo "Number of markers with matching major / minor alleles"
cat ${temp}/exact_matches.beagle | wc -l

# clean up
[ "$remove_temp_files" = "Yes" ] && { rm ${temp}/matched_rows.beagle; }

### --- STEP 5: IDENTIFY MARKERS WHERE MAJ/MIN ALLELES ARE SWAPPED
# Specifically, here I identify markers (columns 1 ) where the major allele (column 2) in one beagle is the minor allele (column 3) in the other beagle and visa versa; save those rows from beagle2 to new file

awk 'FNR==NR {a[$1,$3,$2]=$0; next} ($1,$2,$3) in a' ${temp}/sorted_beagle1_filtered.beagle ${temp}/sorted_beagle2_filtered.beagle | sort -k1,1 > ${temp}/unmatched_rows2.beagle

# clean up
[ "$remove_temp_files" = "Yes" ] && { rm ${temp}/sorted_beagle2_filtered.beagle; }

### --- STEP 6: SWAP MAJ/MIN ALLELES AND GENOTYPE LIKELIHOODS IN THE unmatched_rows2.beagle FILE
# Specifically, first we swap columns 2<->3, then for columns 4+ we swap the homozygous columns for each sample (e.g. 4<->6, 7<->9, etc.)

awk '{tmp = $2; $2 = $3; $3 = tmp; \
for (i = 4; i <= NF; i += 3) {if ((i + 2) <= NF) {tmp = $i; $i = $(i + 2); $(i + 2) = tmp}}; \
print}' OFS='\t' ${temp}/unmatched_rows2.beagle > ${temp}/reordered_unmatched_rows2.beagle

### --- STEP 7: LEFT-JOIN BEAGLE1 AND BEAGLE2 GLS FOR MARKERS WITH SWAPPED ALLELES

# First, filter sorted beagle1 to only include markers that did not match, which we'll join with our reordered beagle2 file
cat ${temp}/unmatched_rows2.beagle | cut -f1 > ${temp}/unmatched-markers.txt
cat ${temp}/sorted_beagle1_filtered.beagle | awk 'NR==FNR {a[$1]; next} $1 in a' ${temp}/unmatched-markers.txt - | sort -k1,1 > ${temp}/unmatched_rows1.beagle

# clean up
[ "$remove_temp_files" = "Yes" ] && { rm ${temp}/sorted_beagle1_filtered.beagle; }

# Perform the Left Join on Columns 1 in beagle files with unmached rows
join --check-order --header -t $'\t' -1 1 -2 1 -o auto ${temp}/unmatched_rows1.beagle ${temp}/reordered_unmatched_rows2.beagle > ${temp}/unmatched_rows.beagle

echo "Number of markers with swapped major / minor alleles"
cat ${temp}/unmatched_rows.beagle | wc -l

# clean up
[ "$remove_temp_files" = "Yes" ] && { rm ${temp}/unmatched_rows1.beagle ${temp}/unmatched_rows2.beagle ${temp}/unmatched-markers.txt ${temp}/reordered_unmatched_rows2.beagle; }


## --- STEP 8: COMBINE GLS FOR MARKERS WITH MATCHED MAJ/MIN WITH GLS FOR MARKERS WITH SWAPPED MAJ/MIN

# Concatenate the exact_matches.beagle (where maj/min alleles matched) with the unmatched_rows.beagle (where maj/min did not agree).
{ cat ${temp}/exact_matches.beagle; cat ${temp}/unmatched_rows.beagle; } > ${temp}/merged-temp1.beagle

# sort concatenated file based on marker name (column 1)
{ cat "${temp}/merged-temp1.beagle" | head -n 1; cat "${temp}/merged-temp1.beagle" | tail -n +2 | sort -k1,1; } > "${temp}/merged-temp2.beagle"

# clean up
cat ${temp}/exact_matches.beagle | cut -f1 > ${temp}/markers_matched-alleles.txt
cat ${temp}/unmatched_rows.beagle | cut -f1 > ${temp}/markers_swapped-alleles.txt
[ "$remove_temp_files" = "Yes" ] && { rm ${temp}/exact_matches.beagle ${temp}/unmatched_rows.beagle ${temp}/merged-temp1.beagle; }


## -- STEP 9: CHECK TO MAKE SURE ALLELE1_1 == ALLELE1_2 AND ALLELE2_2 == ALLELE2_2

awk -F'\t' 'BEGIN {col1_index = 0; col2_index = 0; col3_index = 0; col4_index = 0; all_match = 1} \
NR == 1 {for (i = 1; i <= NF; i++) {if ($i == "allele1_1") col1_index = i; \
if ($i == "allele1_2") col2_index = i; if ($i == "allele2_1") col3_index = i; \
if ($i == "allele2_2") col4_index = i} next} \
{if ($col1_index != $col2_index || $col3_index != $col4_index) all_match = 0} \
END {if (all_match) print "All major and minor alleles match, proceed to next step!"; \
else print "ERROR: There are discrepancies between the alleles."}' ${temp}/merged-temp2.beagle

## STOP - DID YOU GET AN ERROR MESSAGE IN THE LAST STEP?

## -- STEP 10: IF STEP 9 LOOKED GOOD, REMOVE DUPLICATED COLUMNS
# Finally I need to remove "allele1_2" and "allele2_2" columns

# identify column index for columns "allele1_2" and "allele2_2"
index1=$(head -n 1 ${temp}/merged-temp2.beagle | tr '\t' '\n' | grep -n "allele1_2" | cut -d':' -f1)
index2=$(head -n 1 ${temp}/merged-temp2.beagle | tr '\t' '\n' | grep -n "allele2_2" | cut -d':' -f1)

# remove those two columns
cut --complement -f$index1,$index2 ${temp}/merged-temp2.beagle | gzip > ${base}/${merged_beagle_name}.beagle.gz

echo "Total number of sites in resulting joined beagle"
zcat ${base}/${merged_beagle_name}.beagle.gz | wc -l

# Report number of markers that were not included in output 
# add this code 

# clean up
[ "$remove_temp_files" = "Yes" ] && { rm ${temp}/merged-temp2.beagle; }

## -- STEP 11: Remove beagle files with headers renamed to include sample IDs
[ "$remove_rehead_beagles" = "Yes" ] && { rm ${temp}/rehead_beagle1.gz ${temp}/rehead_beagle2.gz; }

