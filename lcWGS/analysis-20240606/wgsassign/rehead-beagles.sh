#!/bin/bash

#SBATCH --job-name=rehead-beagles
#SBATCH --output=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/wgsassign/err-out/relabel-beagles.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 8
#SBATCH -t 2-0:0:0

## Relabel beagle columns script
# Laura Spencer, laura.spencer@noaa.gov or lhs3@uw.edu
# This script relabels columns in a beagle file based on a bamlist.txt that was used to generate the beagle file (in angsd)

## ---  STEP 0: SET DIRECTORIES AND INPUT BEAGLE VARIABLES
base=/home/lspencer/pcod-lcwgs-2023/analysis-20240606

# beagle1 information.  beagle1 should be derived from a set of individuals with more confidence in major/minor allele designations, e.g. reference population
beagle1=pcod-refs_wholegenome_wgassign_filtered
bamlist1=${base}/reference/pcod-refs_filtered_bamslist.txt #this file is what you feed into angsd; it contains just a list of paths to bam files
prefix1=ABLG # This code assumes bamlists have full paths and the sample IDs are in the format "/path/to/bams/string_PREFIX####_string.bam", and extracts "PREFIX####" to label beagle columns

# beagle 2 information.
beagle2=pcod-exp_wholegenome_wgassign
bamlist2=${base}/experimental/pcod-exp_filtered_bamslist.txt
prefix2=GM

rm -r ${base}/wgsassign/rehead-temp
mkdir ${base}/wgsassign/rehead-temp
temp=${base}/wgsassign/rehead-temp

# Define arrays for the input and output file sets
bamlists=("${bamlist1}" "${bamlist2}")
beagle_names=("${temp}/beagle1.names.txt" "${temp}/beagle2.names.txt")
prefixes=("${prefix1}" "${prefix2}")
beagles=("${base}/reference/gls_wgassign/${beagle1}.beagle.gz" "${base}/experimental/gls_wgassign/${beagle2}.beagle.gz")
output_files=("${base}/wgsassign/${beagle1}_rehead.beagle.gz" "${base}/wgsassign/${beagle2}_rehead_beagle.gz")

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

rm -r ${temp}
 
done

