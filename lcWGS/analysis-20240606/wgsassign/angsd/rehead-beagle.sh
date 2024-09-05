#!/bin/bash

#SBATCH --job-name=rehead-beagles
#SBATCH --output=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/wgsassign/err-out/rehead-one-beagle.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 8
#SBATCH -t 2-0:0:0

## Relabel beagle columns script
# Laura Spencer, laura.spencer@noaa.gov or lhs3@uw.edu
# This script relabels columns in a beagle file based on a bamlist.txt that was used to generate the beagle file (in angsd)

base=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/wgsassign/angsd

beagle=${base}/pcod_wholegenome_wgassign.beagle.gz
bamlist=${base}/refs-exp_filtered-bamslist.txt #this file is what you feed into angsd; it contains just a list of paths to bam files
prefix1="ABLG" # This code assumes bamlists have full paths and the sample IDs are in the format "/path/to/bams/string_PREFIX####_string.bam", and extracts "PREFIX####" to label beagle columns
prefix2="GM"
beagle_names_file=${base}/beagle.column.names.txt
output_file=${base}/pcod_wholegenome_wgsassign_rehead.beagle.gz

# Create a temp file that will store new column names
rm -f $beagle_names_file # remove in case it already exists
cat >> $beagle_names_file << END
marker
allele1
allele2
END

# Add sample IDs repeated 3x with the addition of "_AA", "_AB", "_BB" (major homo., hetero., minor homo.)
sed -E "s/.*(_(${prefix1}|${prefix2})[0-9]+)_.*$/\1/" $bamlist | \
sed 's/_//g' | \
awk 'BEGIN {suffixes[1] = "_AA"; suffixes[2] = "_AB"; suffixes[3] = "_BB"} \
{for (i = 1; i <= 3; i++) print $0 suffixes[i]}' >> $beagle_names_file

# Read in new column names
header=$(paste -sd'\t' $beagle_names_file)

# Create new beagle file with renamed columns
{ echo -e "$header"; zcat $beagle | tail -n +2; } | gzip > $output_file

#rm -r ${temp}

