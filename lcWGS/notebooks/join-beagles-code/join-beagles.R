### Create one beagle file with both reference fish AND experimental fish
#I want to use PCA to look for overlaps among our reference fish and experimental fish. So, I need to merge the two beagle files. 

### Here are the final 2 beagle files that contain genotype likelihoods for the same sites.  
#Reference fish:    pcod-refs_wholegenome_wgassign_filtered.beagle.gz  
#Experimental fish: pcod-exp_wholegenome_wgassign.beagle.gz  

# The beagle giles are gzipped, unzip them (NOTE: This removes thes .gz files too!, there's probably an option to keep them)
library(R.utils)
gunzip("pcod-exp_wholegenome_wgassign.beagle.gz")
gunzip("pcod-refs_wholegenome_wgassign_filtered.beagle.gz")

sample.order.exp.beagle <- paste(rep(sample.order.exp, each=3), rep(1:3, times=length(sample.order.exp)), sep="_")
sample.order.refs.beagle <- paste(rep(sample.order.refs, each=3), rep(1:3, times=length(sample.order.refs)), sep="_")

gls.exp <- read_delim("pcod-exp_wholegenome_wgassign.beagle", skip = 1, col_names = c("marker", "allele1", "allele2", sample.order.exp.beagle))
gls.exp[1:10, 1:10] #preview

gls.ref <- read_delim("pcod-refs_wholegenome_wgassign_filtered.beagle", skip = 1, col_names = c("marker", "allele1", "allele2", sample.order.refs.beagle))
gls.ref[1:10, 1:10] #preview

# Are the markers in thes same order and do they exactly match? 
all(gls.exp$marker == gls.ref$marker)

# Join the two beagle dataframes 
gls.all <- left_join(gls.ref, gls.exp, by=c("marker", "allele1", "allele2")) 

# Correct number of sample columns? 
gls.all[-c(1:3)] %>% ncol() == length(sample.order.refs)*3 + length(sample.order.exp)*3

# BUT ISSUE - there are 4,157 SNPs where Allele 1 and Allele 2 are swapped in our reference fish and experimental fish, so they aren't joining correctly (creating NAs)
gls.nas <- gls.all %>% filter(if_any(everything(), is.na))

# Create dataframe with subset of sites were alleles are "corrected", i.e. reference/alternative alleles are swapped in the experimental fish 
gls.allele.swap <- gls.nas[c("marker", "allele1", "allele2")] %>% left_join(gls.ref) %>%
  left_join(
    gls.exp %>% filter(marker %in% c(gls.nas$marker)) %>%
      dplyr::select(marker, paste(rep(sample.order.exp, each=3), rep(3:1, times=length(sample.order.exp)), sep="_")),  # For these markers, swap reference / alternate genotypes in experimental fish
    by="marker")

gls.all.corrected <- rbind(gls.all %>% filter(marker %!in% gls.nas$marker), 
                           gls.allele.swap %>%  `colnames<-`(c("marker", "allele1", "allele2", sample.order.refs.beagle, sample.order.exp.beagle))) %>%
  arrange(match(marker, gls.all$marker)) # reorder markers by chrom_site 


# Manually review GLs for reference homozygous alleles ("GM###_1") and alternate homozygous alleles ("GM###_3")

# Here's the original experimental allele frequencies 
(gls.exp %>% filter(marker %in% c(gls.nas$marker))) %>% head(n=20) %>%
  dplyr::select(marker, last_col(offset = 12 - 1):last_col()) #look at last columns

# Here's the corrected / swapped experimental allele frequencies
gls.all.corrected %>% filter(marker %in% gls.nas$marker) %>% head(n=20) %>%
  dplyr::select(marker, last_col(offset = 12 - 1):last_col()) #look at last columns to see 

# Any NA in corrected dataframe? 
gls.all.corrected %>% filter(if_any(everything(), is.na)) %>% nrow()

# Write beagle to file, then gzip 
write_delim(gls.all.corrected, file = "refs-and-exp.beagle", delim = "\t")
gzip("refs-and-exp.beagle")