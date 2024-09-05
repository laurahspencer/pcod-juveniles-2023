# I executed this R script on Sedna using the sbatch script "dexseq.sh"
# To run interactively using high memory node use `srun --pty -p himem -c 24 /bin/bash`

# Add all required libraries that are installed with install.packages() here
list.of.packages <- c("tidyverse", "plotly", "janitor")
# Add all libraries that are installed using BiocManager here
bioconductor.packages <- c("DEXSeq", "BiocParallel")

# This commented out code need only be run once per machine (installs packages). Don't re-do it, it can take a while.
# Install BiocManager if needed
#if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

## Get names of all required packages that aren't installed
#new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
#new.bioc.packages <- bioconductor.packages[!(bioconductor.packages %in% installed.packages()[, "Package"])]
## Install all new packages
#if(length(new.packages)) install.packages(new.packages)
#if(length(new.bioc.packages)) BiocManager::install(new.bioc.packages)

# Load all required libraries
all.packages <- c(list.of.packages, bioconductor.packages)
lapply(all.packages, FUN = function(X) {
  do.call("require", list(X))
})

source("load_SubreadOutput.R")
`%!in%` = Negate(`%in%`)

# Load in datasets
b <- c(data.frame(read.table("featurecounts_exon_dexseq_filtered", header = T, stringsAsFactors = F, fill = FALSE)) %>%
  dplyr::select(starts_with("sample_")) %>% colnames())

samp <- b %>% as.data.frame() %>%
  set_names("sample_name") %>% mutate(sample_number=gsub("sample_", "", sample_name)) %>%
  left_join(
    read_delim("DESeq2_Sample_Information.txt", delim="\t") %>%
    clean_names() %>%
    mutate(condition=case_when(
      temp_treatment=="9" ~ "a_9",
      temp_treatment=="0" ~ "b_0",
      temp_treatment=="5" ~ "c_5",
      temp_treatment=="16" ~ "d_16")) %>%
    mutate(condition=factor(as.character(condition)),
           sample_number=as.character(sample_number)) %>%
    mutate(sample_name=gsub("RESUB-", "", sample_name))) %>%
  mutate(sample_number=gsub(".G|.S", "", sample_number))%>%
  group_by(sample_number) %>% fill(c(tank, condition), .direction = "down") %>% #fill in missing tank and treatment info
  mutate(tissue_type=case_when(
    grepl(".G", sample_name) ~ "Gill",
    grepl(".S", sample_name) ~ "Spleen",
    TRUE ~ "Liver") %>% as.factor()) %>%
  column_to_rownames("sample_name")

samp.0.9 <- samp %>% filter(temp_treatment %in% c(0,9)) %>% droplevels()
samp.5.9 <- samp %>% filter(temp_treatment %in% c(5,9)) %>% droplevels()
samp.16.9 <- samp %>% filter(temp_treatment %in% c(16,9)) %>% droplevels()

save(samp.0.9, file="samp.0.9")
save(samp.5.9, file="samp.5.9")
save(samp.16.9, file="samp.16.9")

# Create DEXSeq object
dxd.fc <- DEXSeqDataSetFromFeatureCounts(countfile = "featurecounts_exon_dexseq_filtered",
                                         flattenedfile = "/home/lspencer/references/pcod-ncbi/GCF_031168955.1_ASM3116895v1_genomic_flat.gtf",
                                         sampleData = samp %>% dplyr::select(condition))

# Normalization
dxd.fc = estimateSizeFactors( dxd.fc )

# Dispersion estimation
## This is very memory intensive. Here I use multiple cores.
BPPARAM = BiocParallel::MulticoreParam(24)
dxd.fc = DEXSeq::estimateDispersions( dxd.fc, BPPARAM=BPPARAM )

pdf(file = "disp-estimates.pdf", width = 12, height = 9);
plotDispEsts( dxd.fc )
dev.off()

dxd.fc = testForDEU( dxd.fc, BPPARAM=BPPARAM )

# Estimate fold changes
dxd.fc = estimateExonFoldChanges( dxd.fc, fitExpToVar="condition", BPPARAM=BPPARAM )
save(dxd.fc, file="dxd.fc.final")

# Extract results object
dxr1 = DEXSeqResults( dxd.fc )
dxr1
save(dxr1, file="dxr1")

# Explore results
table ( dxr1$padj < 0.1 )
table ( tapply( dxr1$padj < 0.1, dxr1$groupID, any ) )

pdf(file = "volcano-plot.pdf", width = 12, height = 9);
plotMA( dxr1, cex=0.8 )
dev.off()

# Rerun DEXSeq on subsetted objects, one for each treatment (0, 5, and 16, all vs. 9)

dxd.fc <- DEXSeqDataSetFromFeatureCounts(countfile = "featurecounts_exon_dexseq_filtered",
                                         flattenedfile = "/home/lspencer/references/pcod-ncbi/GCF_031168955.1_ASM3116895v1_genomic_flat.gtf",
                                         sampleData = samp %>% dplyr::select(condition))


dxd.fc.0vs9 <- DEXSeqDataSetFromFeatureCounts(countfile = "featurecounts_exon_dexseq_filtered_0.9",
                                         flattenedfile = "/home/lspencer/references/pcod-ncbi/GCF_031168955.1_ASM3116895v1_genomic_flat.gtf",
                                         sampleData = samp.0.9 %>% dplyr::select(condition))

dxd.fc.5vs9 <- DEXSeqDataSetFromFeatureCounts(countfile = "featurecounts_exon_dexseq_filtered_5.9",
                                         flattenedfile = "/home/lspencer/references/pcod-ncbi/GCF_031168955.1_ASM3116895v1_genomic_flat.gtf",
                                         sampleData = samp.5.9 %>% dplyr::select(condition))

dxd.fc.16vs9 <- DEXSeqDataSetFromFeatureCounts(countfile = "featurecounts_exon_dexseq_filtered_16.9",
                                         flattenedfile = "/home/lspencer/references/pcod-ncbi/GCF_031168955.1_ASM3116895v1_genomic_flat.gtf",
                                         sampleData = samp.16.9 %>% dplyr::select(condition))

dxd.fc.treats <- vector(mode = "list", length = 3)
names(dxd.fc.treats) <-c ("dxd.fc.0vs9", "dxd.fc.5vs9", "dxd.fc.16vs9")
dxd.fc.treats[1] <- dxd.fc.0vs9
dxd.fc.treats[2] <- dxd.fc.5vs9
dxd.fc.treats[3] <- dxd.fc.16vs9
dxr.treats <- vector(mode = "list", length = 3)
names(dxr.treats) <- c("dxr.0vs9", "dxr.5vs9", "dxr.16vs9")

for (i in 2:3){

 # Normalization
 dxd.fc.treats[[i]] = estimateSizeFactors( dxd.fc.treats[[i]] )

 # Dispersion estimation
 ## This is very memory intensive. Here I use multiple cores.
 BPPARAM = BiocParallel::MulticoreParam(24)
 dxd.fc.treats[[i]] = DEXSeq::estimateDispersions( dxd.fc.treats[[i]], BPPARAM=BPPARAM )

 pdf(file = paste("disp-estimates_", gsub("dxd.fc.", "", names(dxd.fc.treats[i])),".pdf", sep=""), width = 12, height = 9);
 plotDispEsts( dxd.fc.treats[[i]] )
 dev.off()

 dxd.fc.treats[[i]] = testForDEU( dxd.fc.treats[[i]], BPPARAM=BPPARAM )

 # Estimate fold changes
 dxd.fc.treats[[i]] = estimateExonFoldChanges( dxd.fc.treats[[i]], fitExpToVar="condition", BPPARAM=BPPARAM )

 # Extract results object
 dxr.treats[[i]] = DEXSeqResults( dxd.fc.treats[[i]] )

 pdf(file = paste("volcano-plot_", gsub("dxd.fc.", "", names(dxd.fc.treats[i])),".pdf", sep=""), width = 12, height = 9);
 plotMA(  dxr.treats[[i]], cex=0.8 )
 dev.off()

 save(dxd.fc.treats, file="dxd.fc.treats")
 save(dxr.treats, file="dxr.treats")
 }

# Save again, just in case! 
 save(dxd.fc.treats, file="dxd.fc.treats")
 save(dxr.treats, file="dxr.treats")
