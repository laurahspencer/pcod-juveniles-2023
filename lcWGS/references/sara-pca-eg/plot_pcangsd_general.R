#######################################
### SCRIPT FOR PLOTTING PCOD PCAS
### Sara Michele Schaal
### Feburary 10, 2022
######################################


######################################
### INSTALL PACKAGES & LOAD FUNCTIONS

packages_needed <- c("ggplot2", "plotly", "ggpubr", "tidyverse", "plyr")

for(i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
  library(packages_needed[i], character.only = TRUE)
}

######################################


######################################
## DIRECTORIES AND FILE NAMES

DATADIR <- "C:/Users/sara.schaal/Work/Pacific_Cod/data/"
BAM_LIST <- "pcod20230103_filtered_combo_bamslist.txt"
SAMPLEMETADATA <- "20220518_PCOD_Database.csv" 
#METADATA <- "fst_meta_data.txt"

######################################

######################################
## Read in Data

bam_df <- read.delim(paste0(DATADIR,BAM_LIST), header = F, col.names = "sample_id")
seq_samples <- read.delim(paste0(DATADIR,SAMPLEMETADATA), header = TRUE, sep = ",")

######################################

######################################
## Manipulate Data for Plotting

head(bam_df)

dim(bam_df)
sampleIDS <- NULL
for(i in 1:nrow(bam_df)){
 sampleIDS <-  rbind(sampleIDS, str_extract(bam_df$sample_id[i], "ABLG\\d+"))
}

sampleIDs_df <- as.data.frame(sampleIDS)
sampleNUMs <- NULL
for(i in 1:nrow(sampleIDs_df)){
  sampleNUMs <-  rbind(sampleNUMs, str_extract(sampleIDs_df$V1[i], "\\d+"))
}

sampleNUMs_df <- as.data.frame(sampleNUMs)
sampleNUMs_df$V1 <- as.numeric(sampleNUMs_df$V1)

seq_samples <- read.delim(paste0(DATADIR,SAMPLEMETADATA), header = TRUE, sep = ",")
head(seq_samples)
str(seq_samples)
plotData <- left_join(sampleNUMs_df, seq_samples, by = c("V1" = "ABLG"))[c("V1", "Locality", "Region")]
dim(plotData)
plotData <- plotData[complete.cases(plotData),]

pop.factor.levels <- c("AK Knight-tagged", "Aleutian-tagged", "Vesteraalen-tagged",
                       "Japan", "Korea","Russia", "Zhemchug", "Pervenets Canyon", "Pribilof",
                       "Near Islands", "Tanaga Island", "Amchitka Pass", "Adak",
                       "Unimak", "Shumagins", "West_Kodiak", "Kodiak", 
                       "Cook Inlet", "PWS", "Port_Gravina", "Baikof_Bay", "Lynn_Canal",
                       "Hecate Strait")


 mypalette <- c( "orange", "goldenrod", "yellow",  # tagged - orange
                 "tomato3", "red", "navyblue",  # russia (Japan, Korea)
                "mediumblue", "cornflowerblue", "lightblue", # Bering Sea - blues (Zhemchug, Pribilof)
                "#006d2c", "#31a354", "#74c476", "lightgreen",# AI - greens 
                "antiquewhite", "mistyrose", "#fcc5c0", "palevioletred2",  # AK Pen - pinks (west kodiak)
                "#dd3497", "#ae017e", "firebrick", "darkred", "#7a0177", "purple") # GOA - purples (Port_Gravina, Baikof_Bay)


region.factor.levels <- c("Tagged", "Bering Sea", "Aleutian", "wGOA", "eGOA", "Eastern Pacific")

######################################


######################################
### Genome-wide PCA

pca <- as.matrix(read.table("./data/pcangsd/covMatrices_pcod20230117/pcod20230117_combo_wholegenome-polymorphic.cov"))
pca_e <- eigen(pca)
prop_var <-  round(pca_e$values/sum(pca_e$values)*100,2)
first3pcs <- data.frame(pca_e$vectors[,1:3])
dim(pca)
plotData$Batch <- ifelse(plotData$V1 > 10000, "Batch 2", "Batch 1")
pca_df <- cbind(plotData[c("Locality", "Region")], first3pcs)

colnames(pca_df) <- c( "pop", "region", 
                      "PC1", "PC2", "PC3")

pca_df$pop <- as.factor(pca_df$pop)
pca_df$pop <- factor(pca_df$pop, levels = pop.factor.levels)
pca_df$region <- factor(pca_df$region, levels = region.factor.levels)
pca_df$region <- revalue(pca_df$region, c("Eastern Pacific" = "Western Pacific"))

mytagalpha <- c(0.4,0.4,0.4,0.4,1)
mytagsize <- c(4,4,4,4,6)
myshapes <- c(9, 18, 15, 19, 17, 8)
genomewidePlot <- ggplot(data = pca_df, aes(x = PC1, y = PC2, color = pop, shape = region)) + 
  geom_point(aes(x = PC1 , y = PC2, color = pop, shape = region), size = 3) +
  ggtitle(paste0("Pacific cod aligned to pacific cod genome\nwhole genome PCA ")) +
  scale_color_manual(values = mypalette) +
  scale_shape_manual(values = myshapes) +
  xlab(paste0("PC1 - ", format(round(prop_var[1], 2), nsmall = 2), "%")) +
  ylab(paste0("PC2 - ", format(round(prop_var[2], 2), nsmall = 2), "%")) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

genomewidePlot

myplot.legend <- g_legend(chromPlot)
genomewidePlot.noleg <- genomewidePlot + theme(legend.position = "none")

pdf(paste0("./figures/pcas/20230103/pcod20230103_genome-wide_PCA.pdf"), width = 15, height = 7 )
print(genomewidePlot)
dev.off()

jpeg(paste0("./figures/pcas/gtSEQ_panel/pcodGMAC_postThirdDropFromGtscore_final_panel.jpg"), width = 15, height = 7, res = 150, units = "in")
print(genomewidePlot)
dev.off()


genomewidePlot <- ggplot(data = pca_df, aes(x = PC1, y = PC3, color = pop, shape = region)) + 
  geom_point(aes(x = PC1 , y = PC3, color = pop, shape = region), size = 3) +
  ggtitle(paste0("Pacific cod Genome-wide PCA PC1 vs PC3")) +
  # scale_alpha_manual(values = mytagalpha) + 
  scale_color_manual(values = mypalette) +
  scale_shape_manual(values = myshapes) +
  #  scale_size_manual(values = mytagsize) +
  xlab(paste0("PC1 - ", format(round(prop_var[1], 2), nsmall = 2), "%")) +
  ylab(paste0("PC3 - ", format(round(prop_var[3], 2), nsmall = 2), "%")) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

genomewidePlot

######################################

