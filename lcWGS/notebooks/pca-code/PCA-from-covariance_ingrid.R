### Explore population structure in reference fish 
# I used `pcangsd` to generate covariance matrix from genome-wide genotype likelihoods with script reference/scripts/pcod-refs_wholegenome_pcangsd.sh. 
# Before using WGSassign to identify population-of-origin for my experimental fish, I want to explore my reference fish population structure. I will do this using PCA. 

require(tidyverse)
require(plotly)
require(readxl)
require(janitor)
source("../scripts/biostats.R")
`%!in%` = Negate(`%in%`)

sample.info.ingrid <- 
  # read in experimental fish metadata 
  read_delim("lcWGS/notebooks/pca-code/allcod_filtered_bamslist_meta_rmasia.csv") %>% clean_names() %>% #read in reference fish metadata
  mutate_at(c("marine_region", "marine_region2", "location1"), factor) %>% #convert location columns to factor
  mutate(group="reference", bay=NA) %>%  #Add a "Group" columns, add Bay column with NA values
  dplyr::select(group, number, order, marine_region, marine_region2, location1, bay) %>% #select only desired columns 

  # add reference fish metadata
  rbind(
    read_delim("lcWGS/notebooks/pca-code/selected_juvPcod_lcWGS_metadata.csv") %>% clean_names() %>% #read in juvenile fish metadata
      mutate(number=paste("JUV", x2, sep=""), group="juveniles", marine_region=NA, marine_region2=NA) %>%  #Add column with "JUV" pasted with number ID, add "Group" & MarineRegion columns
      dplyr::rename("location1"="region") %>% 
      mutate_at(c("bay", "location1"), factor) %>%  # Convert to factor
      dplyr::select(group, number, order, marine_region, marine_region2, location1, bay)) 

# Vector of sample numbers in order for covariance matrix   
sample.order.refs.ingrid <- (sample.info.ingrid %>% filter(group=="reference") %>% arrange(order))$number

# Vector of sample numbers in order for covariance matrix   
sample.order.juvs.ingrid <- (sample.info.ingrid %>% filter(group=="juveniles") %>% arrange(order))$number

# Read in covariance matrix and add sample IDs based on order listed in .bam file list
# outlier=c("") #If you want to remove specific samples, add their numbers here (e.g. "ABLG2518")

genome.cov.refs.ingrid <- read_delim(file="lcWGS/notebooks/pca-code/allcod_wholegenome-polymorphic_rmasia.cov", col_names = sample.order.refs.ingrid) %>% 
  #    dplyr::select(-outlier) %>%  # Use this line to remove specific outlier samples defined above
  as.matrix() %>%
  `rownames<-`(sample.order.refs.ingrid)

# Get sample IDs for populations we'd like to include in PCAs
pops.refs.ingrid <- (sample.info.ingrid %>% 
#                      filter(location1 %!in% c("westKodiak")) %>% # If you want to remove some of the reference "location1" groupings, use this line too 
                       filter(group=="reference"))$number 

# Run PCA
pca.refs.ingrid <- prcomp(genome.cov.refs.ingrid[,pops.refs.ingrid], scale=F) #scale=F for variance-covariance matrix
#pca.eigenval(pca.princomp) #The Proporation of Variance = %variance 
pc.percent.ingrid <- pca.eigenval(pca.refs.ingrid)[2,1:6]*100 #PC % for axes 1-6
screeplot(pca.refs.ingrid, bstick=FALSE)  #inspect scree plot, which axes influential? 
pc.percent.ingrid[1:2] %>% sum() # total percent explained by PCs 1 & 2

#### Generate dataframe with prcomp results 
pc.scores.refs.ingrid <- data.frame(number = colnames(genome.cov.refs.ingrid[,pops.refs.ingrid]),
                             PC1 = pca.refs.ingrid$rotation[,1],    # the first eigenvector
                             PC2 = pca.refs.ingrid$rotation[,2],    # the second eigenvector
                             PC3 = pca.refs.ingrid$rotation[,3],    # the third eigenvector
                             PC4 = pca.refs.ingrid$rotation[,4],    # the fourth eigenvector
                             PC5 = pca.refs.ingrid$rotation[,5],    # the fourth eigenvector
                             PC6 = pca.refs.ingrid$rotation[,6],    # the fourth eigenvector
                             stringsAsFactors = FALSE)
#shapiro.test(pca.princomp$x) #Shapiro test 
#hist(pca.princomp$x) #Distribution normal? 

# Add metadata
pc.scores.refs.ingrid <- left_join(pc.scores.refs.ingrid, sample.info.ingrid) %>% droplevels() 

# Create dataframe to loop through for PC biplots 
axes <- data.frame("pc.x"=c(1,1,1,1,1,2,2,2,2,3,3,3,4,4,5), 
                   "pc.y"=c(2,3,4,5,6,3,4,5,6,4,5,6,5,6,6)) %>%
  mutate(pc.x=paste("PC", pc.x, sep=""), pc.y=paste("PC", pc.y, sep=""))

# Variance explained by each PC axis
variance.ingrid <- pc.percent.ingrid %>% as.data.frame() %>% set_names("variance") %>% rownames_to_column("axis") %>% 
  mutate(axis=as.numeric(gsub("PC", "", axis))) 

# PLOTS
pcas.genome.ingrid <- list(1:nrow(axes))
for (i in 1:nrow(axes)){
  pcas.genome.ingrid[[i]] <- 
    ggplotly( #this enables interactive plots; remove to plot to PDF figure 
      ggplot(pc.scores.refs.ingrid,
             aes(col=location1, 
#             aes(col=marine_region2, # Alternative line for the grouping variable  
                    text=number)) + 
        geom_point(aes_string(x=axes[i,"pc.x"], y=axes[i,"pc.y"]), size=2, alpha=0.85) +  
        theme_minimal() + ggtitle("Global gene expression PC1xPC2") + 
        ylab(paste(axes[i, "pc.y"], " (", round(variance.ingrid[variance.ingrid$axis==axes[i, "pc.y"], "variance"], digits = 2), "%)", sep="")) + 
        xlab(paste(axes[i, "pc.x"], " (", round(variance.ingrid[variance.ingrid$axis==axes[i, "pc.x"], "variance"], digits = 2), "%)", sep="")) + 
        theme(legend.position = "right", legend.text=element_text(size=8), legend.title=element_text(size=9)) + 
        ggtitle(paste(axes[i, "pc.y"], "x", axes[i, "pc.x"], sep=" ")), tooltip = c("location1", "number")) 
}

#pdf('pcas-genomewide.pdf', width = 8.5, height = 7) #remove "ggplotly" and comment out to save plots to PDF
pcas.genome.ingrid #call plots 
#dev.off() #comment out to save plots to PDF
