## Pacific cod juvenile temperature experiment 
This is an experimental and data summary. See the Google doc [Pcod juvenile temp bioinformatics](https://docs.google.com/document/d/1gPv5Vgm6pS3Q3O4KS7C7OvbgQm7-C63xDKtxU6bkhv8/edit?usp=sharing) for details. 

### Experiment summary
Broad questions: 
- Which biological processes are affected by 2-month exposure to varying temperatures in Pacific cod?  
- What is the adaptive potential of Pacific cod to temperature changes? 
- Are there processes, gene expression levels, or gene variants that predict “high performance” in or “resilience” to suboptimal temperatures? 
    - We need a “high performance” index derived from lipids, growth, hsi data. 

#### Design:
- Juveniles collected in early fall 2022 off Kodiak   
- Four temperature treatments: 0, 5, 9, 16  
- Exposures conducted Winter 2022-2023, fish euthanized February 2023  
- Biometrics collected throughout on per-fish basis (fish ID'd using tags) - length, weight, liver weight   
- Tissues collected: Livers & muscle for lipid analysess; Liver, gill, caudal fin, blood, & spleen for 'omics analysis     

#### Genetic analysis, lcWGS
- lcWGS performed using fin clips from all 160 fish (including three (four?) that died).
  - lcWGS data archived on SednaGold (see [AFSC Genomics data inventory](https://docs.google.com/spreadsheets/d/1RU_aFByYbsXMEIrMgU30YVLwf2AEpjKs8W2vWBT79R8/edit?usp=sharing))
  - Analysis code and output located in [lcWGS](../lcWGS). Any work in R is located in [lcWGS-analysis.Rmd](../lcWGS/notebooks/lcWGS-analysis.Rmd).  
  - AFSC lcWGS pipeline used to analyze genotype data to determine whether more than one source population contributed to experimental fish, and possibly to identify spawning source population(s) 
  - Ran lcWGS pipeline with reference fish of known origin to produce [PCAs for each chromosome/whole genome](../lcWGS/analysis-20230922/pca), but batch effects were a problem, did not allow for confident population assignment. PCA did not indicate any major genetic differences among experimental treatments, which is good
    - FYI: PCA indicates three “outliers” (sample prep issue?)  -  GM1, GM2, & GM121. 
  - Did a quick-n-dirty ZP3 haplotype call from lcWGS (using ANGSD) and RNASeq data (manually, using IGV), compared haplotypes from each approach and found good agreement. Based on these ZP3 haplotypes it seems there could be some fish from the Eastern GOA.  
  - B/C ZP3 haplotypes indicate possible population diversity, I am using the AFSC lcWGS pipeline + WGSassign along with lcWGS data from hundreds of "reference" fish (known origins) to "assign" each experimental fish to a population

#### Gene expression analysis, RNASeq  
- RNASeq of liver tissue conducted on 20 fish from each temperature treatment
  - I ran my RNASeq pipeline.
  - MORE HERE SOON  
  - I also looked for ZP3 haplotypes in the three exonic variant loci - details below! 
