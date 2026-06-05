rm(list=ls())
setwd("C:/Users/aniaf/Projects/BlackSpruce/945_GEA_lfmm")
dir()

library(tidyverse)
library(cowplot)
library(readxl)
library(stringr)
library(LEA)



#------------#
#    META    #
#------------#

dmeta <- read.csv("../DATA_intermediate/23_filter_EPN_indiv_metrics.tsv", sep="\t", header=T)
dmeta$POP <- as.character(dmeta$POP)
dmeta$POP_SITE <- paste0(dmeta$POP, "_", dmeta$SITE_ID)
dmeta_sub <- dmeta %>% dplyr::select(POP, STRATA, SITE_ID, POP_SITE, lat, lon) %>% distinct() %>% arrange(POP_SITE)
head(dmeta_sub)


dclim <- read.csv("../34_climate_transfer_distance/03_climate_distance_PCA.tsv", sep="\t", header=T)
dclim3 <- dclim[c("Row.names","PC1","PC2","PC3")]
dim(dclim3)
head(dclim3)




#-------------------#
#    ONE EXAMPLE    #
#-------------------#

# 1. Read file with pop frequencies
# 2. A function to create input files for lfmm

datafile <- "01_frequencies_POP_Height.tsv"


### One file example

dataset_name <- "Height"
i <- 1

# Reading frequencies
df <- read.csv(datafile, sep="\t", header=T, row.names=1)
print(dim(df))
df[1:5, 1:5]

# Filtering SNPs with MAF < 0.05 and > 0.95
# Average Freq per SNP
vec.freqs <- colSums(df)/(nrow(df))
hist(vec.freqs, breaks=50)
MAF <- 0.05
comb.freqs <- c(which(vec.freqs < MAF), which(vec.freqs > (1-MAF)))
names(comb.freqs)

df_filtered <- df %>% select(-all_of(names(comb.freqs)))
dim(df)
dim(df_filtered)

### Saving SNP names
dsnps <- data.frame("snp" = colnames(df_filtered))
write.table(dsnps, file = paste0("02_lea_inputs_",dataset_name,".snps"), sep="\t", row.names = F, col.names = T, quote=F, append=F)

### Saving frequencies in lfmm format
#df2[is.na(df2)] <- 9
write.table(df, file = paste0("02_lea_inputs_",dataset_name,".tsv"), sep="\t", row.names = T, col.names = T, quote=F, append=F)


### Getting metadata

samples <- data.frame("sample"=row.names(df_filtered))
samples$num <- c(1:length(row.names(df_filtered)))
smeta <- merge(samples, dmeta_sub, by.x = "sample", by.y = "POP_SITE", sort=F) 
cmeta <- merge(smeta, dclim3, by.x = "POP", by.y = "Row.names", sort=F)
write.table(cmeta, file = paste0("02_lea_inputs_metadata_",dataset_name,".tsv"), sep="\t", row.names = F, col.names = T, quote=F, append=F)





#-------------------#
#        LOOP       #
#-------------------#


make_inputs <- function(datafile, dataset_name, seed) {
  
  # Reading frequencies
  df <- read.csv(datafile, sep="\t", header=T, row.names=1)
  print(dim(df))
  
  # Filtering SNPs with MAF < 0.05 and > 0.95
  vec.freqs <- colSums(df)/(nrow(df))
  MAF <- 0.05
  comb.freqs <- c(which(vec.freqs < MAF), which(vec.freqs > (1-MAF)))
  df_filtered <- df %>% select(-all_of(names(comb.freqs)))
  
  ### Saving SNP names
  dsnps <- data.frame("snp" = colnames(df_filtered))
  write.table(dsnps, file = paste0("02_lea_inputs_",dataset_name,".snps"), sep="\t", row.names = F, col.names = T, quote=F, append=F)
  
  ### Saving frequencies in lfmm format
  #df[is.na(df)] <- 9
  write.table(df_filtered, file = paste0("02_lea_inputs_",dataset_name,".tsv"), sep="\t", row.names = T, col.names = T, quote=F, append=F)

  ### Getting metadata
  
  samples <- data.frame("sample"=row.names(df_filtered))
  samples$num <- c(1:length(row.names(df_filtered)))
  smeta <- merge(samples, dmeta_sub, by.x = "sample", by.y = "POP_SITE", sort=F) 
  cmeta <- merge(smeta, dclim3, by.x = "POP", by.y = "Row.names", sort=F)
  write.table(cmeta, file = paste0("02_lea_inputs_metadata_",dataset_name,".tsv"), sep="\t", row.names = F, col.names = T, quote=F, append=F)
  
  
  
}


#-----------------#
#     FILES       #
#-----------------#



datafile <- "01_frequencies_POP_Height.tsv"
i <- 1
for (trait in c("Height", "Biomass_Increment", "Biomass_Increment_1980", "Biomass_Increment_1985", 
                "Biomass_Increment_1990", "Biomass_Increment_1995", "Biomass_Increment_2000", 
                "Biomass_Increment_2005", "Biomass_Increment_2010", "Biomass_Increment_2015",
                "Average_Ring_Density", "DBH", "Rs", "Rr", "Rl", "Rc")) {
  
  filename <- paste0("01_frequencies_POP_",trait,".tsv")
  print(filename)
  print(i)
  make_inputs(filename, trait, i)
  i = i+1
}










#-----------------#
#       END       #
#-----------------#



