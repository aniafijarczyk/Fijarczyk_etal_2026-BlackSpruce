rm(list=ls())
setwd("C:/Users/aniaf/Projects/BlackSpruce/947_GEA_RDA")
dir()

library(tidyverse)
library(cowplot)
library(readxl)
library(stringr)


#------------#
#    DATA    #
#------------#

trait_name <- "Height"
filename <- paste0("01_rda_results_",trait_name,"_noPCA.tsv")
df <- read.csv(filename, sep="\t", header=T)
head(df)  
snps <- df %>% filter(adj.p<0.05) %>% dplyr::select(snp) %>% pull()
length(snps)

Y <- read.csv(paste0("../945_GEA_lfmm/02_lea_inputs_",trait_name,".tsv"), header=T, sep="\t")
Y[1:5,1:5]
dim(Y)
Y.sub <- Y[,snps]
dim(Y.sub)
Y.sub[1:5,1:5]

write.table(Y.sub, file = paste0("02_outlier_genotypes_RDA_",trait_name,".tsv"), sep="\t", row.names = T, col.names = T, quote=F, append=F)


filename2 <- paste0("01_rda_results_",trait_name,"_PCAcorrected.tsv")
df2 <- read.csv(filename2, sep="\t", header=T)
snps2 <- df2 %>% filter(adj.p<0.05) %>% dplyr::select(snp) %>% pull()
Y.sub2 <- Y[,snps2]
dim(Y.sub2)
Y.sub2[1:5,1:5]

write.table(Y.sub2, file = paste0("02_outlier_genotypes_RDAcorrected_",trait_name,".tsv"), sep="\t", row.names = T, col.names = T, quote=F, append=F)



#------------#
#    LOOP    #
#------------#

write_inputs <- function(trait_name) {
  
  ### Data
  Y <- read.csv(paste0("../945_GEA_lfmm/02_lea_inputs_",trait_name,".tsv"), header=T, sep="\t")
  filename <- paste0("01_rda_results_",trait_name,"_noPCA.tsv")
  filename2 <- paste0("01_rda_results_",trait_name,"_PCAcorrected.tsv")
  
  ### Subset
  df <- read.csv(filename, sep="\t", header=T)
  snps <- df %>% filter(adj.p<0.05) %>% dplyr::select(snp) %>% pull()
  Y.sub <- Y[,snps]

  ### Write
  write.table(Y.sub, file = paste0("02_outlier_genotypes_RDA_",trait_name,".tsv"), sep="\t", row.names = T, col.names = T, quote=F, append=F)
  
  ### Subset2
  df2 <- read.csv(filename2, sep="\t", header=T)
  snps2 <- df2 %>% filter(adj.p<0.05) %>% dplyr::select(snp) %>% pull()
  Y.sub2 <- Y[,snps2]
  
  ### Write
  write.table(Y.sub2, file = paste0("02_outlier_genotypes_RDAcorrected_",trait_name,".tsv"), sep="\t", row.names = T, col.names = T, quote=F, append=F)
  
}



##########

for (trait in c("Height", "Biomass_Increment", "Biomass_Increment_1980", "Biomass_Increment_1985", 
                "Biomass_Increment_1990", "Biomass_Increment_1995", "Biomass_Increment_2000", 
                "Biomass_Increment_2005", "Biomass_Increment_2010", "Biomass_Increment_2015",
                "Average_Ring_Density", "DBH", "Rs", "Rr", "Rl", "Rc")) {
  
  print(trait)
  write_inputs(trait)
}





#-----------------#
#       END       #
#-----------------#



