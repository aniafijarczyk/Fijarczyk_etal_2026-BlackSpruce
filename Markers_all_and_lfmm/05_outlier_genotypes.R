rm(list=ls())
setwd("C:/Users/aniaf/Projects/BlackSpruce/945_GEA_lfmm")
dir()

library(tidyverse)
library(cowplot)
library(readxl)
library(stringr)


#------------#
#    DATA    #
#------------#

trait_name <- "Height"
filename <- paste0("03_lfmm_results_",trait_name,".tsv")
df <- read.csv(filename, sep="\t", header=T)
head(df)  
df %>% filter(adj.p.PC1<0.05 | adj.p.PC2<0.05 | adj.p.PC3<0.05)
snps <- df %>% filter(adj.p.PC1<0.05 | adj.p.PC2<0.05 | adj.p.PC3<0.05) %>% dplyr::select(snp) %>% pull()
length(snps)


Y <- read.csv(paste0("../945_GEA_lfmm/02_lea_inputs_",trait_name,".tsv"), header=T, sep="\t")
Y[1:5,1:5]
dim(Y)
Y.sub <- Y[,snps]
dim(Y.sub)
Y.sub[1:5,1:5]

write.table(Y.sub, file = paste0("05_outlier_genotypes_lfmm_",trait_name,".tsv"), sep="\t", row.names = T, col.names = T, quote=F, append=F)


#------------#
#    LOOP    #
#------------#

write_inputs <- function(trait_name, seed) {
  
  ### Data
  Y <- read.csv(paste0("../945_GEA_lfmm/02_lea_inputs_",trait_name,".tsv"), header=T, sep="\t")
  filename <- paste0("03_lfmm_results_",trait_name,".tsv")

  ### Subset
  df <- read.csv(filename, sep="\t", header=T)
  snps <- df %>% filter(adj.p.PC1<0.05 | adj.p.PC2<0.05 | adj.p.PC3<0.05) %>% dplyr::select(snp) %>% pull()
  Y.sub <- Y[,snps]

  ### Write
  write.table(Y.sub, file = paste0("05_outlier_genotypes_lfmm_",trait_name,".tsv"), sep="\t", row.names = T, col.names = T, quote=F, append=F)

  ### Getting random subsets
  all_snps <- colnames(Y)
  
  set.seed(778343+seed)
  random_snp <- sample(all_snps, 10000)
  Y.random <- Y[,random_snp]
  print(dim(Y.random))
  write.table(Y.random, paste0("05_outlier_genotypes_10000_",trait_name,".tsv"), sep="\t", row.names = T, col.names = T, quote=F, append=F)
  
  set.seed(778343+seed+100)
  random_snp <- sample(all_snps, 1000)
  Y.random <- Y[,random_snp]
  print(dim(Y.random))
  write.table(Y.random, paste0("05_outlier_genotypes_1000_",trait_name,".tsv"), sep="\t", row.names = T, col.names = T, quote=F, append=F)
  
  set.seed(778343+seed+200)
  random_snp <- sample(all_snps, 100)
  Y.random <- Y[,random_snp]
  print(dim(Y.random))
  write.table(Y.random, paste0("05_outlier_genotypes_100_",trait_name,".tsv"), sep="\t", row.names = T, col.names = T, quote=F, append=F)
  
  
  
}



##########

i <- 1
for (trait in c("Height", "Biomass_Increment", "Biomass_Increment_1980", "Biomass_Increment_1985", 
                "Biomass_Increment_1990", "Biomass_Increment_1995", "Biomass_Increment_2000", 
                "Biomass_Increment_2005", "Biomass_Increment_2010", "Biomass_Increment_2015",
                "Average_Ring_Density", "DBH", "Rs", "Rr", "Rl", "Rc")) {
  
  print(trait)
  print(i)
  write_inputs(trait, i)
  i = i+1
  
}





#-----------------#
#       END       #
#-----------------#



