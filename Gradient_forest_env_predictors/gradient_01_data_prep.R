rm(list=ls())
setwd("C:/Users/aniaf/Projects/BlackSpruce/932_gradient_forest_env")

library(gradientForest)
library(tidyverse)
library(cowplot)
library(ggrepel)
library(readxl)
library(RColorBrewer)
library(dplyr)
library(stringr)

help(gradientForest)
# The data should not include NAs.




#########################
#-------  AFs ----------#
#########################

### Reading allele frequencies
dgeno <- read.csv("../926_MoransI/01_moran_Pop_AFs.csv", sep=",", header=T)
dgeno[1:10, 1:10]
dim(dgeno)
# Removing invariant snps
zero.snps <- colnames(dgeno)[colSums(dgeno)==0]
one.snps <- colnames(dgeno)[colSums(dgeno)==nrow(dgeno)]
dgeno <- dgeno %>% dplyr::select(!all_of(zero.snps))
dgeno <- dgeno %>% dplyr::select(!all_of(one.snps))
dim(dgeno)



pops <- row.names(dgeno)
snps <- colnames(dgeno)
pops
snps
length(snps)
length(pops)

########################
#--SUBSET RANDOM 1000--#
########################

set.seed(21102024)
dgeno.sub <- dgeno[,sample(ncol(dgeno), 1000)]
dim(dgeno.sub)
dgeno.sub[1:10, 1:10]

# Saving AFs for input
write.table(dgeno.sub, "./inputs/gradient_01_data_prep_input_AF_random.csv",sep=',',col.names = TRUE, row.names = TRUE, quote = FALSE, append = FALSE)



########################
#--ENVIRONMENTAL DATA--#
########################

dclim <- read.csv("../31_climate/02_garden_climate_with_gardens.tab", sep="\t", header=TRUE)
head(dclim)
dclim <- dclim %>% dplyr::filter(KeyID %in% pops)
rownames(dclim) <- dclim$KeyID
dclim <- dclim %>% dplyr::select(-Latitude, -Longitude, -Elevation, -KeyID, -POP_GR, -dataset)
dclim %>% head()
row.names(dclim)


##########
#--MEMs--#
##########

dmem <- read.csv("../926_MoransI/01_moran_MEMs.tsv", sep="\t", header=T)
dmem %>% head()
rownames(dmem)
dim(dmem)
dmem <- dmem %>% dplyr::select(MEM1, MEM2, MEM3)
head(dmem)
row.names(dmem)


# Combining environmental vars and memes
row.names(dmem)
row.names(dclim)

env.mem <- cbind(dclim, dmem)
dim(env.mem)
head(env.mem)


# Saving variables for input
write.table(env.mem, "./inputs/gradient_01_data_prep_input_variables.csv", sep=',',col.names = TRUE, row.names = TRUE, quote = FALSE, append = FALSE)







