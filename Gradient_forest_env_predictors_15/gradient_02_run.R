rm(list=ls())
setwd("C:/Users/aniaf/Projects/BlackSpruce/933_gradient_forest_env_15")

library(gradientForest)
library(tidyverse)
library(cowplot)
library(ggrepel)
library(readxl)
library(RColorBrewer)
library(dplyr)

help(gradientForest)
# The data should not include NAs.





################
#---  DATA  ---#
################

df.sub <- read.csv("./inputs/gradient_01_data_prep_input_AF_random.csv", sep=",", header=TRUE)
dim(df.sub)

env.mem <- read.csv("./inputs/gradient_01_data_prep_input_variables.csv", sep=",", header=TRUE)
head(env.mem)
dim(env.mem)


names(which(colSums(is.na(df.sub))>0))
sum(colSums(is.na(df.sub)))
rownames(df.sub) == rownames(env.mem)



####################################################################
### Examining compositional change along environmental gradients ###
####################################################################


####################################################################
# Gradient Forest random SNPs

# Setting max level
nSites <- dim(df.sub)[1]
nSpecs <- dim(df.sub)[2]
lev <- floor(log2(nSites * 0.368/2))
# Running gradient forest
gf <- gradientForest(data = cbind(env.mem, df.sub), predictor.vars = colnames(env.mem),
                     response.vars = colnames(df.sub), ntree = 500, transform = NULL,
                     corr.threshold = 0.5, compact = T, nbin = 201, maxLevel = lev)


gf

######################
#---  SAVING  -------#
######################


save(gf, file = "gradient_02_run_results_random.RData")





