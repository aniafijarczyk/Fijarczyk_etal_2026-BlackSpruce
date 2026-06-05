rm(list=ls())
setwd("C:/Users/aniaf/Projects/BlackSpruce/945_GEA_lfmm")
dir()

library(tidyverse)
library(cowplot)
library(readxl)
library(stringr)
library(LEA)



#------------#
#    DATA    #
#------------#

# Example with one file

trait <- "Height"
dmeta <- read.csv(paste0("02_lea_inputs_metadata_",trait,".tsv"), sep="\t", header=T)
head(dmeta)
dim(dmeta)

dsnps <- read.csv(paste0("02_lea_inputs_",trait,".snps"), sep="\t", header=T)
head(dsnps)
dim(dsnps)



#------------#
#   LFMM2    #
#------------#


# Preparing input

Y <- read.csv(paste0("02_lea_inputs_",trait,".tsv"), header=T, sep="\t")
Y[1:5, 1:5]

X <- read.csv(paste0("02_lea_inputs_metadata_",trait,".tsv"), header=T, sep="\t")
head(X)
X$sample == rownames(Y)
print(sum(X$sample == rownames(Y))/length(X$sample))

rownames(X) <- X$sample
X.env <- X[c("PC1","PC2","PC3")]
head(X.env)

# Fitting an LFMM with K = 3 factors
# lfmm2() command generates an object of class lfmm2Class which
# contains estimated factors (mod2@U) and loadings (mod2@V) for being intro-
# duced as correction factors in genome-wide association tests.

mod <- lfmm2(input = Y, env = X.env, K = 3)

# Latent factors
mod@U %>% head()
# Loadings
mod@V %>% head()
dim(mod@V)


# Computing P-values and plotting their minus log10 values
# The function returns a vector of p-values for association between loci 
# and environmental variables adjusted for latent factors computed by lfmm2.
# full = TRUE - p-values for the full set of environmental variables
# full = FALSE - (p-values) for each environmental variable

pv <- lfmm2.test(object = mod,
                 input = Y,
                 env = X.env,
                 linear = TRUE,
                 )

pv$pvalues[1:3, 1:10]
pv$zscores[1:3, 1:10]
pv$adj.r.squared
as.vector(pv$zscores[1,])
#  genomic inflation factors computed for each environmental variable
pv$gif




### Full - all environmental variables included in the model

pv.full <- lfmm2.test(object = mod,
                 input = Y,
                 env = X.env,
                 linear = TRUE,
                 full = TRUE
)

pv.full$pvalues[1:10]
pv.full$fscores
#  genomic inflation factors computed for each environmental variable
pv.full$gif
# a vector of R squared values or variances explained by all environmental variables for all loci.
pv.full$adj.r.squared

sum(pv.full$adj.r.squared)
hist(pv.full$adj.r.squared)
pv.full$adj.r.squared[pv.full$adj.r.squared>0.3]


### Create output table
colnames(Y)

dres <- data.frame("snp" = dsnps$snp,
                   "snp_id" = colnames(Y),
                   "pvalues" = as.vector(pv.full$pvalues),
                   "fscores" = as.vector(pv.full$fscores),
                   "adj.r.sq" = as.vector(pv.full$adj.r.squared),
                   "pvals.PC1" = as.vector(pv$pvalues[1,]),
                   "pvals.PC2" = as.vector(pv$pvalues[2,]),
                   "pvals.PC3" = as.vector(pv$pvalues[3,]),
                   "zscores.PC1" = as.vector(pv$zscores[1,]),
                   "zscores.PC2" = as.vector(pv$zscores[2,]),
                   "zscores.PC3" = as.vector(pv$zscores[3,])                   
                   )
head(dres)

dres$adj.p <- p.adjust(dres$pvalues, method = "bonferroni")
dres$adj.p.PC1 <- p.adjust(dres$pvals.PC1, method = "bonferroni")
dres$adj.p.PC2 <- p.adjust(dres$pvals.PC2, method = "bonferroni")
dres$adj.p.PC3 <- p.adjust(dres$pvals.PC3, method = "bonferroni")
head(dres)


#write.table(dres, "03_lfmm_results.tsv", sep="\t", row.names = F, col.names = T, quote=F, append=F)



#------------#
#    LOOP    #
#------------#


run_lfmm <- function(trait_name) {
  
  print(trait_name)
  
  ### read input files
  dmeta <- read.csv(paste0("02_lea_inputs_metadata_",trait_name,".tsv"), sep="\t", header=T)
  dsnps <- read.csv(paste0("02_lea_inputs_",trait_name,".snps"), sep="\t", header=T)
  Y <- read.csv(paste0("02_lea_inputs_",trait,".tsv"), header=T, sep="\t")
  X <- read.csv(paste0("02_lea_inputs_metadata_",trait,".tsv"), header=T, sep="\t")
  # Check if names of pops are in order, should be 1
  print(sum(X$sample == rownames(Y))/length(X$sample))
  rownames(X) <- X$sample
  X.env <- X[c("PC1","PC2","PC3")]
  
  ### run lfmms
  mod <- lfmm2(input = Y, env = X.env, K = 3)
  
  ### pvalues - each PC separately
  pv <- lfmm2.test(object = mod,
                   input = Y,
                   env = X.env,
                   linear = TRUE,
  )
  
  ### pvalues - full - variables together
  pv.full <- lfmm2.test(object = mod,
                        input = Y,
                        env = X.env,
                        linear = TRUE,
                        full = TRUE
  )
  
  ### Creating output table
  dres <- data.frame("snp" = dsnps$snp,
                     "snp_id" = colnames(Y),
                     "pvalues" = as.vector(pv.full$pvalues),
                     "fscores" = as.vector(pv.full$fscores),
                     "adj.r.sq" = as.vector(pv.full$adj.r.squared),
                     "pvals.PC1" = as.vector(pv$pvalues[1,]),
                     "pvals.PC2" = as.vector(pv$pvalues[2,]),
                     "pvals.PC3" = as.vector(pv$pvalues[3,]),
                     "zscores.PC1" = as.vector(pv$zscores[1,]),
                     "zscores.PC2" = as.vector(pv$zscores[2,]),
                     "zscores.PC3" = as.vector(pv$zscores[3,])                   
  )

  dres$adj.p <- p.adjust(dres$pvalues, method = "bonferroni")
  dres$adj.p.PC1 <- p.adjust(dres$pvals.PC1, method = "bonferroni")
  dres$adj.p.PC2 <- p.adjust(dres$pvals.PC2, method = "bonferroni")
  dres$adj.p.PC3 <- p.adjust(dres$pvals.PC3, method = "bonferroni")
  head(dres)
  
  ### Writing output
  write.table(dres, paste0("03_lfmm_results_",trait_name,".tsv"), sep="\t", row.names = F, col.names = T, quote=F, append=F)
  

}



for (trait in c("Height", "Biomass_Increment", "Biomass_Increment_1980", "Biomass_Increment_1985", 
                "Biomass_Increment_1990", "Biomass_Increment_1995", "Biomass_Increment_2000", 
                "Biomass_Increment_2005", "Biomass_Increment_2010", "Biomass_Increment_2015",
                "Average_Ring_Density", "DBH", "Rs", "Rr", "Rl", "Rc")) {
  
  print(trait)
  run_lfmm(trait)
}









#------------#
#    END     #
#------------#
