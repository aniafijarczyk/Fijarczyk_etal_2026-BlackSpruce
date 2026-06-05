rm(list=ls())
setwd("C:/Users/aniaf/Projects/BlackSpruce/947_GEA_RDA")
dir()

library(tidyverse)
library(cowplot)
library(readxl)
library(stringr)
library(vegan)
library(robust)
library(ade4)




#------------#
#    DATA    #
#------------#

# genotypes 

trait <- "Height"
dmeta <- read.csv(paste0("../945_GEA_lfmm/02_lea_inputs_metadata_",trait,".tsv"), sep="\t", header=T)
head(dmeta)
dim(dmeta)

# SNPs filtered for GEA
dsnps <- read.csv(paste0("../945_GEA_lfmm/02_lea_inputs_",trait,".snps"), sep="\t", header=T)
head(dsnps)
dim(dsnps)



#-----------------------#
#    RDA one example    #
#-----------------------#

# Dataset with snps
Y <- read.csv(paste0("../945_GEA_lfmm/02_lea_inputs_",trait,".tsv"), header=T, sep="\t")
Y[1:5, 1:5]

# Checking if MAF is ok - should be ok becouse low freq SNPs were filtered out
frequences <- colSums(Y)/(nrow(Y))
hist(frequences)
MAF <- 0.05
length(which(frequences < MAF))
length(which(frequences > (1-MAF)))

# Metadata
X <- read.csv(paste0("../45_GEA_lfmm/02_lea_inputs_metadata_",trait,".tsv"), header=T, sep="\t")
head(X)
X$sample == rownames(Y)
print(sum(X$sample == rownames(Y))/length(X$sample))

rownames(X) <- X$sample
X.env <- X[c("PC1","PC2","PC3")]
head(X.env)

dim(Y)
dim(X)

# Fitting an LFMM with K = 3 factors
# lfmm2() command generates an object of class lfmm2Class which
#contains estimated factors (mod2@U) and loadings (mod2@V) for being intro-
#  duced as correction factors in genome-wide association tests.
mod <- rda(Y ~ X.env$PC1 + X.env$PC2 + X.env$PC3, X.env)
mod
summary(mod)
# Our constrained ordination (environment) explains about 20% of the variation
RsquareAdj(mod)
# The eigenvalues for the constrained axes reflect the variance explained by each canonical axis
summary(eigenvals(mod, model = "constrained"))
screeplot(mod)

# check our RDA model for significance using formal tests
# We can assess both the full model and each constrained axis using F-statistics (Legendre et al, 2010)
signif.full <- anova.cca(mod, parallel=getOption("mc.cores")) # default is permutation=999
signif.full

# For this test, each constrained axis is tested using all previous constrained axes as conditions
signif.axis <- anova.cca(mod, by="axis", parallel=getOption("mc.cores"))
signif.axis

# RDA1       1    84.15 28.4495  0.001 ***
# RDA2       1     9.74  3.2924  0.001 ***
# RDA3       1     3.73  1.2619  0.065 .  
# Residual 111   328.34                   


# checking Variance Inflation Factors for the predictor variables used in the model:
vif.cca(mod)
plot(mod, scaling=3)
plot(mod, choices = c(1, 3), scaling=3)

mod$CCA$v %>% head()
mod$CA$v

ggplot() +
  geom_line(aes(x=c(1:length(mod$CCA$eig)), y=as.vector(mod$CCA$eig)), linetype="dotted",
            size = 1.5, color="darkgrey") +
  geom_point(aes(x=c(1:length(mod$CCA$eig)), y=as.vector(mod$CCA$eig)), size = 3,
             color="darkgrey") +
  scale_x_discrete(name = "Ordination axes", limits=c(1:9)) +
  ylab("Inertia") +
  theme_bw()


rdadapt<-function(rda,K)
{
  zscores<-rda$CCA$v[,1:as.numeric(K)]
  resscale <- apply(zscores, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  #qval <- qvalue(reschi2test)
  qval <- p.adjust(reschi2test, method = "bonferroni")
  q.values_rdadapt<-qval
  return(data.frame(p.values=reschi2test, q.values=q.values_rdadapt))
}

res_rdadapt <- rdadapt(mod, 3)
res_rdadapt %>% head()
res_rdadapt$p.values

as.data.frame(mod$CCA$v)$RDA1
as.vector(res_rdadapt$p.values)
as.vector(res_rdadapt$q.values)
mod$CCA$v %>% head()

dres <- data.frame("snp" = colnames(Y),
                   "pvalues" = as.vector(res_rdadapt$p.values),
                   "adj.p" = as.vector(res_rdadapt$q.values),
                   "zscores.RDA1" = as.vector(as.data.frame(mod$CCA$v)$RDA1),
                   "zscores.RDA2" = as.vector(as.data.frame(mod$CCA$v)$RDA2),
                   "zscores.RDA3" = as.vector(as.data.frame(mod$CCA$v)$RDA3)
)
head(dres)
cutoff <- 0.01/length(colnames(Y))
cutoff
dres %>% filter(pvalues < cutoff) %>% dim()
dres %>% filter(adj.p < 0.01) %>% dim()

ggplot() +
  geom_point(aes(x=mod$CCA$v[,1], y=mod$CCA$v[,2]), col = "gray86") +
  geom_point(aes(x=mod$CCA$v[which(res_rdadapt[,2] < 0.05),1],
                 y=mod$CCA$v[which(res_rdadapt[,2] < 0.05),2]), col = "orange") +
  geom_segment(aes(xend=mod$CCA$biplot[,1]/10, yend=mod$CCA$biplot[,2]/10, x=0, y=0),
               colour="black", size=0.5, linetype=1,
               arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(aes(x=1.2*mod$CCA$biplot[,1]/10, y=1.2*mod$CCA$biplot[,2]/10,
                label = colnames(X.env[,1:3]))) +
  xlab("RDA 1") + ylab("RDA 2") +
  theme_bw() +
  theme(legend.position="none")


########### With structure
Y[1:5, 1:5]
pca.snps <- dudi.pca(Y, scale = FALSE, scannf = FALSE, nf = 6)
pca.2.axes <- pca.snps$li[c(1,2)]
head(pca.2.axes)
plot(pca.2.axes$Axis1, pca.2.axes$Axis2)

head(X.env)
X.env$Axis1 <- pca.2.axes$Axis1
X.env$Axis2 <- pca.2.axes$Axis2

mod.2 <- rda(Y ~ X.env$PC1 + X.env$PC2 + X.env$PC3 + Condition(Axis1 + Axis2), X.env)
mod.2
summary(mod.2)
# Our constrained ordination (environment) explains about 20% of the variation
RsquareAdj(mod.2)
# The eigenvalues for the constrained axes reflect the variance explained by each canonical axis
summary(eigenvals(mod.2, model = "constrained"))
screeplot(mod.2)

anova.cca(mod.2)




################################\
#------------#
#    LOOP    #
#------------#


run_rda <- function(trait_name) {
  
  print(trait_name)
  
  ### read input files
  X <- read.csv(paste0("../945_GEA_lfmm/02_lea_inputs_metadata_",trait_name,".tsv"), header=T, sep="\t")
  Y <- read.csv(paste0("../945_GEA_lfmm/02_lea_inputs_",trait_name,".tsv"), header=T, sep="\t")
  # Check if names of pops are in order, should be 1
  print(sum(X$sample == rownames(Y))/length(X$sample))
  rownames(X) <- X$sample
  X.env <- X[c("PC1","PC2","PC3")]

  ### run rda
  mod <- rda(Y ~ X.env$PC1 + X.env$PC2 + X.env$PC3, X.env)
  
  ### pvalues - each PC separately
  res_rdadapt <- rdadapt(mod, 3)

  ### Creating output table
  dres <- data.frame("snp" = colnames(Y),
                     "pvalues" = as.vector(res_rdadapt$p.values),
                     "adj.p" = as.vector(res_rdadapt$q.values),
                     "adj.r.squared" = RsquareAdj(mod),
                     "zscores.RDA1" = as.vector(as.data.frame(mod$CCA$v)$RDA1),
                     "zscores.RDA2" = as.vector(as.data.frame(mod$CCA$v)$RDA2),
                     "zscores.RDA3" = as.vector(as.data.frame(mod$CCA$v)$RDA3)
  )

  ### Writing output
  write.table(dres, paste0("01_rda_results_",trait_name,"_noPCA.tsv"), sep="\t", row.names = F, col.names = T, quote=F, append=F)
  
  
}




for (trait in c("Height", "Biomass_Increment", "Biomass_Increment_1980", "Biomass_Increment_1985", 
                "Biomass_Increment_1990", "Biomass_Increment_1995", "Biomass_Increment_2000", 
                "Biomass_Increment_2005", "Biomass_Increment_2010", "Biomass_Increment_2015",
                "Average_Ring_Density", "DBH", "Rs", "Rr", "Rl", "Rc")) {
  
  print(trait)
  run_rda(trait)
}





#---------------------#
#    LOOP CORRECTED   #
#---------------------#

run_rda_corrected <- function(trait_name) {
  
  print(trait_name)
  
  ### read input files
  X <- read.csv(paste0("../45_GEA_lfmm/02_lea_inputs_metadata_",trait_name,".tsv"), header=T, sep="\t")
  Y <- read.csv(paste0("../45_GEA_lfmm/02_lea_inputs_",trait_name,".tsv"), header=T, sep="\t")
  
  # Check if names of pops are in order, should be 1
  print(sum(X$sample == rownames(Y))/length(X$sample))
  rownames(X) <- X$sample
  X.env <- X[c("PC1","PC2","PC3")]

  # Structure
  pca.snps <- dudi.pca(Y, scale = FALSE, scannf = FALSE, nf = 6)
  pca.2.axes <- pca.snps$li[c(1,2)]
  X.env$Axis1 <- pca.2.axes$Axis1
  X.env$Axis2 <- pca.2.axes$Axis2
  
  ### run rda
  mod.2 <- rda(Y ~ X.env$PC1 + X.env$PC2 + X.env$PC3 + Condition(Axis1 + Axis2), X.env)
  
  ### pvalues - each PC separately
  res_rdadapt <- rdadapt(mod.2, 3)
  
  ### Creating output table
  dres <- data.frame("snp" = colnames(Y),
                     "pvalues" = as.vector(res_rdadapt$p.values),
                     "adj.p" = as.vector(res_rdadapt$q.values),
                     "adj.r.squared" = RsquareAdj(mod.2),
                     "zscores.RDA1" = as.vector(as.data.frame(mod.2$CCA$v)$RDA1),
                     "zscores.RDA2" = as.vector(as.data.frame(mod.2$CCA$v)$RDA2),
                     "zscores.RDA3" = as.vector(as.data.frame(mod.2$CCA$v)$RDA3)
  )
  
  ### Writing output
  write.table(dres, paste0("01_rda_results_",trait_name,"_PCAcorrected.tsv"), sep="\t", row.names = F, col.names = T, quote=F, append=F)
  
  
}



for (trait in c("Height", "Biomass_Increment", "Biomass_Increment_1980", "Biomass_Increment_1985", 
                "Biomass_Increment_1990", "Biomass_Increment_1995", "Biomass_Increment_2000", 
                "Biomass_Increment_2005", "Biomass_Increment_2010", "Biomass_Increment_2015",
                "Average_Ring_Density", "DBH", "Rs", "Rr", "Rl", "Rc")) {
  
  print(trait)
  run_rda_corrected(trait)
}


for (trait in c("Biomass_Increment_1980")) {
  
  print(trait)
  run_rda_corrected(trait)
}




#################################################################################
