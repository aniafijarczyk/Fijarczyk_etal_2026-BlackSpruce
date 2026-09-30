rm(list=ls())
setwd("C:/Users/aniaf/Projects/BlackSpruce/930_RDA")

library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(pals)
library(ade4)
library(adespatial)
library(adegraphics)
library(spdep)
library(dartR)
library(vegan)
library(sp)
library(stringr)





#===================#
#   Genetic data    #
#===================#

# AFs for populations 
dgeno.W <- read.csv("04_AF_by_cluster_West.csv", sep=",", header=T)
dgeno.W[c(1:5), c(1:5)]
rownames(dgeno.W)
dim(dgeno.W)
pops <- rownames(dgeno.W)
df.pops.W <- data.frame("pop" = pops)
head(df.pops.W)

dgeno.C <- read.csv("04_AF_by_cluster_Central.csv", sep=",", header=T)
dgeno.C[c(1:5), c(1:5)]
rownames(dgeno.C)
dim(dgeno.C)
pops <- rownames(dgeno.C)
df.pops.C <- data.frame("pop" = pops)
head(df.pops.C)

dgeno.E <- read.csv("04_AF_by_cluster_East.csv", sep=",", header=T)
dgeno.E[c(1:5), c(1:5)]
rownames(dgeno.E)
dim(dgeno.E)
pops <- rownames(dgeno.E)
df.pops.E <- data.frame("pop" = pops)
head(df.pops.E)


dim(dgeno.W)
dim(dgeno.C)
dim(dgeno.E)



#======================#
#  Environmental data  #
#======================#

# Data
df.env <- read.csv("../31_climate/02_garden_climate_PCA.tsv", sep="\t", header=TRUE)
head(df.env)

# Filter populations
env.pops <- merge(df.pops.W, df.env, by.x = "pop", by.y = "Row.names", sort=F)
row.names(env.pops) <- env.pops$pop
env.W <- env.pops %>% dplyr::select(PC1, PC2, PC3)
head(env.W)
dim(env.W)


env.pops.C <- merge(df.pops.C, df.env, by.x = "pop", by.y = "Row.names", sort=F)
row.names(env.pops.C) <- env.pops.C$pop
env.C <- env.pops.C %>% dplyr::select(PC1, PC2, PC3)
head(env.C)
dim(env.C)


env.pops.E <- merge(df.pops.E, df.env, by.x = "pop", by.y = "Row.names", sort=F)
row.names(env.pops.E) <- env.pops.E$pop
env.E <- env.pops.E %>% dplyr::select(PC1, PC2, PC3)
head(env.E)
dim(env.E)



#======================#
#       Top MEMs       #
#======================#

#===================#
#     Metadata      #
#===================#

# Selecting pops
dmeta <- read.csv("../DATA_intermediate/23_filter_EPN_indiv_metrics.tsv", sep="\t", header=T)
dpops <- dmeta %>% filter(SPECIES_ID == 'EPN') %>% group_by(POP) %>% dplyr::summarise(n=n())
dpops_sub <- dpops %>% filter(n>=5)
fmeta <- dmeta %>% filter(POP %in% dpops_sub$POP)
head(fmeta)
xy <- fmeta %>% dplyr::select(POP, lon, lat) %>% distinct() %>% arrange(POP)
head(xy)
rownames(xy) <- xy$POP


fmeta.W <- merge(df.pops.W, xy, by.x = "pop", by.y = "POP", sort=F)
row.names(fmeta.W) <- fmeta.W$pop
dim(fmeta.W)
head(fmeta.W)

fmeta.C <- merge(df.pops.C, xy, by.x = "pop", by.y = "POP", sort=F)
row.names(fmeta.C) <- fmeta.C$pop
dim(fmeta.C)
head(fmeta.C)

fmeta.E <- merge(df.pops.E, xy, by.x = "pop", by.y = "POP", sort=F)
row.names(fmeta.E) <- fmeta.E$pop
dim(fmeta.E)
head(fmeta.E)



#==========================#
#     Geographic data      #
#==========================#


get_mems <- function(coords) {
  
  mxy <- as.matrix(coords[c(2,3)])
  rownames(mxy) <- NULL
  nb2 <- chooseCN(coordinates(mxy), type = 6, k = 5, plot.nb = FALSE)
  listwgab <- nb2listw(nb2, style = 'W', zero.policy = TRUE)
  mem.gab <- mem(listwgab)
  row.names(mem.gab) <- rownames(coords)
  df.mem <- cbind(mxy, mem.gab)
  df.gat.mem <- df.mem[c(3:5)]
  return(df.gat.mem)
  
}



mems.W <- get_mems(fmeta.W)
head(mems.W)
dim(mems.W)

mems.C <- get_mems(fmeta.C)
head(mems.C)
dim(mems.C)

mems.E <- get_mems(fmeta.E)
head(mems.E)
dim(mems.E)



#======================#
#      PCA on AF       #
#======================#

pca.snps.W <- dudi.pca(dgeno.W, scale = FALSE, scannf = FALSE, nf = 6)
pca.W <- pca.snps.W$li[c(1,2)]
head(pca.W)
dim(pca.W)

pca.snps.C <- dudi.pca(dgeno.C, scale = FALSE, scannf = FALSE, nf = 6)
pca.C <- pca.snps.C$li[c(1,2)]
head(pca.C)
dim(pca.C)

pca.snps.E <- dudi.pca(dgeno.E, scale = FALSE, scannf = FALSE, nf = 6)
pca.E <- pca.snps.E$li[c(1,2)]
head(pca.E)
dim(pca.E)





################################################################################

#============================================#
#    Full model: genet + climate + spatial   #
#============================================#


run_full_model <- function(allele_freqs, vars, model_name) {
  
  #Running full model
  rda.mod <- rda(allele_freqs, vars)

  # The model’s explanatory power - adjusted R2
  mod.r2 <- RsquareAdj(rda.mod)$adj.r.squared

  # Model test
  mod.anova <- anova.cca(rda.mod, step = 1000)

  df.out <- data.frame("model" = c(model_name),
                     "R2" = c(mod.r2),
                     "variance" = c(mod.anova$`Variance`[1]),
                     "df" = c(mod.anova$`Df`[1]),
                     "F" = c(mod.anova$`F`[1]),
                     "P-value" = c(mod.anova$`Pr(>F)`[1]))
 return(df.out)
}


#============================================#
#    Climate: F ~ climate | genet + spatial  #
#============================================#



run_full_climate <- function(allele_freqs, vars, model_name) {
  
  #Running full model
  rda.mod <- rda(allele_freqs ~ PC1 + PC2 + PC3 + 
                   Condition(Axis1 + Axis2 + MEM1 + MEM2 + MEM3),
                 data = vars)
  
  # The model’s explanatory power - adjusted R2
  mod.r2 <- RsquareAdj(rda.mod)$adj.r.squared
  
  # Model test
  mod.anova <- anova.cca(rda.mod, step = 1000)
  
  df.out <- data.frame("model" = c(model_name),
                       "R2" = c(mod.r2),
                       "variance" = c(mod.anova$`Variance`[1]),
                       "df" = c(mod.anova$`Df`[1]),
                       "F" = c(mod.anova$`F`[1]),
                       "P-value" = c(mod.anova$`Pr(>F)`[1]))
  return(df.out)
}


run_climate <- function(allele_freqs, vars, model_name) {
  
  #Running full model
  rda.mod <- rda(allele_freqs ~ PC1 + PC2 + PC3,
                 data = vars)
  
  # The model’s explanatory power - adjusted R2
  mod.r2 <- RsquareAdj(rda.mod)$adj.r.squared
  
  # Model test
  mod.anova <- anova.cca(rda.mod, step = 1000)
  
  df.out <- data.frame("model" = c(model_name),
                       "R2" = c(mod.r2),
                       "variance" = c(mod.anova$`Variance`[1]),
                       "df" = c(mod.anova$`Df`[1]),
                       "F" = c(mod.anova$`F`[1]),
                       "P-value" = c(mod.anova$`Pr(>F)`[1]))
  return(df.out)
}

#============================================#
#    Spatial: F ~ spatial | climate + genet  #
#============================================#

run_full_spatial <- function(allele_freqs, vars, model_name) {
  
  #Running full model
  rda.mod <- rda(allele_freqs ~ MEM1 + MEM2 + MEM3 +
                   Condition(Axis1 + Axis2 + 
                               PC1 + PC2 + PC3),
                 data = vars)
  
  # The model’s explanatory power - adjusted R2
  mod.r2 <- RsquareAdj(rda.mod)$adj.r.squared
  
  # Model test
  mod.anova <- anova.cca(rda.mod, step = 1000)
  
  df.out <- data.frame("model" = c(model_name),
                       "R2" = c(mod.r2),
                       "variance" = c(mod.anova$`Variance`[1]),
                       "df" = c(mod.anova$`Df`[1]),
                       "F" = c(mod.anova$`F`[1]),
                       "P-value" = c(mod.anova$`Pr(>F)`[1]))
  return(df.out)
}



#============================================#
#    Genet : F ~ genet | climate + spatial   #
#============================================#

run_full_genet <- function(allele_freqs, vars, model_name) {
  
  #Running full model
  rda.mod <- rda(allele_freqs ~ Axis1 + Axis2 + 
                   Condition(PC1 + PC2 + PC3 +
                               MEM1 + MEM2 + MEM3),
                 data = vars)
  
  # The model’s explanatory power - adjusted R2
  mod.r2 <- RsquareAdj(rda.mod)$adj.r.squared
  
  # Model test
  mod.anova <- anova.cca(rda.mod, step = 1000)
  
  df.out <- data.frame("model" = c(model_name),
                       "R2" = c(mod.r2),
                       "variance" = c(mod.anova$`Variance`[1]),
                       "df" = c(mod.anova$`Df`[1]),
                       "F" = c(mod.anova$`F`[1]),
                       "P-value" = c(mod.anova$`Pr(>F)`[1]))
  return(df.out)
}


################################################################################


#======================#
#     Partial RDA      #
#======================#

# WEST
rownames(mems.W)
rownames(pca.W)
rownames(env.W)
variables.all <- cbind(pca.W,env.W,mems.W)
head(variables.all)

df.1 <- run_full_model(dgeno.W, variables.all, "full")
df.2 <- run_full_climate(dgeno.W, variables.all, "climate")
df.3 <- run_full_spatial(dgeno.W, variables.all, "spatial")
df.4 <- run_full_genet(dgeno.W, variables.all, "genet")
df.5 <- run_climate(dgeno.W, variables.all, "climate_only")
dm <- rbind(df.1, df.2, df.3, df.4, df.5)
dm$p.adjust <- p.adjust(dm$P.value)
head(dm)
write.table(dm, "05_partialRDA_West.tsv", sep="\t", col.names = T, row.names = F, quote=F, append=F)



# CENTRAL
rownames(mems.C)
rownames(pca.C)
rownames(env.C)
variables.all <- cbind(pca.C,env.C,mems.C)
head(variables.all)

df.1 <- run_full_model(dgeno.C, variables.all, "full")
df.2 <- run_full_climate(dgeno.C, variables.all, "climate")
df.3 <- run_full_spatial(dgeno.C, variables.all, "spatial")
df.4 <- run_full_genet(dgeno.C, variables.all, "genet")
df.5 <- run_climate(dgeno.C, variables.all, "climate_only")
dm <- rbind(df.1, df.2, df.3, df.4, df.5)
dm$p.adjust <- p.adjust(dm$P.value)
head(dm)
write.table(dm, "05_partialRDA_Central.tsv", sep="\t", col.names = T, row.names = F, quote=F, append=F)


# EAST
rownames(mems.E)
rownames(pca.E)
rownames(env.E)
variables.all <- cbind(pca.E,env.E,mems.E)
head(variables.all)

df.1 <- run_full_model(dgeno.E, variables.all, "full")
df.2 <- run_full_climate(dgeno.E, variables.all, "climate")
df.3 <- run_full_spatial(dgeno.E, variables.all, "spatial")
df.4 <- run_full_genet(dgeno.E, variables.all, "genet")
df.5 <- run_climate(dgeno.E, variables.all, "climate_only")
dm <- rbind(df.1, df.2, df.3, df.4, df.5)
dm$p.adjust <- p.adjust(dm$P.value)
head(dm)
write.table(dm, "05_partialRDA_East.tsv", sep="\t", col.names = T, row.names = F, quote=F, append=F)






### Venn diagram of contributions of each class
head(sel.env)
vp1 <- varpart(dgeno.W, pca.W, env.W, mems.W)
vp1


png("01_partialRDA_Venn.png",w=700,h=700,res=150)
plot(vp1, Xnames = c("genetic","climate", "spatial"))
dev.off()



