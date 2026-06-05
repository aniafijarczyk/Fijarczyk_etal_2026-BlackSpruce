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





#===================#
#   Genetic data    #
#===================#

# AFs for populations 
dgeno <- read.csv("../926_MoransI/01_moran_Pop_AFs.csv", sep=",", header=T)
dgeno[c(1:5), c(1:5)]
rownames(dgeno)
dim(dgeno)
pops <- rownames(dgeno)
pops
df.pops <- data.frame("pop" = pops)
head(df.pops)



#======================#
#  Environmental data  #
#======================#

# Data
df.env <- read.csv("../31_climate/02_garden_climate_PCA.tsv", sep="\t", header=TRUE)
head(df.env)

# Filter populations
env.pops <- merge(df.pops, df.env, by.x = "pop", by.y = "Row.names", sort=F)
dim(env.pops)
head(env.pops)
row.names(env.pops) <- env.pops$pop

env <- env.pops %>% dplyr::select(PC1, PC2, PC3)
head(env)
dim(env)




#======================#
#       Top MEMs       #
#======================#

# Data
mems <- read.csv("../926_MoransI/01_moran_MEMs.tsv", sep="\t", header=TRUE)
head(mems)
mems$pop <- rownames(mems)
df.mems <- merge(df.pops, mems, by = "pop", sort=F)
head(df.mems)
rownames(df.mems) <- df.mems$pop
mems.3 <- df.mems %>% dplyr::select(MEM1, MEM2, MEM3)
head(mems.3)
dim(mems.3)




#======================#
#      PCA on AF       #
#======================#

pca.snps <- dudi.pca(dgeno, scale = FALSE, scannf = FALSE, nf = 6)
pca.2.axes <- pca.snps$li[c(1,2)]
head(pca.2.axes)



#======================#
#     Partial RDA      #
#======================#

rownames(mems.3)
rownames(pca.2.axes)
rownames(env)


# Combining datasets
cols.gen <- colnames(pca.2.axes)
cols.env <- colnames(env)
cols.mem <- colnames(mems.3)
variables.all <- cbind(pca.2.axes,env,mems.3)
head(variables.all)



#============================================#
#    Full model: genet + climate + spatial   #
#============================================#

#Running full model
spe.partial.rda.full <- rda(dgeno, variables.all)
summary(spe.partial.rda.full)

# The model’s explanatory power - adjusted R2
mod.r2 <- RsquareAdj(spe.partial.rda.full)$adj.r.squared
mod.r2

# Model test
mod.anova <- anova.cca(spe.partial.rda.full, step = 1000)
mod.anova

df.1 <- data.frame("model" = c("full"),
                   "R2" = c(mod.r2),
                   "variance" = c(mod.anova$`Variance`[1]),
                   "df" = c(mod.anova$`Df`[1]),
                   "F" = c(mod.anova$`F`[1]),
                   "P-value" = c(mod.anova$`Pr(>F)`[1]))
df.1

#============================================#
#    Climate: F ~ climate | genet + spatial  #
#============================================#

spe.partial.rda.climate <- rda(dgeno ~ 
                                 PC1 + PC2 + PC3 +
                                 Condition(Axis1 + Axis2 + MEM1 + MEM2 + MEM3),
                               data = variables.all)
spe.partial.rda.climate

mod.r2 <- RsquareAdj(spe.partial.rda.climate)$adj.r.squared
mod.r2

mod.anova <- anova.cca(spe.partial.rda.climate, step = 1000)
mod.anova

df.2 <- data.frame("model" = c("climate"),
                   "R2" = c(mod.r2),
                   "variance" = c(mod.anova$`Variance`[1]),
                   "df" = c(mod.anova$`Df`[1]),
                   "F" = c(mod.anova$`F`[1]),
                   "P-value" = c(mod.anova$`Pr(>F)`[1]))



#============================================#
#    Spatial: F ~ spatial | climate + genet  #
#============================================#

spe.partial.rda.geo <- rda(dgeno ~ MEM1 + MEM2 + MEM3 +
                             Condition(Axis1 + Axis2 + 
                                         PC1 + PC2 + PC3),
                           data = variables.all)
spe.partial.rda.geo

mod.r2 <- RsquareAdj(spe.partial.rda.geo)$adj.r.squared
mod.r2

mod.anova <- anova.cca(spe.partial.rda.geo, step = 1000)
mod.anova

df.3 <- data.frame("model" = c("spatial"),
                   "R2" = c(mod.r2),
                   "variance" = c(mod.anova$`Variance`[1]),
                   "df" = c(mod.anova$`Df`[1]),
                   "F" = c(mod.anova$`F`[1]),
                   "P-value" = c(mod.anova$`Pr(>F)`[1]))



#============================================#
#    Genet : F ~ genet | climate + spatial   #
#============================================#

spe.partial.rda.structure <- rda(dgeno ~ Axis1 + Axis2 + 
                                   Condition(PC1 + PC2 + PC3 +
                                               MEM1 + MEM2 + MEM3),
                                 data = variables.all)
summary(spe.partial.rda.structure)

mod.r2 <- RsquareAdj(spe.partial.rda.structure)$adj.r.squared
mod.r2

mod.anova <- anova.cca(spe.partial.rda.structure, step = 1000)
mod.anova


df.4 <- data.frame("model" = c("genet"),
                   "R2" = c(mod.r2),
                   "variance" = c(mod.anova$`Variance`[1]),
                   "df" = c(mod.anova$`Df`[1]),
                   "F" = c(mod.anova$`F`[1]),
                   "P-value" = c(mod.anova$`Pr(>F)`[1]))


dm <- rbind(df.1, df.2, df.3, df.4)
head(dm)
write.table(dm, "01_partialRDA.tsv", sep="\t", col.names = T, row.names = F, quote=F, append=F)




### Venn diagram of contributions of each class
head(sel.env)
vp1 <- varpart(dgeno, pca.2.axes, env, mems.3)
vp1


png("01_partialRDA_Venn.png",w=700,h=700,res=150)
plot(vp1, Xnames = c("genetic","climate", "spatial"))
dev.off()



