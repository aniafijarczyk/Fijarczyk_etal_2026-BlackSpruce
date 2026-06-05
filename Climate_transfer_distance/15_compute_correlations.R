rm(list=ls())
setwd("C:/Users/aniaf/Projects/BlackSpruce/943_climate_transfer_distance")


library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape)
library(cowplot)
library(grid)
library(RColorBrewer)
library(SpatialPack) # modified.ttest


# Calculating correlation between phenotype and climate transfer distance
# Accounting for autocorrelation

# Performs a modified version of the t test to assess the correlation between two spatial processes.
#  Spearman correlation
# Dutilleul test

#===========#
#   DATA    #
#===========#

dclim <- read.csv("./13_plot_climate_distance_input.tsv", sep="\t", header=T)
dclim %>% head()
class(dclim$POP_ID)

# Getting coordinates
dcoords <- read.csv("../METADATA/analyzeTables_pop_summary.tsv", sep="\t", header=T)
dcoords <- dcoords %>% dplyr::select(-n)
class(dcoords$pop)

df <- merge(dclim, dcoords, by.x = "POP_ID", by.y = "pop", sort=F, all.x=T)
head(df)
df[is.na(df$EuclDist),]


#===================#
#   One example     #
#===================#

head(df)
dx <- df %>% filter(SITE_ID == "AC") %>% filter(Trait == "Height")
dim(dx)
coords <- as.matrix(dx[c("lat","lon")])
coords <- dx[c("lat","lon")]
coords <- dx[c("lon","lat")]
xr <- rank(dx$mean)    # Spearman => use ranks
yr <- rank(dx$EuclDist)
cor.test(xr, yr)

dtest <- modified.ttest(xr, yr, coords, nclass = NULL)
dtest
dtest$corr
dtest$dof
dtest$p.value
dtest$upper.bounds
dtest$Fstat
dtest$ESS
dtest$dims

#===================#
#     FUNCTION      #
#===================#

calculate_dtest <- function(dataset, X, Y) {
  xr <- rank(dataset[[X]])
  yr <- rank(dataset[[Y]])
  coords <- dataset[c("lat","lon")] 
  dtest <- modified.ttest(xr, yr, coords, nclass = NULL)
  output <- c(dtest$corr, dtest$p.value, dtest$Fstat, dtest$dof, length(xr))
  return(output)
}


out <- calculate_dtest(dx, "mean", "EuclDist")
out


#===================#
#       LOOP        #
#===================#

traits <- unique(df$Trait)
traits
gardens <- unique(df$SITE_ID)
gardens

results <- list()
i <- 1
for (trait in traits) {
  for (garden in gardens) {
    
    df_ <- df %>% filter(Trait == trait) %>% filter(SITE_ID == garden)
    cor_test <- calculate_dtest(df_, "mean", "EuclDist")
    cor_test[6] <- trait
    cor_test[7] <- garden
    results[[i]] <- cor_test
    i = i + 1
  }
}
results


dR <- as.data.frame(do.call(rbind, results))
colnames(dR) <- c("rho","p.value","Fstat","dof","n_trees","Trait","SITE_ID")
dR$adj.p <- p.adjust(dR$p.value, method = "BH")
head(dR)

head(df)
dM <- merge(df, dR, by=c("Trait","SITE_ID"), sort=F, all.x=T)
head(dM)

### SAVE

write.table(dM, file = "15_compute_correlations.tsv", sep="\t", col.names = T, row.names = F, quote=F, append=F)




