rm(list=ls())
setwd("C:/Users/aniaf/Projects/BlackSpruce/930_RDA")

library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape)
library(gridExtra)
library(grid)





#===================#
#  Calculating AF   #
#===================#

# Getting per population allele frequencies
load("../DATA_intermediate/23_filter_EPN_genlight.Rdata")
genl@other$group

popNames(genl)
genl@other %>% head()
as.character(genl@other$group)

gel.clust <- genl[which(genl@other$group == "ME"), ]
gel.clust <- gl.filter.monomorphs(gel.clust)
gel.clust <- gl.recalc.metrics(gel.clust)
pop.names <- popNames(gel.clust)
loc.names <- locNames(gel.clust)
pop.names
pop = "4360"
gel.pop <- gl.keep.pop(gel.clust, pop, as.pop="POP")
gel.pop

# A function to calculate allele frequency of allele 1 for each population
# Columns (SNPs) with NaNs are removed


calculate_AF1 <- function(x) {
  CPOPS = list()
  cluster.names <- unique(as.character(x@other$group))
  
  for (clust in cluster.names) {
    print(clust)
    gel.clust <- x[which(x@other$group == clust), ]
    gel.clust <- gl.filter.monomorphs(gel.clust)
    gel.clust <- gl.recalc.metrics(gel.clust)
    pop.names <- popNames(gel.clust)
    loc.names <- locNames(gel.clust)
    
    DF <- vector()
    pop.names.list <- vector()
    for (pop in pop.names) {
      print(pop)
      gel.pop <- gl.keep.pop(gel.clust, pop, as.pop="POP")
      
      if (length(gel.pop$ind.names)>=5) {
        af.pop <- gl.alf(gel.pop)$alf1
        DF <- cbind(DF, af.pop)
        pop.names.list[pop] <- pop
      } else {
        DF <- DF
      }
    }
    
    colnames(DF) <- pop.names.list
    rownames(DF) <- locNames(gel.clust)
    af.df <- as.data.frame(t(DF))
    af.df.nonan <- af.df[, colSums(is.na(af.df)) == 0, drop = FALSE]
    
    CPOPS[[clust]] <- af.df.nonan
  }

  return(CPOPS)
}


###

# Calculate AF1
df.genet <- calculate_AF1(genl)
names(df.genet)
df.genet$Central
df.genet[c(1:10), c(1:10)]
dim(df.genet)
df.genet$ME
cluster <- "ME"
dsub <- df.genet[[cluster]]
dsub[1:5, 1:5]


for (cluster in names(df.genet)) {
  # Write down table
  print(cluster)
  dsub <- df.genet[[cluster]]
  write.table(dsub, file = paste0("04_AF_by_cluster_",cluster,".csv"), sep=",", col.names = TRUE, row.names = TRUE, append=FALSE)
}




