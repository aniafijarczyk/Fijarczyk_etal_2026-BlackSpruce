rm(list=ls())
setwd("C:/Users/aniaf/Projects/BlackSpruce/926_MoransI")


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
library(vegan)
library(sp)
library(stringr)




#===================#
#     Metadata      #
#===================#

# Selecting pops
fmeta <- read.csv("../31_climate/02_garden_climate_with_gardens.tab", sep="\t", header=T)
head(fmeta)
dim(fmeta)
class(fmeta$KeyID)

pops <- fmeta %>% arrange(KeyID) %>% pull(KeyID) %>% unique()
pops
length(pops)



#==========================#
#     Geographic data      #
#==========================#

# Samples list
xy <- fmeta %>% dplyr::select(KeyID, Longitude, Latitude) %>% distinct() %>% arrange(KeyID)
rownames(xy) <- xy$KeyID
xy[c(2,3)]
mxy <- as.matrix(xy[c(2,3)])
head(mxy)
dim(mxy)

rownames(mxy) <- NULL
dim(mxy)
head(mxy)
s.label(mxy, ppoint.pch = 20, ppoint.col = "darkseagreen4")
ggplot(xy) + aes(x = Longitude, y = Latitude, label = KeyID) + geom_point() + geom_text()



#===================#
#     NB & SWM      #
#===================#

# Spatial neighborhoods (NB) are connectivity matrices which can be represented 
# by an unweighted graph
# A spatial weighting matrices (SWM) is computed by a transformation of 
# a spatial neighborhood

listw.explore()

# selecting neighbor matrix based on k-neighbors
nb2 <- chooseCN(coordinates(mxy), type = 6, k = 8, plot.nb = FALSE)
s.Spatial(mxy, nb = nb2, plabel.cex = 0, pnb.edge.col = 'red')
listwgab <- nb2listw(nb2, style = 'W', zero.policy = TRUE)
listwgab

s.Spatial(mxy, nb = nb2, plabel.cex = 0, pnb.edge.col = 'red')




#===================#
#       MEMs        #
#===================#

# MEMs - spatial predictors (eigenfunctions)
# orthogonal vectors with a unit norm that maximize Moran’s coefficient 
# of spatial autocorrelation


row.names(xy)

mem.gab <- mem(listwgab)
row.names(mem.gab)
row.names(xy)

row.names(mem.gab) <- rownames(xy)
head(mem.gab)
attr(mem.gab, "values")
attr(mem.gab, "weights")


plot(mem.gab[,c(1, 2, 3, 4, 5, 6, 7)], SpORcoords = mxy)
s.value(mxy, mem.gab[,c(1, 2, 3, 4, 5, 6, 7)], symbol = "circle", ppoint.cex = 0.6)

# Saving MEMs
write.table(mem.gab, file = "02_moran_with_gardens_MEMs.tsv", sep="\t", col.names = TRUE, row.names = TRUE, append=FALSE)

