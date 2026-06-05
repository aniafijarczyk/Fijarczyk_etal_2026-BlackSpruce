rm(list=ls())
setwd("C:/Users/aniaf/Projects/BlackSpruce/34_climate_transfer_distance")


library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape)
library(gridExtra)
library(grid)
library(RColorBrewer)






########################
###---    DATA    ---###
########################


# Getting list of populations
dmeta <- read.csv("../31_climate/02_garden_climate_with_gardens.tab", sep="\t", header=TRUE)
head(dmeta)
dim(dmeta)
rownames(dmeta) <- dmeta$KeyID

dclim <- dmeta %>% dplyr::select(-KeyID, -Latitude, -Longitude, -Elevation, -dataset, -POP_GR)
head(dclim)
dim(dclim)

###################
#----  PCA  ------#
###################

pca.clim <- prcomp(dclim, scale=TRUE, center=TRUE, retx=TRUE)
pca.clim$rotation
pca.clim$x
biplot(pca.clim)
var <- pca.clim$sdev^2
prop_var <- var / sum(var)
sum(prop_var[1:4])
# First four PCs explain ~86% of variance

#############################
###   CLIMATE DISTANCES   ###
#############################

# PCA data
pca.run <- pca.clim
pca.x <- as.data.frame(pca.run$x)
head(pca.x)
pca.x.meta <- merge(pca.x, dmeta[c("KeyID","Latitude","Longitude","dataset","POP_GR")], by.x = 0, by.y = 'KeyID', sort=FALSE)
head(pca.x.meta)
#write.table(pca.x.meta, "03_climate_distance_PCA.tsv", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, append=FALSE)

r0 <- ggplot(pca.x.meta) + aes(x = PC1, y = PC2) +
  geom_point(aes(color = as.factor(dataset))) +
  #scale_fill_manual(values = c(brewer.pal(9, "BrBG"))) +
  labs(x = "PC1 climate", y = "PC2 climate") +
  theme(panel.background = element_rect(fill=NA, color="black"),
        panel.grid = element_blank(),
        legend.position="none",
        axis.text = element_text(size=14),
        axis.title = element_text(size=18))
r0


euc_dist <- function(x1, x2){
  return(sqrt(sum((x1 - x2)^2)))
}

x <- c(pca.x.meta[1,2], pca.x.meta[1,3])
y <- c(pca.x.meta[2,2], pca.x.meta[2,3])
euc_dist(x,y)


pca.x.meta$Row.names
common_gardens <- c("Chibougamau","Valcartier","MontLaurier","Acadia","PeaceRiver")
provenances <- pca.x.meta %>% filter(!Row.names %in% common_gardens) %>% pull(Row.names)
length(provenances)

Mat <- matrix(ncol = 3)
Mat
for (cg in 1:length(common_gardens)) {
  common_garden.1 <- pca.x.meta %>% filter(Row.names == common_gardens[cg]) %>% dplyr::select(PC1) %>% as.numeric()
  common_garden.2 <- pca.x.meta %>% filter(Row.names == common_gardens[cg]) %>% dplyr::select(PC2) %>% as.numeric()
  common_garden.3 <- pca.x.meta %>% filter(Row.names == common_gardens[cg]) %>% dplyr::select(PC3) %>% as.numeric()
  common_garden.4 <- pca.x.meta %>% filter(Row.names == common_gardens[cg]) %>% dplyr::select(PC4) %>% as.numeric()
  
  common_garden.points <- c(common_garden.1, common_garden.2, common_garden.3, common_garden.4)
  print(common_garden.points)
  
  mat <- matrix(nrow = length(provenances), ncol = 3)
  for (prov in 1:length(provenances)) {
    prov.1 <- pca.x.meta %>% filter(Row.names == provenances[prov]) %>% dplyr::select(PC1) %>% as.numeric()
    prov.2 <- pca.x.meta %>% filter(Row.names == provenances[prov]) %>% dplyr::select(PC2) %>% as.numeric()
    prov.3 <- pca.x.meta %>% filter(Row.names == provenances[prov]) %>% dplyr::select(PC3) %>% as.numeric()
    prov.4 <- pca.x.meta %>% filter(Row.names == provenances[prov]) %>% dplyr::select(PC4) %>% as.numeric()
    
    dist <- euc_dist(c(common_garden.1, common_garden.2, common_garden.3, common_garden.4), c(prov.1, prov.2, prov.3, prov.4))
    mat[prov, 1] <- common_gardens[cg]
    mat[prov, 2] <- provenances[prov]
    mat[prov, 3] <- dist
  }
  Mat <- rbind(Mat, mat)
}

Mat
dM <- as.data.frame(Mat) %>% drop_na()
colnames(dM) <- c('CommonGarden','Provenance','EuclDist')
dM %>% head()




### Formatting table

dM %>% head()
dim(dM)

clim_dist <- dM %>% spread(CommonGarden, EuclDist)
clim_dist <- clim_dist %>% dplyr::rename(ECD_AC = Acadia,
                            ECD_CH = Chibougamau,
                            ECD_ML = MontLaurier,
                            ECD_PR = PeaceRiver,
                            POP_ID = Provenance)
head(clim_dist)


write.table(clim_dist, "03_climate_distance.tsv", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, append=FALSE)



