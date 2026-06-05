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
library(dartR)
library(vegan)
library(sp)
library(stringr)
library(sf)
library(raster)
library(spData)
library(tmap)



#===================#
#       DATA        #
#===================#

# Importing population AFs 
dgeno <- read.csv("../26_MoransI/01_moran_Pop_AFs.csv", sep=",", header=T)
dgeno[c(1:5), c(1:5)]
rownames(dgeno)
dim(dgeno)


# Metadata with geographic info
dmeta <- read.csv("../DATA_intermediate/23_filter_EPN_indiv_metrics.tsv", sep="\t", header=T)
dpops <- dmeta %>% filter(SPECIES_ID == 'EPN') %>% group_by(POP) %>% dplyr::summarise(n=n())
hist(dpops$n)
dpops_sub <- dpops %>% filter(n>=5)
fmeta <- dmeta %>% filter(POP %in% dpops_sub$POP)
dim(fmeta)
head(fmeta)



#============================#
# Spatial connectivity model #
#============================#

# Samples list
xy <- fmeta %>% dplyr::select(POP, lon, lat) %>% distinct() %>% arrange(POP)
head(xy)
rownames(xy) <- xy$POP
xy[c(2,3)]
mxy <- as.matrix(xy[c(2,3)])
head(mxy)
dim(mxy)
rownames(mxy) <- NULL
dim(mxy)
head(mxy)


#listw.explore()
nb2 <- chooseCN(coordinates(mxy), type = 6, k = 8, plot.nb = FALSE)
s.Spatial(mxy, nb = nb2, plabel.cex = 0, pnb.edge.col = 'grey')
listwgab <- nb2listw(nb2, style = 'W', zero.policy = TRUE)
listwgab


#===================#
#       sPCA        #
#===================#

# Here we run spatial autocorrelation on pop allelic frequencies of randomly selected SNPs;
# in contrast earlier it was done on PCs derived from pop allelic frequencies

# Get coordinates
head(xy)
xy_coords <- xy %>% dplyr::select(lon, lat)

# Selecting radom 5000 snps
set.seed(101)
dgeno_sub <- sample(dgeno, 5000)
dim(dgeno_sub)

# Interactive run showing the screeplot to select the number of global and local axes
spca(dgeno_sub, xy=xy_coords, cn=listwgab)


# Run sPCA
dgeno.spca1 <- spca(dgeno_sub, xy=xy_coords, cn=listwgab, scannf=FALSE,
                     nfposi=4, nfnega=0)
dgeno.spca1


# Saving output
save(dgeno.spca1, file = "03_spca.Rdata")




#===================#
#      Plots        #
#===================#

# Import data if necessary
load("03_spca.Rdata")

# Plotting eigenvalues
png("03_spca_EIG_barplot.png",w=1000,h=800,res=150)
barplot(dgeno.spca1$eig, col=rep(c("red","grey"), c(4,1000)),
        main="BS dataset - sPCA eigenvalues")
dev.off()


# Screeplot - takes time
#screeplot(dgeno.spca1)

# Mapping the sPCA results using s.value and lagged scores ($ls) instead of the PC ($li),
# which are a 'denoisied' version of the PCs.
dgeno.spca1$xy
s.value(dgeno.spca1$xy, dgeno.spca1$ls[,1])
s.value(dgeno.spca1$xy, dgeno.spca1$ls[,2])
s.value(dgeno.spca1$xy, dgeno.spca1$ls[,3])
s.value(dgeno.spca1$xy, dgeno.spca1$ls[,4])


# Calculating Moran's I
MC.snps <- moran.randtest(dgeno.spca1$ls, listwgab, nrepet = 9999)
MC.snps

df.MC <- data.frame("Axis" = MC.snps$names,
                    "Obs" = as.numeric(MC.snps$obs),
                    "Alter" = MC.snps$alter,
                    "Pvalue" = MC.snps$adj.pvalue,
                    "nrep" = as.numeric(MC.snps$rep))
df.MC

# Getting table
df.spca <- dgeno.spca1$ls
df.spca[c("x","y")] <- dgeno.spca1$xy
head(df.spca)

dg.spca <- df.spca %>% gather(key = "Axis", value = "value", c(1:4))
dg.spca$fill <- ifelse(dg.spca$value <=0, "N", "P")
dg.spca$abs_value <- abs(dg.spca$value)
head(dg.spca)

p4 <- ggplot(dg.spca) + aes(x = x, y = y, size = abs_value, fill = fill) + 
  geom_point(pch=21) +
  scale_fill_manual(values = c("black","khaki")) +
  facet_wrap(~Axis, ncol=2) +
  theme(panel.background = element_rect(fill=NA, colour="black"),
        panel.grid = element_blank(),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=12, face = "bold"))
p4




colorplot(dgeno.spca1$xy, dgeno.spca1$ls, axes=1:3, cex=3)

#===================#
#    Nice plot      #
#===================#

world_NA = world[world$continent == "North America", ]
my_sf <- read_sf("C:/Users/aniaf/Projects/BlackSpruce/METADATA/maps/BlackSpruce/data/commondata/data0/picemari.shp")
species_distribution_outline <- st_union(my_sf)

# Moran I values
df.MC$lon <- -170
df.MC$lat <- 42
df.MC

p5 <- ggplot() + 
  geom_sf(data = world_NA, fill = "floralwhite") +
  geom_sf(data = species_distribution_outline, fill = "white") +
  geom_point(data = dg.spca, aes(x = x, y = y, size = abs_value, fill = fill), pch=21) +
  
  geom_text(data = df.MC, aes(x = lon, y = lat, label = paste0("Moran's I=", round(Obs,2),"; p=",round(Pvalue,4))),
            color="black", hjust=0) +
  
  scale_fill_manual(values = c("black","khaki")) +
  scale_size_continuous(name="Absolute\nvalue") +
  scale_x_continuous(limits = c(-170, -52)) +
  scale_y_continuous(limits = c(40, 70)) +
  facet_wrap(~Axis, ncol=2) +
  labs(x = "Longitude", y = "Latitude") +
  guides(fill = "none") +
  theme(panel.background = element_rect(fill=NA, colour="black"),
        panel.grid = element_blank(),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=12, face = "bold"))
p5



png("03_spca_4Axes.png",w=2400,h=1400,res=300)
p5
dev.off()






