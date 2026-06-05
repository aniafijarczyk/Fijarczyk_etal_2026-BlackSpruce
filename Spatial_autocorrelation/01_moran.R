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
dmeta <- read.csv("../DATA_intermediate/23_filter_EPN_indiv_metrics.tsv", sep="\t", header=T)
head(dmeta)
dim(dmeta)
class(dmeta$POP)

dpops <- dmeta %>% filter(SPECIES_ID == 'EPN') %>% group_by(POP) %>% dplyr::summarise(n=n())
hist(dpops$n)
dpops_sub <- dpops %>% filter(n>=5)
fmeta <- dmeta %>% filter(POP %in% dpops_sub$POP)
dim(fmeta)
head(fmeta)

# Sorted list of pops (numerically)
pops <- fmeta %>% arrange(POP) %>% pull(POP) %>% unique()
pops
length(pops)




#==========================#
#     Geographic data      #
#==========================#

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
s.label(mxy, ppoint.pch = 20, ppoint.col = "darkseagreen4")
ggplot(xy) + aes(x = lon, y = lat, label = POP) + geom_point() + geom_text()



#===================#
#     NB & SWM      #
#===================#

# Spatial neighborhoods (NB) are connectivity matrices which can be represented 
# by an unweighted graph
# A spatial weighting matrices (SWM) is computed by a transformation of 
# a spatial neighborhood

#listw.explore()

# Network based on distance is tricky; when increasing the distance (d2), it creates
# a very dense mesh, and when decreasing it (d2) some locations are left as islands

# This is a neighbor matrix based on distance network
nb <- chooseCN(coordinates(mxy), type = 5, d1 = 0, d2 = 13, plot.nb = FALSE)
distnb <- nbdists(nb, mxy)
fdist <- lapply(distnb, function(x) 1 - x/max(dist(mxy)))
listwgab <- nb2listw(nb, style = 'W', glist = fdist, zero.policy = TRUE)
s.Spatial(mxy, nb = nb, plabel.cex = 0, pnb.edge.col = 'black')

# Checking other option with smaller distance
nb <- chooseCN(coordinates(mxy), type = 5, d1 = 0, d2 = 3, plot.nb = FALSE)
s.Spatial(mxy, nb = nb, plabel.cex = 0, pnb.edge.col = 'black')

nb <- chooseCN(coordinates(mxy), type = 5, d1 = 0, d2 = 8, plot.nb = FALSE)
s.Spatial(mxy, nb = nb, plabel.cex = 0, pnb.edge.col = 'black')


# This is network based on K-neighbors
nb2 <- chooseCN(coordinates(mxy), type = 6, k = 5, plot.nb = FALSE)
s.Spatial(mxy, nb = nb2, plabel.cex = 0, pnb.edge.col = 'red')

nb2 <- chooseCN(coordinates(mxy), type = 6, k = 6, plot.nb = FALSE)
s.Spatial(mxy, nb = nb2, plabel.cex = 0, pnb.edge.col = 'red')

nb2 <- chooseCN(coordinates(mxy), type = 6, k = 7, plot.nb = FALSE)
s.Spatial(mxy, nb = nb2, plabel.cex = 0, pnb.edge.col = 'red')

nb2 <- chooseCN(coordinates(mxy), type = 6, k = 8, plot.nb = FALSE)
s.Spatial(mxy, nb = nb2, plabel.cex = 0, pnb.edge.col = 'red')


# Creating swm
nb2 <- chooseCN(coordinates(mxy), type = 6, k = 8, plot.nb = FALSE)
s.Spatial(mxy, nb = nb2, plabel.cex = 0, pnb.edge.col = 'red')
listwgab <- nb2listw(nb2, style = 'W', zero.policy = TRUE)
listwgab

# Saving network plot
png("01_moran_network_distance.png",w=1000,h=800,res=150)
s.Spatial(mxy, nb = nb2, plabel.cex = 0, pnb.edge.col = 'red')
dev.off()


#===================#
#       MEMs        #
#===================#

# MEMs - spatial predictors (eigenfunctions)
# orthogonal vectors with a unit norm that maximize Moran’s coefficient 
# of spatial autocorrelation

mem.gab <- mem(listwgab)
head(mem.gab)
row.names(mem.gab) <- rownames(xy)
head(mem.gab)
attr(mem.gab, "values")
attr(mem.gab, "weights")

png("01_moran_MEMs_barplot.png",w=1000,h=800,res=150)
barplot(attr(mem.gab, "values"), 
        main = "Eigenvalues of the spatial weighting matrix", cex.main = 0.7)
dev.off()

# Calculating Moran's I - only global (positive) mems are significant
moranI <- moran.randtest(mem.gab, listwgab, nrepet=9999)
moranI

df.mc <- data.frame("MEM" = moranI$names,
                    "Obs" = as.numeric(moranI$obs),
                    "Alter" = moranI$alter,
                    "Pvalue" = moranI$adj.pvalue)

# Checking which mems are significant, only global (positive) or local too?
df.mc %>% filter(Pvalue < 0.05)

plot(mem.gab[,c(1:8)], SpORcoords = mxy)
s.value(mxy, mem.gab[,c(1:8)], symbol = "circle", ppoint.cex = 0.6)

png("01_moran_MEMs.png",w=1600,h=1400,res=150)
s.value(mxy, mem.gab[,c(1:6)], symbol = "circle", ppoint.cex = 0.6)
dev.off()


# Saving MEMs
write.table(mem.gab, file = "01_moran_MEMs.tsv", sep="\t", col.names = TRUE, row.names = TRUE, append=FALSE)


#===================#
#    MEMs plots     #
#===================#

### Plotting four main spatial autocorrelation axes
df.mem <- cbind(mxy, mem.gab)
head(df.mem)
df.gat.mem <- df.mem[c(1:8)] %>% gather(key = "MEM", value = "value", c(3:8))
df.gat.mem$fill <- ifelse(df.gat.mem$value <=0, "N", "P")
df.gat.mem$abs_value <- abs(df.gat.mem$value)
head(df.gat.mem)


world_NA = world[world$continent == "North America", ]
my_sf <- read_sf("C:/Users/aniaf/Projects/BlackSpruce/METADATA/maps/BlackSpruce/data/commondata/data0/picemari.shp")
species_distribution_outline <- st_union(my_sf)


df.gat.mem.sub <- df.gat.mem %>% filter(MEM %in% c("MEM1","MEM2","MEM3","MEM4"))

# morans coefficients table
df.mc$lon <- -170
df.mc$lat <- 42
df.mc <- df.mc %>% filter(MEM %in% c("MEM1","MEM2","MEM3","MEM4"))
head(df.mc)


p3 <- ggplot() + 
  geom_sf(data = world_NA, fill = "floralwhite") +
  geom_sf(data = species_distribution_outline, fill = "white") +
  geom_point(data = df.gat.mem.sub, aes(x = lon, y = lat, size = abs_value, fill = fill), pch=21) +
  
  geom_text(data = df.mc, aes(x = lon, y = lat, label = paste0("Moran's I=", round(Obs,2),"; p=",round(Pvalue,4))),
            color="black", hjust=0) +
  
  scale_fill_manual(values = c("black","khaki")) +
  scale_size_continuous(name="Absolute\nvalue") +
  scale_x_continuous(limits = c(-170, -52)) +
  scale_y_continuous(limits = c(40, 70)) +
  facet_wrap(~MEM, ncol=2) +
  labs(x = "Longitude", y = "Latitude") +
  guides(fill = "none") +
  theme(panel.background = element_rect(fill=NA, colour="black"),
        panel.grid = element_blank(),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=12, face = "bold"))
p3

png("01_moran_MEMs_4.png",w=2400,h=1400,res=300)
p3
dev.off()




#==========================#
#    Genetic structure     #
#==========================#

#===================#
#  Calculating AF   #
#===================#

# Getting per population allele frequencies
gel <- gl.load("../DATA_intermediate/23_filter_EPN.Rdata")
gel
length(pops)
as.character(pops)
popNames(gel)

af.pop <- gl.alf(gel)
head(af.pop)

# A function to calculate allele frequency of allele 1 for each population
# Columns (SNPs) with NaNs are removed
calculate_AF1 <- function(x) {
  pop.names <- as.character(pops)
  loc.names <- locNames(x)
  DF <- vector()
  for (pop in pop.names) {
    print(pop)
    #gel.pop <- x[popNames(x) == pop]
    gel.pop <- gl.keep.pop(x, pop, as.pop="POP")
    af.pop <- gl.alf(gel.pop)$alf1
    DF <- cbind(DF, af.pop)
  }
  colnames(DF) <- pop.names
  rownames(DF) <- loc.names
  af.df <- as.data.frame(t(DF))
  # Columns with NaN 
  #af.df.nonan <- af.df %>% mutate_all(~ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x))
  af.df.nonan <- af.df %>% dplyr::select_if(~ !any(is.na(.)))
  return(af.df.nonan)
}


# Calculate AF1
df.genet <- calculate_AF1(gel)
df.genet[c(1:10), c(1:10)]
dim(df.genet)

# Write down table
write.table(df.genet, file = "01_moran_Pop_AFs.csv", sep=",", col.names = TRUE, row.names = TRUE, append=FALSE)




#===================#
#        PCA        #
#===================#
# Calculating PCA based on population allele frequencies

# Import pre-calculated AFs
df.genet <- read.csv("./01_moran_Pop_AFs.csv", sep=",", header=T)
df.genet[c(1:5), c(1:5)]
rownames(df.genet)
dim(df.genet)

# Run PCa with dudi.pca
pca.snps <- dudi.pca(df.genet, scale = FALSE, scannf = FALSE, nf = 6)
dim(pca.snps$li)
head(pca.snps$li)
barplot(pca.snps$eig)
pca.snps$tab[1:10,1:10]
class(pca.snps)

# Check PCA plot if it's what you expect
plot(pca.snps$li[,1], pca.snps$li[,2])
pca.eiv <- pca.snps$li[,1:4]
pca.eiv$POP <- rownames(pca.eiv)
pop.info <- fmeta %>% dplyr::select(POP, lat, lon, STRATA) %>% distinct()
pca.out <- merge(pca.eiv, pop.info, by = "POP", sort=F)
head(pca.out)

ggplot(pca.out) + aes(x = Axis1, y = Axis2, fill=STRATA) +
  scale_fill_manual(values = brewer.pal(9, "BrBG")) +
  geom_point(size=3, pch=21)


### Checking if rownames order is correct
row.names(df.genet)
row.names(pca.snps$li)
row.names(mem.gab)





#==========================#
# Spatial patterns of AF   #
#==========================#

# Spatial autocorrelations associated with pop structure axes
# (spatial autocorrelation of allele frequencies)
# A centered principal component analysis on the snp data

# MC (Morans Coefficient) can be computed for PCA scores and the associated 
# spatial structures can be visualized on a map
listwgab

pc1.mctest <- moran.mc(pca.snps$li[,1], listwgab, 999)
plot(pc1.mctest)

pc2.mctest <- moran.mc(pca.snps$li[,2], listwgab, 999)
plot(pc2.mctest)

pc3.mctest <- moran.mc(pca.snps$li[,3], listwgab, 999)
plot(pc3.mctest)

pc4.mctest <- moran.mc(pca.snps$li[,4], listwgab, 999)
plot(pc4.mctest)


# Moran's Coefficient - the level of spatial autocorrelation of a quantitative variable
# obs = observed Moran’s I index calculated from the actual data
# Close to 1: Strong positive spatial autocorrelation (similar values are clustered together)
# Close to 0: Random spatial distribution (no pattern)
# Close to -1: Negative spatial autocorrelation (dissimilar values are clustered; a "checkerboard" pattern)
MC.snps <- moran.randtest(pca.snps$li, listwgab, nrepet = 9999)
MC.snps

MC.snps.tab <- data.frame("Test" = MC.snps$names,
                          "Obs" = as.numeric(MC.snps$obs),
                          "Alter" = MC.snps$alter,
                          "Pvalue" = MC.snps$adj.pvalue,
                          "nrep" = as.numeric(MC.snps$rep))
head(MC.snps.tab)

write.table(MC.snps.tab, file = "01_moran_MCI_results_PCA.tsv", sep="\t", col.names = TRUE, row.names = TRUE, append=FALSE, quote=F)


# Plotting moran coefficients
png("01_moran_PCA_MC.png",w=1200,h=600,res=150)
mc.bounds <- moran.bounds(listwgab)
snps.maps <- s1d.barchart(MC.snps$obs, labels = MC.snps$names, plot = FALSE, xlim = 1.1 * mc.bounds, paxes.draw = TRUE, pgrid.draw = FALSE)
addline(snps.maps, v = mc.bounds, plot = TRUE, pline.col = 'red', pline.lty = 3)
dev.off()


#=============================#
# Plotting AF autocorrelation #
#=============================#

# Check
# Size of the Circle: Represents the magnitude of the score for that specific 
# individual/pop on that PC axis. A large circle means that individual has a strong 
# influence (high absolute value) on that axis
# Color represents the sign (positive or negative) of the score
# Sharp contrast: A clear "line" where colors change suggests a genetic break or a barrier to gene flow.
# Gradual shift: A smooth transition from small white to large green circles 
# suggests Isolation by Distance (IBD) or a continuous environmental cline.
# Clusters: Patches of the same color in specific areas suggest refugia or distinct population groups.
# Uniformity: If all circles are small and colors are randomly mixed, 
# that specific PC axis does not have a strong spatial component.

s.value(mxy, pca.snps$li[c(1:4)], symbol = "circle", col = c("white", "palegreen4"), ppoint.cex = 0.6)

# Getting dataframe
snp.m <- as.data.frame(MC.snps$obs) %>% arrange(-MC.snps$obs) %>% head(4) %>% rownames()
snp.m.renamed <- str_replace(snp.m, "\\.statistic", "")

df.mem.snp <- cbind(mxy, pca.snps$li[snp.m.renamed])
head(df.mem.snp)
df.gat.mem.snp <- df.mem.snp %>% gather(key = "PCA_Axes", value = "value", c(3:5))
df.gat.mem.snp$fill <- ifelse(df.gat.mem.snp$value <=0, "N", "P")
df.gat.mem.snp$avalue <- abs(df.gat.mem.snp$value)
head(df.gat.mem.snp)

# morans coefficients table
df.MC <- MC.snps.tab
df.MC$PCA_Axes <- df.MC$Test
df.MC$lon <- -170
df.MC$lat <- 42
df.MC <- df.MC %>% filter(PCA_Axes %in% c("Axis1","Axis2","Axis3"))
head(df.MC)


#####################

p4 <- ggplot() + 
  geom_sf(data = world_NA, fill = "floralwhite") +
  geom_sf(data = species_distribution_outline, fill = "white") +
  geom_point(data = df.gat.mem.snp, aes(x = lon, y = lat, size = avalue, fill = fill), pch=21) +
  
  geom_text(data = df.MC, aes(x = lon, y = lat, label = paste0("Moran's I=", round(Obs,2),"; p=",round(Pvalue,4))),
            color="black", hjust=0) +
  
  scale_fill_manual(values = c("black","khaki")) +
  scale_size_continuous(name="Absolute\nvalue") +
  scale_x_continuous(limits = c(-170, -52)) +
  scale_y_continuous(limits = c(40, 70)) +
  facet_wrap(~PCA_Axes, ncol=2) +
  labs(x = "Longitude", y = "Latitude") +
  guides(fill = "none") +
  theme(panel.background = element_rect(fill=NA, colour="black"),
        panel.grid = element_blank(),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=12, face = "bold"))
p4


png("01_moran_SNPs_MEMs.png",w=2400,h=1400,res=300)
p4
dev.off()






