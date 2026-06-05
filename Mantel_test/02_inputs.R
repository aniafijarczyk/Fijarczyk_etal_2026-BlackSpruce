rm(list=ls())
setwd("C:/Users/aniaf/Projects/BlackSpruce/925_Mantel")


library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)
library(sf)

library(RColorBrewer)
library(readxl)
library(units)

library(ade4) # mantel.rtest
library(vegan) # mantel



#######################
#------- DATA---------#
#######################

### GEO

dmeta <- read.csv("../DATA_intermediate/23_filter_clean_indiv_metrics.tsv", sep="\t")
dpops <- dmeta %>% filter(SPECIES_ID == 'EPN') %>% group_by(POP) %>% dplyr::summarise(n=n())
hist(dpops$n)
dpops %>% arrange(n)
dpops_sub <- dpops %>% filter(n>=5)
fmeta <- dmeta %>% filter(POP %in% dpops_sub$POP)
dim(fmeta)
head(fmeta)
# Sort by longitude
fmeta <- fmeta %>% arrange(lon)
head(fmeta)

fmeta_pops <- fmeta %>% dplyr::select(POP, lon, lat) %>% distinct() %>% as_tibble() %>% arrange(lon)
dim(fmeta_pops)
head(fmeta_pops)





### Fst

load("01_pairwise_fst_run.RData")

pop_order <- read.csv("01_pairwise_fst_pop_order.tsv", sep="\t", header=T)
head(pop_order)


# Read matrix and change pop names
mat <- pairwise.fst.wc$pairwise.fst.full.matrix
dim(mat)
length(pop_order$pop)
colnames(mat) <- pop_order$pop
rownames(mat) <- pop_order$pop
head(mat)

# Change to dataframe
dmat <- as.data.frame(mat)
head(dmat)
class(dmat)

# Change order of columns & rows
dmat$POP <- row.names(dmat)
head(dmat)
dm_mat <- merge(fmeta_pops, dmat, by.x = "POP", by.y = "POP", sort=F)
rownames(dm_mat) <- dm_mat$POP
head(dm_mat)
ord_mat <- dm_mat %>% dplyr::select(all_of(dm_mat$POP))
head(ord_mat)

# Change values to numeric and back to matrix
nmat <- mutate_all(ord_mat, function(x) as.numeric(as.character(x)))
head(nmat)
newmat <- as.matrix(nmat)
head(newmat)





#########################
#--Geographic distance--#
#########################

### Function for getting distance matrix in km from meta data frame

get.geo.distance <- function(x) {
  # Calculating distance in km
  pop.sf <- st_as_sf(x, coords = c("lon", "lat"), crs = 4326)
  pop.gdists <- st_distance(pop.sf, pop.sf) %>% set_units(m) %>% set_units(km)
  # Adding names to the matrix
  colnames(pop.gdists) <- pop.sf$POP
  rownames(pop.gdists) <- pop.sf$POP
  geo.dist <- as.dist(as.matrix(pop.gdists))
  return(geo.dist)
}


### All pops

geo.dist <- get.geo.distance(fmeta_pops)
geo.dist
row.names(as.matrix(geo.dist))
geo.dist
write.csv(as.matrix(geo.dist), "02_inputs_geo_matrix.csv")


#########################
#--        FST        --#
#########################

#Fst
head(newmat)
class(newmat)

# Calculate Fst/ (1 - Fst)

fst.mat <- newmat / (1 - newmat)
newmat[1:5, 1:5]
fst.mat[1:5, 1:5]
fst.dist <- as.dist(fst.mat)
class(fst.dist)
fst.dist

row.names(as.matrix(fst.dist))
write.csv(as.matrix(fst.dist), "02_inputs_transfst_matrix.csv")



#########################
#--     TRIPLETS      --#
#########################

### series of three adjecent strata

# Getting strata triplets: calculating Mantel's test on each set of 3 adjecent strata
strata_pops <- fmeta %>% dplyr::select(POP, lon, lat, STRATA) %>% distinct() %>% as_tibble() %>% arrange(lon)
head(strata_pops)

strata <- unique(fmeta$STRATA)
ord.strata <- strata[order(strata)]
ord.strata

S <- list()
for (i in 1:(length(ord.strata)-2)) {
  #print(i)
  group <- c(ord.strata[i], ord.strata[i+1], ord.strata[i+2])
  strata_ <- strata_pops[strata_pops$STRATA %in% group,]
  strata_ <- strata_ %>% arrange(lon)
  strata_$group <- i
  strata.dist <- get.geo.distance(strata_)
  write.csv(as.matrix(strata.dist), paste0("02_inputs_geo_matrix_strata_",i,".csv"))
  
  S[[i]] <- strata.dist
}

S




FST <- list()
for (i in 1:(length(ord.strata)-2)) {
  #print(i)
  group <- c(ord.strata[i], ord.strata[i+1], ord.strata[i+2])
  strata_ <- strata_pops[strata_pops$STRATA %in% group,]
  strata_ <- strata_ %>% arrange(lon)
  fst_ <- fst.mat[strata_$POP, strata_$POP]
  write.csv(fst_, paste0("02_inputs_transfst_matrix_strata_",i,".csv"))
  FST[[i]] <- as.dist(fst_)
}

FST




