rm(list=ls())
#setwd("C:/Users/aniaf/Projects/BlackSpruce/23_filter")
dir()


library(tidyverse)
library(dartR)
library(cowplot)
library(readxl)
library(stringr)
library(adegenet)





############
### DATA ###
############


### Removing mixed individuals, keeping only 1495

# Genotypes
gel <- gl.load("../DATA_intermediate/06_filtering_TOT_CR50_nosec.Rdata")
class(gel)
gel

gel@ind.names

# Table with samples
df <- read.csv("01_remove_mixed_admixture_metadata_clean.tsv", sep="\t", header=T)
dim(df)
head(df)
unique(df$SPECIES_ID)
samples.to.keep <- df$sample
samples.to.keep


# Table with strata
dadm <- read.csv("../23_filter/01_remove_mixed_admixture_metadata_clean.tsv", sep="\t", header=T)
head(dadm)
sadm <- dadm %>% filter(id %in% samples.to.keep) %>% dplyr::select(id,SPECIES_ID,POP,STRATA,BestK4,BestK6)
head(sadm)
# genetic cluster
sadm$cluster <- ifelse(sadm$BestK6 == "K6.RedSpruce",
                     gsub("K4.","",sadm$BestK4),
                     gsub("K6.","",sadm$BestK6))
# genetic group
pop_clusters <- sadm %>% group_by(POP,cluster) %>% reframe(n=n()) %>% arrange(POP, -n) %>% group_by(POP) %>% summarise(first(cluster))
colnames(pop_clusters) <- c("POP","group")
head(pop_clusters)
madm <- merge(sadm, pop_clusters, by = "POP", sort=F)
madm <- madm %>% dplyr::select(id, cluster, group)
dim(madm)
head(madm)



##########################
#  Filtering genotypes   #
##########################

gel.clean <- gl.keep.ind(gel, ind.list = samples.to.keep, recalc = TRUE, mono.rm = TRUE)
gel.clean
## 1,495 genotypes,  30,071 SNPs

# Keeping balck spruce only
gel.bs <- gl.keep.pop(gel.clean, c("EPN"), as.pop="SPECIES_ID", recalc = TRUE, mono.rm = TRUE)
gel.bs
# 1,467 genotypes,  29,977 SNPs

## Saving in dart format
gl.save(gel.clean, file="../DATA_intermediate/23_filter_clean.Rdata")
gl.save(gel.bs, file="../DATA_intermediate/23_filter_EPN.Rdata")

# open these files like this
#gel <- gl.load("../DATA_intermediate/23_filter_EPN.Rdata")



# Saving metadata

dmeta.clean <- gel.clean@other$loc.metrics
imeta.clean <- gel.clean@other$ind.metrics
dim(imeta.clean)
imeta.clean <- merge(imeta.clean, madm, by = "id", sort=F)
dim(imeta.clean)
head(imeta.clean)



head(dmeta.clean)
write.table(dmeta.clean, "../DATA_intermediate/23_filter_clean_locus_metrics.tsv", sep="\t", row.names = TRUE, col.names = TRUE, quote=FALSE, append=FALSE)
write.table(imeta.clean, "../DATA_intermediate/23_filter_clean_indiv_metrics.tsv", sep="\t", row.names = TRUE, col.names = TRUE, quote=FALSE, append=FALSE)


dmeta.bs <- gel.bs@other$loc.metrics
imeta.bs <- gel.bs@other$ind.metrics
dim(imeta.bs)
imeta.bs <- merge(imeta.bs, madm, by = "id", sort=F)
dim(imeta.bs)
head(dmeta.bs)
write.table(dmeta.bs, "../DATA_intermediate/23_filter_EPN_locus_metrics.tsv", sep="\t", row.names = TRUE, col.names = TRUE, quote=FALSE, append=FALSE)
write.table(imeta.bs, "../DATA_intermediate/23_filter_EPN_indiv_metrics.tsv", sep="\t", row.names = TRUE, col.names = TRUE, quote=FALSE, append=FALSE)





### Create genelight object
# https://www.rdocumentation.org/packages/adegenet/versions/2.1.10/topics/genlight-class

# Genotypes and metadata
gts <- data.frame(gel.clean)
gts[c(1:10), c(1:10)]
dmeta.clean <- gel.clean@other$loc.metrics
imeta.clean <- gel.clean@other$ind.metrics
imeta.clean <- merge(imeta.clean, madm, by = "id", sort=F)
imeta.clean %>% head()
dmeta.clean %>% head()

dmeta.clean$locus <- row.names(dmeta.clean)
dmeta.clean[c('CloneID', 'TagPos', 'Alleles')] <- str_split_fixed(dmeta.clean$locus, '-', 3)
head(dmeta.clean)

strata_list <- list('region' = imeta.clean$STRATA,
                    'group' = imeta.clean$group,
                    'cluster' = imeta.clean$cluster)
strata_list

genl <- new("genlight", gts, ploidy=2,
            ind.names=imeta.clean$id,
            loc.names=dmeta.clean$locus,
            loc.all=dmeta.clean$Alleles,
            chromosome=dmeta.clean$CloneID,
            position=dmeta.clean$TagPos,
            pop=imeta.clean$POP,
            other=strata_list)


genl@other


save(genl, file="../DATA_intermediate/23_filter_clean_genlight.Rdata")







### No red spruce

# Genotypes and metadata
gts <- data.frame(gel.bs)
gts[c(1:10), c(1:10)]
dim(gts)
dmeta <- gel.bs@other$loc.metrics
imeta <- gel.bs@other$ind.metrics
imeta <- merge(imeta, madm, by = "id", sort=F)
imeta %>% head()
dmeta %>% head()

dmeta$locus <- row.names(dmeta)
dmeta[c('CloneID', 'TagPos', 'Alleles')] <- str_split_fixed(dmeta$locus, '-', 3)
head(dmeta)

strata_list <- list('region' = imeta$STRATA,
                    'group' = imeta$group,
                    'cluster' = imeta$cluster)
strata_list


genl <- new("genlight", gts, ploidy=2,
            ind.names=imeta$id,
            loc.names=dmeta$locus,
            loc.all=dmeta$Alleles,
            chromosome=dmeta$CloneID,
            position=dmeta$TagPos,
            pop=imeta$POP,
            other=strata_list)

genl@other
genl



save(genl, file="../DATA_intermediate/23_filter_EPN_genlight.Rdata")



load("../DATA_intermediate/23_filter_EPN_genlight.Rdata")
genl


######### Saving to csv

genl@gen[[1]]
mat <- as.matrix(genl)
mat[1:10,1:10]
write.csv(mat, file = "../DATA_intermediate/23_filter_EPN.csv", sep=",", append = F, quote = F, row.names = T, col.names = T)



###########


other=list(1:3, letters,data.frame(2:4))
other


other <- list()
other[["a"]] <- c(1,2,3)
other[["b"]] <- c(3,4,5)
other[["c"]] <- c(5,6,9)
other
names(other)
data.frame(other)



