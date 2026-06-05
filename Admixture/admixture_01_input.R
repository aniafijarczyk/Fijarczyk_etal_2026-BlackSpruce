rm(list=ls())


library(ggplot2)
library(tidyr)
library(dplyr)
library(dartR)
library(tibble)
library(sf)
library(raster)
library(spData)
library(RColorBrewer)
library(readxl)
library(units)
library(vegan)
library(hierfstat)
library(SNPRelate)
library(stringr)


#######################
#------- DATA---------#
#######################

### Getting genotypes and samples (SNPs excluding secondary SNPs)

gel <- gl.load("../DATA_intermediate/06_filtering_TOT_CR50_nosec.Rdata")

# Samples list
dmeta <- gel@other$ind.metrics
samples <- gel@other$ind.metrics$id

df.loci <- gel$other$loc.metrics
df.loci$locus <- locNames(gel)
df.loci[c('CloneID', 'TagPos', 'Alleles')] <- str_split_fixed(df.loci$locus, '-', 3)
head(df.loci)



### Filter samples with > 50% missing data

gel.flt <- gl.filter.callrate(x = gel,
                                method = "ind",
                                threshold = 0.5,
                                mono.rm = TRUE,
                                recalc = TRUE,
                                recursive = TRUE)

gel.flt2 <- gl.filter.callrate(x = gel.flt,
                              method = "loc",
                              threshold = 0.85,
                              mono.rm = TRUE,
                              recalc = TRUE,
                              recursive = TRUE)


gel
# 2,322 genotypes,  31,294 SNPs
gel.flt
# 2,056 genotypes,  31,285 SNPs
gel.flt2
# 2,056 genotypes,  19,176 SNPs

# First set the order for old chroms

#chroms <- paste0("chr", df.loci$CloneID)
#chroms <- c(1 : length(df.loci$CloneID))


head(df.loci)
#chroms <- rep("1", dim(df.loci)[1])
chroms <- df.loci$CloneID
chroms

positions <- df.loci$TagPos
#positions <- c(1 : dim(df.loci)[1])
#positions <- positions*33
positions

#df.chr <- data.frame("chrom" = paste0("contig_", chroms))
#df.chr$pos <- positions
#df.chr$order <- c(1 : length(df.chr$chroms))
#head(df.chr)

# Chromosome new names
gel$chromosome <- as.factor(chroms)
gel$chromosome

gel$position <- as.factor(positions)
gel$position

dmeta <- gel.flt2@other$ind.metrics
pops <- dmeta$POP 
npops <- replace(pops, is.na(pops), 'popX')

gel$pop <- as.factor(npops)
gel$pop

gel$position

# Saving as plink

gl2plink(
  gel,
  plink_path = getwd(),
  bed_file = FALSE,
  outfile = "admixture_01_input",
  outpath = getwd(),
  chr_format = "character",
  pos_cM = "0",
  ID_dad = "0",
  ID_mom = "0",
  sex_code = "unknown",
  phen_value = "0",
  verbose = NULL
)

# If bed is not created run plink --make-bed to create one


### bim, fam, and ped

##################### 
# Renaming unmapped SNPs does not really work, because plink 
# does not create bed file, so I output plink as it is, 
# and then convert all the chrom names into sequential integers

# bim is ok
bim <- read.csv("admixture_01_input.bim", sep="\t", header=FALSE)
head(bim)
bim$V4 <- positions
write.table(bim, "admixture_01_input.bim", sep="\t", col.names=FALSE, row.names = FALSE, append=FALSE, quote = FALSE)


# fam create
dmeta <- gel.flt2@other$ind.metrics
head(dmeta)
#fam <- dmeta  %>% dplyr::select(POP, id)
fam <- read.csv("admixture_01_input.fam", sep=" ", header=FALSE)
head(fam)

#fam$V3 <- "0"
#fam$V4 <- "0"
fam$V5 <- "unknown"
fam$V6 <- "0"
head(fam)
write.table(fam, "admixture_01_input.fam", sep=" ", col.names=FALSE, row.names = FALSE, append=FALSE, quote = FALSE)


# ped
ped <- read.csv("admixture_01_input.ped", sep=" ", header=FALSE)
ped$V1 <- dmeta$POP
write.table(ped, "admixture_01_input.ped", sep=" ", col.names=FALSE, row.names = FALSE, append=FALSE, quote = FALSE)





# Formatting map file
map <- read.csv("admixture_01_input.map", sep=" ", header=FALSE)
head(map)
map["V1"] <- c(1:length(map$V1))
head(map)
dim(map)
#map["V4"][map['V4'] == 0] <- 1
#head(map)
write.table(map, "admixture_01_input_chr.map", sep=" ", col.names=FALSE, row.names = FALSE, append=FALSE, quote = FALSE)

# Formatting bim file
map_names <- map %>% dplyr::select(V1,V2)
head(map_names)
colnames(map_names) <- c("chr", "snp")

bim <- read.csv("admixture_01_input.bim", sep="\t", header=FALSE)
head(bim)
bim.map <- merge(bim, map_names, by.x = "V2", by.y = "snp", sort=FALSE)
head(bim)
head(bim.map)
nbim <- bim.map %>% dplyr::select(chr,V2,V3,V4,V5,V6)
head(nbim)
#nbim["V4"][nbim['V4'] == 0] <- 1
head(nbim)
nbim <- nbim %>% arrange(chr)
write.table(nbim, "admixture_01_input_chr.bim", sep="\t", col.names=FALSE, row.names = FALSE, append=FALSE, quote = FALSE)

# Now run plink --make-bed to create bed file