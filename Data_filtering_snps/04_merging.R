rm(list=ls())
#setwd("//wsl.localhost/Ubuntu/home/BlackSpruce/01_snp_treatment")
dir()


library(tidyverse)
library(dartR)
library(cowplot)
library(readxl)
library(stringr)

############
### DATA ###
############

gl1 <- gl.read.dart("../DATA/DartSeq_290323/Report_DSpr20-4991_2_moreOrders_SNP_mapping_1.csv")
gl1@n.loc
gl1@ind.names
head(gl1@other$ind.metrics)
head(gl1@other$loc.metrics)
as.data.frame(gl1)[1:40, 1:10]
# n ind = 2347
# n.loc = 36112

gl2 <- gl.read.dart("../DATA/DartSeq_290323/Report_DSpr20-4991_2_moreOrders_SNP_mapping_2.csv")
gl2@n.loc
gl2@ind.names
head(gl2@other$ind.metrics)
head(gl2@other$loc.metrics)
gl2
# n ind = 2347
# n.loc = 40090

gl3 <- gl.read.dart("../DATA/DartSeq_290323/Report_DSpr20-4991_2_moreOrders_SNP_mapping_3.csv")
gl3@n.loc
gl3@ind.names
head(gl3@other$ind.metrics)
head(gl3@other$loc.metrics)
gl3

# n ind = 2347
# n.loc = 42450


########################
### OVERLAPPING SNPS ###
########################


### Checking overlap of SNPs (contig locations) between datasets
# There is overlap between and within datasets (mappings)
# I'm skipping this filtering step for now

locNames(gl1)

# Selecting contig coordinates of SNPs
loc1 <- gl1@other$loc.metrics %>% dplyr::select(CloneID, SNP, CallRate, Chrom_Spruce_v1, ChromPosSnp_Spruce_v1)
loc1$locName <- paste0(loc1$CloneID, '-', str_replace(loc1$SNP, ":","-"))
loc1$locPos <- paste0(loc1$Chrom_Spruce_v1, '_', loc1$ChromPosSnp_Spruce_v1)
loc1$mapping <- 'M1'
loc1 <- loc1 %>% filter(locPos != '_0')
loc1 <- loc1 %>% dplyr::select(locName, locPos, CallRate, mapping)
loc1 %>% head()

loc2 <- gl2@other$loc.metrics %>% dplyr::select(CloneID, SNP, CallRate, Chrom_Spruce_v1, ChromPosSnp_Spruce_v1)
loc2$locName <- paste0(loc2$CloneID, '-', str_replace(loc2$SNP, ":","-"))
loc2$locPos <- paste0(loc2$Chrom_Spruce_v1, '_', loc2$ChromPosSnp_Spruce_v1)
loc2$mapping <- 'M2'
loc2 <- loc2 %>% filter(locPos != '_0')
loc2 <- loc2 %>% dplyr::select(locName, locPos, CallRate, mapping)

loc3 <- gl3@other$loc.metrics %>% dplyr::select(CloneID, SNP, CallRate, Chrom_Spruce_v1, ChromPosSnp_Spruce_v1)
loc3$locName <- paste0(loc3$CloneID, '-', str_replace(loc3$SNP, ":","-"))
loc3$locPos <- paste0(loc3$Chrom_Spruce_v1, '_', loc3$ChromPosSnp_Spruce_v1)
loc3$mapping <- 'M3'
loc3 <- loc3 %>% filter(locPos != '_0')
loc3 <- loc3 %>% dplyr::select(locName, locPos, CallRate, mapping)

# Merging 3 mappings
locs <- rbind(loc1, loc2, loc3)
head(locs)

# SNPs with identical coordinates across mappings
locs_grouped <- locs %>% group_by(locPos) %>% summarize(n = n(), maps = toString(mapping))
locs_multiple <- locs_grouped %>% arrange(-n) %>% filter(n > 1)
locs_multiple




##################
### READ DEPTH ###
##################

# Removing SNPs with read depth > 30

# Function to drop SNPs with read depth higher than maxRD

filterRdepth <- function(gel, maxRD) {
  gel_meta <- gel@other$loc.metrics
  gel_meta$locus <- locNames(gel)
  gel_rdepth_30 <- gel_meta %>% filter(rdepth > maxRD) %>% pull(locus)
  
  ngel <- gl.drop.loc(x = gel, gel_rdepth_30)
  ngel <- gl.filter.monomorphs(ngel)
  ngel <- gl.recalc.metrics(ngel)
  
  return(ngel)
}

gl1_rd <- filterRdepth(gl1, 30)
gl2_rd <- filterRdepth(gl2, 30)
gl3_rd <- filterRdepth(gl3, 30)

gl1_rd@n.loc 
#36089
gl2_rd@n.loc 
#40065
gl3_rd@n.loc
#42394


###############
### MERGING ###
###############

# Merging datasets
glM <- cbind(gl1_rd, gl2_rd, gl3_rd)
glM
# 2,347 genotypes,  118,548 SNPs , size: 473.5 Mb

# Recalculating metrics
glM <- gl.filter.monomorphs(glM)
glM <- gl.recalc.metrics(glM)
glM
# 2,347 genotypes,  118,548 SNPs , size: 557.7 Mb

# Running compliance function
glM <- gl.compliance.check(glM)


glM@other$ind.metrics %>% head()
glM@other$ind.metrics$SPECIES_ID %>% unique()
glM@n.loc
glM@ind.names %>% length()








##################
### REPLICATES ###
##################

# A function to obtain per sample call rate
ind_call_rate = function(row, output){
  x = length(row[!is.na(row)])/length(row)
  return(x)
}


# Replicates - In the spreadsheet "OUR_REPLICATES=TRUE" (in general) are doubled in the dartSeq (second name has added ".1")
# Selecting one with better quality
# Identifying replicates from the sample name ending
# Renaming samples - ending ".1" removed

removeReplicates <- function(gel) {
  
  replicates.2 <- gel@ind.names[grep("\\.1$",gel@ind.names)]
  replicates.1 <- gsub("\\.1","", replicates.2)
  df.reps <- data.frame("id" = c(replicates.1, replicates.2),"Sample_ID" = c(replicates.1, replicates.1))
  
  # Extracting data fro replicates only
  gel_reps <- gl.keep.ind(gel, df.reps$id)
  gel_reps <- gl.filter.monomorphs(gel_reps)
  gel_reps <- gl.recalc.metrics(gel_reps)
  
  # Calculating call rate per sample
  gel_reps_df <- as.data.frame(gel_reps)
  gel_reps_df$call_rate <- apply(gel_reps_df, 1, ind_call_rate)
  df_call_rate <- gel_reps_df["call_rate"]
  df_call_rate$id <- rownames(df_call_rate)
  
  # Finding samples with lower call rate
  repl_merged <- merge(df.reps, df_call_rate, by = "id", sort=FALSE) 
  repl_merged %>% dplyr::select(id,Sample_ID, call_rate) %>% arrange(Sample_ID)
  replicates_bad <- repl_merged %>% group_by(Sample_ID) %>% dplyr::arrange(call_rate, .by_group = TRUE) %>% 
    summarize(id = first(id), call_rate = first(call_rate)) %>% pull(id)
  
  # Removing replicated samples with lower call rate.
  gel.noreps <- gl.drop.ind(gel, replicates_bad)
  gel.noreps <- gl.filter.monomorphs(gel.noreps)
  gel.noreps <- gl.recalc.metrics(gel.noreps)
  
  # Renaming samples
  new_names <- gsub("\\.1","", gel.noreps@ind.names)
  df.temp <- data.frame("id" = gel.noreps@ind.names, "new_names" = new_names)
  write.table(df.temp, file = "dataPrep_01_overview_new_names_temp.csv", sep=",", row.names = FALSE, col.names = FALSE, append=FALSE)
  gel.noreps.renamed <- gl.recode.ind(gel.noreps, ind.recode="dataPrep_01_overview_new_names_temp.csv", verbose=3)
  gel.noreps.renamed@other$ind.metrics$id <- gel.noreps.renamed@ind.names
  return(gel.noreps.renamed)
  
} 




########################
### Removing replicates

glM.norep <- removeReplicates(glM)
glM.norep
# 2,323 genotypes,  118,517 SNPs , size: 487 Mb

glM.norep@other$ind.metrics %>% head()
glM.norep@other$ind.metrics %>% dim()
glM.norep





#######################
### ADDING METADATA ###
#######################

###  PREPARING METADATA  ###


prepareMetadata <- function(gel, metadata_file) {
  
  # Taking sample list from metadata
  dmeta <- read_excel(metadata_file, sheet = "Strata_file")
  dmeta.EPR <- read_excel(metadata_file, sheet = "EPR (Picea rubens)")
  # Renaming one individual in Red spruce
  # "EPR-NC-RS01-2" into "EPR-nc-RS01-2"
  dmeta.EPR <- dmeta.EPR %>% mutate(Sample = ifelse(Sample == "EPR-NC-RS01-2", "EPR-nc-RS01-2", Sample))
  
  # Removing duplicated samples in metadata
  dmeta <- dmeta %>% dplyr::select(-targetid) %>% distinct(genotype, Sample_ID, .keep_all = TRUE)
  
  # Check intersection with dartSeq samples
  dartseq <- data.frame("id" = gel@ind.names)
  dartseq$source <- "dartSeq"
  dmeta.1 <- merge(dmeta, dartseq, by.x = "genotype", by.y = "id", all.x = TRUE, all.y = TRUE, sort=FALSE)
  
  # Missing samples (not present in dmeta)
  missing.inds <- dmeta.1[is.na(dmeta.1$Sample_ID), 1]
  dartseq.1 <- data.frame("id" = missing.inds)
  dartseq.1$source <- "dartSeq"
  
  # Check intersection of the remaining dartseq samples with Red Spruce
  dmeta.EPR.1 <- merge(dmeta.EPR, dartseq.1, by.x = "Sample", by.y = "id", all.x = TRUE, all.y = TRUE, sort=FALSE)
  # Missing samples (not present in dmeta.EPR)
  missing.inds.2 <- dmeta.EPR.1[is.na(dmeta.EPR.1$POP), 1]
  dartseq.2 <- data.frame("genotype" = missing.inds.2)
  
  # Formatting dmeta.EPR
  colnames(dmeta.EPR) <- c("genotype","POP_ID","SampleNUM", "POP_NAME", "lat", "lon", "elev")
  dmeta.EPR$Sample_ID <- dmeta.EPR$genotype
  dmeta.EPR$DART_DNA_PLATE <- NA
  dmeta.EPR$DART_col <- NA
  dmeta.EPR$DART_row <- NA
  dmeta.EPR$SPECIES_ID <- "EPR"
  dmeta.EPR$SITE_ID <- NA
  dmeta.EPR$OUR_REPLICATES <- FALSE
  dmeta.EPR$DART_REPLICATES <- FALSE
  dmeta.EPR$BATCH_ID <- NA
  dmeta.EPR$DNA_EXTR_PLATE <- NA
  dmeta.EPR$TISSUE_FOR_DNA <- NA
  dmeta.EPR$altBATCH_ID <- NA
  dmeta.EPR$POP_GR <- "R"
  colnames(dmeta.EPR)
  colnames(dmeta)
  dmeta.EPR.filled <- dmeta.EPR %>% dplyr::select(colnames(dmeta))
  
  # Formatting other unknown samples
  dartseq.2$Sample_ID <- dartseq.2$genotype
  dartseq.2$DART_DNA_PLATE <- NA
  dartseq.2$DART_col <- NA
  dartseq.2$DART_row <- NA
  dartseq.2$SPECIES_ID <- "UNK"
  dartseq.2$SITE_ID <- NA
  dartseq.2$OUR_REPLICATES <- FALSE
  dartseq.2$DART_REPLICATES <- FALSE
  dartseq.2$BATCH_ID <- NA
  dartseq.2$DNA_EXTR_PLATE <- NA
  dartseq.2$TISSUE_FOR_DNA <- NA
  dartseq.2$altBATCH_ID <- NA
  dartseq.2$POP_ID <- NA
  dartseq.2$POP_GR <- "X"
  dartseq.2$lat <- NA
  dartseq.2$lon <- NA
  
  # Combining 3 meta files
  dmeta.all <- rbind(dmeta, dmeta.EPR.filled, dartseq.2)
  dmeta.all <- dmeta.all %>% dplyr::rename("id" = "genotype")
  dmeta.all <- dmeta.all %>% dplyr::rename("POP" = "POP_ID")
  dmeta.all <- dmeta.all %>% dplyr::rename("STRATA" = "POP_GR")
  
  # Saving metadata
  #write.table(df.temp, file = "dataPrep_01_overview_METADATA.csv", sep=",", row.names = TRUE, col.names = FALSE, append=FALSE)
  
  # Order metadata like in gl
  dmeta.all.ordered <- dmeta.all[match(indNames(gel), dmeta.all$id),]
  
  return(dmeta.all.ordered)
  
}



########################
### Preparing metadata

glM.dmeta <- prepareMetadata(glM.norep, "../METADATA/SampleSheet_DARTseq_BILAN_08-2022.xlsx")
glM.dmeta %>% head()
glM.dmeta %>% dim()

# Adding info to gl
glM.norep$other$ind.metrics <- merge(glM.norep$other$ind.metrics, glM.dmeta, by = "id", sort=FALSE)
glM.norep$other$ind.metrics %>% head()
glM.norep$other$ind.metrics %>% dim()

# Saving clean dataset with meta
gl.save(glM.norep, file="../DATA_intermediate/Dartseq_123_merged.Rdata")
#glM.norep <- gl.load("../DATA_intermediate/Dartseq_123_merged.Rdata")







############################
### REMAINING DUPLICATES ###
############################

### (repeated Sample_ID names, in the spreadsheet, the remaining "OUR_REPLICATES=TRUE")


df_inds <- glM.norep@other$ind.metrics %>% dplyr::select(id, Sample_ID, OUR_REPLICATES)
head(df_inds)
df_inds$REP <- as.integer(as.logical(df_inds$OUR_REPLICATES))
head(df_inds)
dg_inds <- df_inds %>% group_by(Sample_ID) %>% summarize(n = n(), reps = sum(REP)) %>% arrange(-n)
dg_inds

# These were already filtered out in the previous step - they showed up with .1
one_reps <- dg_inds %>% filter(reps == 1) %>% pull(Sample_ID) %>% as.vector()
one_reps
df_inds %>% filter(Sample_ID %in% one_reps)

# These are the remaining replicates
double_reps <- dg_inds %>% filter(reps == 2) %>% pull(Sample_ID) %>% as.vector()
double_reps

df_replicates <- df_inds %>% filter(Sample_ID %in% double_reps) 
dim(df_replicates)
df_replicates

replicates <- as.vector(df_replicates$id)
replicates
length(replicates)

# Getting CallRate for all duplicates
gel_repl <- gl.keep.ind(glM.norep, replicates)
gel_repl <- gl.filter.monomorphs(gel_repl)
gel_repl <- gl.recalc.metrics(gel_repl)
gel_repl_df <- as.data.frame(gel_repl)
gel_repl_df[1:2, 1:10]

gel_repl_df$call_rate <- apply(gel_repl_df, 1, ind_call_rate)
df_call_rate <- gel_repl_df["call_rate"]
df_call_rate$id <- rownames(df_call_rate)
df_call_rate

# Getting id names for the sample with worse callrate
repl_merged <- merge(df_replicates, df_call_rate, by = "id", sort=FALSE) 
repl_merged

replicates_bad <- repl_merged %>% group_by(Sample_ID) %>% dplyr::arrange(call_rate, .by_group = TRUE) %>% 
  summarize(id = first(id), call_rate = first(call_rate)) %>% pull(id)
replicates_bad

print("Replicated samples to remove")
length(replicates_bad)

# Filtering out duplicated samples with worse call rate
glM.norep.clean <- gl.drop.ind(glM.norep, replicates_bad)
glM.norep.clean <- gl.filter.monomorphs(glM.norep.clean)
glM.norep.clean <- gl.recalc.metrics(glM.norep.clean)
glM.norep.clean
#  2,322 genotypes,  118,516 SNPs , size: 487.1 Mb


glM.norep.clean
glM.norep.clean@n.loc
glM.norep.clean@ind.names %>% length()


# Saving clean dataset with meta
gl.save(glM.norep.clean, file="../DATA_intermediate/04_merging.Rdata")
#glM.norep.clean <- gl.load("../DATA_intermediate/04_merging.Rdata")


