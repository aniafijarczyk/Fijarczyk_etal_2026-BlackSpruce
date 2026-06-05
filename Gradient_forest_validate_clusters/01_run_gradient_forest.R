rm(list=ls())
setwd("C:/Users/aniaf/Projects/BlackSpruce/952_garden_offset_run_clusters")

library(gradientForest)
library(tidyverse)
library(cowplot)
library(ggrepel)
library(readxl)
library(RColorBrewer)
library(dplyr)
library(pdist)




#=================#
#      DATA       #
#=================#

### Environment - Climate PCs for 70 provenances + 5 common gardens
dclim <- read.csv("../34_climate_transfer_distance/03_climate_distance_PCA.tsv", sep="\t", header=T)
head(dclim)
dim(dclim)

# Three climate PCs for 4 common gardens (n=4)
dclim_garden <- dclim %>% filter(Row.names %in% c("Chibougamau","Mont-Laurier","Acadia","Peace_River")) %>% dplyr::select(Row.names, PC1, PC2, PC3)
dclim_garden$SITE_ID <- c("CH","ML","AC","PR")
dclim_garden

# Three climate PCs for 70 provenances (n=70)
dclim_pop <- dclim %>% filter(!Row.names %in% c("Valcartier","Chibougamau","Mont-Laurier","Acadia","Peace_River"))
dclim_pop$POP <- as.integer(dclim_pop$Row.names)
dclim_pop <- dclim_pop %>% dplyr::select(POP, PC1, PC2, PC3)
dclim_pop %>% head()


### Train & test samples
# All filtered individuals (n=1467)
dmeta <- read.csv("../DATA_intermediate/23_filter_EPN_indiv_metrics.tsv", sep="\t", header=T)
dmeta$POP_ID <- dmeta$POP
dim(dmeta)
#dmeta$POP_SITE <- paste0(dmeta$POP, "_", dmeta$SITE_ID)

# All filtered provenances (n=66) 
dmeta <- dmeta %>% dplyr::select(POP_ID, group) %>% distinct()
dim(dmeta)
head(dmeta)

# All provenances-gardens used in train and test for each trait
dt_ <- read.csv("../44_train_and_test/00_combine_samples_means.tsv", sep="\t", header=T)
dim(dt_)

# Adding group information - will be absent from some of the test provenances
dt <- merge(dt_, dmeta, by = "POP_ID", sort=F, all.x=T)
dt$sample <- dt$POP_SITE
dim(dt)
head(dt)
unique(dt$group)


# Adding climate info to TRAIN provenances (n=1738)
#dt_TRAIN <- dt %>% dplyr::select(-Trait_name) %>% filter(SET == "TRAIN") %>% dplyr::select(POP_SITE, POP_ID, SITE_ID, group) %>% distinct()
dt_TRAIN <- dt %>% filter(SET == "TRAIN") %>% dplyr::select(Trait_name, sample, POP_ID, SITE_ID, group)
dm_TRAIN <- merge(dt_TRAIN, dclim_pop, by.x = "POP_ID", by.y = "POP", sort=F)
dim(dm_TRAIN)
head(dm_TRAIN)
dm_TRAIN  %>% group_by(Trait_name, SITE_ID) %>% summarize(n=n()) %>% filter(Trait_name == "Height")


# Adding climate info to TEST provenances (n=2048)
dt_TEST <- dt %>% filter(SET == "TEST") %>% dplyr::select(Trait_name, sample, POP_ID, SITE_ID, group)
dm_TEST <- merge(dt_TEST, dclim_pop, by.x = "POP_ID", by.y = "POP", sort=F)
dm_TEST$lp <- c(1:length(dm_TEST$POP_ID))
dim(dm_TEST)
head(dm_TEST)
dm_TEST %>% group_by(Trait_name, SITE_ID) %>% summarize(n=n()) %>% filter(Trait_name == "Height")

# Adding garden information for TEST provenances (n=2048)
dfut_TEST <- merge(dm_TEST[,c("sample","SITE_ID","Trait_name","lp","group")], dclim_garden, by = "SITE_ID", sort=F)
dfut_TEST <- dfut_TEST %>% arrange(lp)
dim(dfut_TEST)
head(dfut_TEST)
dfut_TEST %>% group_by(Trait_name, SITE_ID) %>% summarize(n=n()) %>% filter(Trait_name == "Height")


### Samples (provenances) for calculating genetic offset (with current and future/garden climate)
meta_current <- dm_TEST %>% dplyr::select(Trait_name, sample, SITE_ID, group, PC1, PC2, PC3)
meta_future <- dfut_TEST %>% dplyr::select(Trait_name, sample, SITE_ID, group, PC1, PC2, PC3)
dim(meta_current)
dim(meta_future)

# Are rownames consistent?
#dim(meta_future)
#dim(meta_current)
sum(meta_future$Trait_name == meta_current$Trait_name)
sum(meta_future$sample == meta_current$sample)




dt %>% filter(Trait_name == "Height") %>% filter(SITE_ID == "AC") %>% filter(SET=="TEST")
dt %>% filter(Trait_name == "Height") %>% filter(SITE_ID == "AC") %>% filter(SET=="TEST") %>% dim()
dt %>% filter(Trait_name == "Height") %>% filter(SITE_ID == "AC") %>% filter(SET=="TRAIN") %>% dim()
dt %>% filter(Trait_name == "Height") %>% filter(SITE_ID == "AC") %>% filter(POP_SITE == "1530_AC")


#=================#
#     FUNCIONS    #
#=================#


#---------------------#
#     Input preps     #
#---------------------#

prepare_inputs <- function(trait, filepath, garden, group_name) {
  
  filename_snps <-  filepath
  
  # Subset of TRAIN samples from the garden, trait; only in the group
  df.snps <- read.csv(filename_snps, sep="\t", header=TRUE, row.names=1)

  dt.sub <- dm_TRAIN %>% 
    filter(sample %in% rownames(df.snps)) %>% 
    filter(SITE_ID == garden) %>%
    filter(Trait_name == trait) %>% 
    filter(group %in% group_name)
  samples_train <- dt.sub %>% pull(sample)
  
  # Filtering SNP data
  dff.snps <- df.snps %>% filter(rownames(df.snps) %in% samples_train)
  print(dim(dff.snps))

  # MAF filter 0.05
  snps.freqs <- colSums(dff.snps)/(nrow(dff.snps))
  snps.maf <- c(which(snps.freqs < 0.05), which(snps.freqs > (1-0.05)))
  dff.filtered <- dff.snps %>% select(-all_of(names(snps.maf)))
  
  # Sorting
  df.sort <- data.frame("sample" = rownames(dff.filtered))
  # current train samples
  df.sort.current.train <- merge(df.sort, dt.sub, by = "sample", sort=F, all.x = TRUE)
  print(dim(df.sort.current.train))

  # Check if order matches
  print("Checking order of metadata")
  print(sum(rownames(dff.filtered) == df.sort.current.train$sample)/length(df.sort.current.train$sample))

  output <- list()
  output[["snps"]] <- dff.filtered
  output[["meta"]] <- df.sort.current.train

  write.table(output[["meta"]], file = paste0("./samples/01_run_gradient_forest_TRAIN_subset_",garden,"_",group_name,"_",trait,".tsv"),
              sep="\t", row.names = F, col.names = T, quote=F, append=F)
  
  
  return(output)
}



#---------------------#
#   Gradient forest   #
#---------------------#


run_gradient_forest <- function(datasets_meta, datasets_snps, garden, group_name, trait) {
  df.env <- datasets_meta %>% dplyr::select(PC1, PC2, PC3)
  
  ### 100
  nSites <- dim(datasets_snps)[1]
  nSpecs <- dim(datasets_snps)[2]
  lev <- floor(log2(nSites * 0.368/2))
  gf <- gradientForest(data = cbind(df.env, datasets_snps), predictor.vars = c("PC1","PC2", "PC3"),
                           response.vars = colnames(datasets_snps), ntree = 500, transform = NULL,
                           corr.threshold = 0.5, compact = T, nbin = 201, maxLevel = lev)
  save(gf, file = paste0("./data/01_run_gradient_forest_",garden,"_",group_name,"_",trait,".RData"))
  return(gf)
}




#---------------------#
#   Genetic offset    #
#---------------------#


calculate_genetic_offset <- function(gradient_forest_output, garden, group_name, trait) {

  current_variables <- meta_current[((meta_current$Trait_name == trait) & (meta_current$SITE_ID == garden)),]
  future_variables <- meta_future[((meta_future$Trait_name == trait) & (meta_future$SITE_ID == garden)),]
  rownames(current_variables) <- current_variables$sample
  rownames(future_variables) <- future_variables$sample
  sum(rownames(current_variables) == rownames(future_variables))/length(future_variables$PC1)
  
  # Prediction results - allelic turnover
  predicted.now <- predict(gradient_forest_output, current_variables[c("PC1","PC2","PC3")])
  predicted.new <- predict(gradient_forest_output, future_variables[c("PC1","PC2","PC3")])
  
  # Changing columns names
  colnames(predicted.now) <- c("now_at1","now_at2","now_at3")
  colnames(predicted.new) <- c("new_at1","new_at2","new_at3")
  
  #dout <- df.env
  dout <- cbind(current_variables, predicted.now)
  dout <- cbind(dout, predicted.new)
  
  get_distance = function(x, output){
    A = x[c("now_at1","now_at2","now_at3")]
    B = x[c("new_at1","new_at2","new_at3")]
    ed <- pdist(A,B)@dist
    return(ed)
  }
  
  dout$offset <- apply(dout, 1, get_distance)
  
  # Saving table
  write.table(dout, file = paste0("./results/01_run_gradient_forest_",garden,"_",group_name,"_",trait,".tsv"),
              sep="\t", row.names = F, col.names = T, quote=F, append=F)
  
  return(dout)
}



#--------------------------#
#      ONE EXAMPLE         #
#--------------------------#

trait = "Height"
path = "../945_GEA_lfmm/"
filecore = "05_outlier_genotypes_"
dset = "1000"
filename = paste0(path, filecore, dset, "_", trait, ".tsv")
print(filename)

garden <- "CH"
group_name <- "Central"


run.inputs <- prepare_inputs(trait, filename, garden, group_name)
head(run.inputs)
run.inputs$meta
dim(run.inputs$meta)

run.gradient <- run_gradient_forest(run.inputs[["meta"]], run.inputs[['snps']], garden, group_name, trait)
run.offset <- calculate_genetic_offset(run.gradient, garden, group_name, trait)
dim(run.offset)
run.offset


#--------------------------#
#   RUN GRADIENT FOREST    #
#--------------------------#


#run_gradient_forest
traits <- c("Height", "Biomass_Increment", "Biomass_Increment_1980", "Biomass_Increment_1985", 
           "Biomass_Increment_1990", "Biomass_Increment_1995", "Biomass_Increment_2000", 
           "Biomass_Increment_2005", "Biomass_Increment_2010", "Biomass_Increment_2015",
           "Average_Ring_Density", "DBH", "Rr", "Rs", "Rl", "Rc")



#################################################################################

### Relatively quick runs

for (trait in traits) {
  print(trait)
  
  dset = "1000"
  path = "../945_GEA_lfmm/"
  filecore = "05_outlier_genotypes_"
  filename = paste0(path, filecore, dset, "_", trait, ".tsv")
  print(filename)

  for (garden in c("AC","CH","ML","PR")) {
    for (group_name in c("East", "Central", "West")) {
      print(paste(trait, garden, group_name))
      run.inputs <- prepare_inputs(trait, filename, garden, group_name)

      run.gradient <- try(run_gradient_forest(run.inputs[["meta"]], run.inputs[['snps']], garden, group_name, trait), silent=TRUE)
      if (inherits(run.gradient, "try-error")) {
        message("Skipped element due to error")
        next
      }
      
      run.offset <- calculate_genetic_offset(run.gradient, garden, group_name, trait)
 
    }
  }
}



#################################################################################



### Relatively quick runs
# 
# for (trait in traits) {
#   print(trait)
#   
#   dset = "1000"
#   path = "../45_GEA_lfmm/"
#   filecore = "05_outlier_genotypes_"
#   filename = paste0(path, filecore, dset, "_", trait, ".tsv")
#   print(filename)
#   for (group_name in c("NoneExcluded")) {
#     print(group_name)
#     run.inputs <- prepare_inputs(trait, filename, group_name)
#     run.gradient <- run_gradient_forest(run.inputs[["meta"]], run.inputs[['snps']], trait, dset, group_name)
#     run.offset <- calculate_genetic_offset(run.gradient, trait, dset, group_name)
#     
#   }
# }



# 
# for (trait in traits) {
#   print(trait)
#   
#   for (dset in c("RDA", "RDAcorrected")) {
#     print(dset)
#     
#     path = "../46_GEA_RDA/"
#     filecore = "02_outlier_genotypes_"
#     filename = paste0(path, filecore, dset, "_", trait, ".tsv")
#     print(filename)
#     
#     run.inputs <- prepare_inputs(trait, filename)
#     run.gradient <- run_gradient_forest(run.inputs[["meta"]], run.inputs[['snps']], trait, dset)
#     run.offset <- calculate_genetic_offset(run.gradient, trait, dset)
#   }
# }


#################################################################################

### Long runs

# 
# for (trait in traits) {
#   print(trait)
#   
#   for (dset in c("10000")) {
#     print(dset)
#     
#     path = "../45_GEA_lfmm/"
#     filecore = "05_outlier_genotypes_"
#     filename = paste0(path, filecore, dset, "_", trait, ".tsv")
#     print(filename)
#     run.inputs <- prepare_inputs(trait, filename)
#     run.gradient <- run_gradient_forest(run.inputs[["meta"]], run.inputs[['snps']], trait, dset)
#     run.offset <- calculate_genetic_offset(run.gradient, trait, dset)
#   }
# }
# 
# 
# for (trait in traits) {
#   print(trait)
#   
#   for (dset in c("all")) {
#     print(dset)
#     
#     path = "../45_GEA_lfmm/"
#     filecore = "02_lea_inputs_"
#     filename.2 = paste0(path, filecore, "_", trait, ".tsv")
#     print(filename.2)
#     
#     run.inputs <- prepare_inputs(trait, filename.2)
#     run.gradient <- run_gradient_forest(run.inputs[["meta"]], run.inputs[['snps']], trait, dset)
#     run.offset <- calculate_genetic_offset(run.gradient, trait, dset)
#   }
# }
# 

















