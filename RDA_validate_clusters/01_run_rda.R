rm(list=ls())
setwd("C:/Users/aniaf/Projects/BlackSpruce/992_run_rda_clusters")
dir()

library(tidyverse)
library(cowplot)
library(readxl)
library(stringr)
library(robust)
library(ade4)

library(vegan)
library(gplots)
library(LEA)
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
head(meta_current)
head(meta_future)

# Are rownames consistent?
#dim(meta_future)
#dim(meta_current)
sum(meta_future$Trait_name == meta_current$Trait_name)
sum(meta_future$sample == meta_current$sample)




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
#   Genetic offset    #
#---------------------#


calculate_genetic_offset <- function(datasets_meta, datasets_snps, trait, garden, group_name) {
  
  
  # calculate the model using current climate and allele frequencies
  rda_model <- rda(datasets_snps ~ PC1 + PC2 + PC3, data = datasets_meta)
  weights <- rda_model$CCA$eig/sum(rda_model$CCA$eig)
  
  # getting current and future variables
  current_variables <- meta_current[((meta_current$Trait_name == trait) & (meta_current$SITE_ID == garden)),]
  future_variables <- meta_future[((meta_future$Trait_name == trait) & (meta_future$SITE_ID == garden)),]
  rownames(current_variables) <- current_variables$sample
  rownames(future_variables) <- future_variables$sample
  sum(rownames(current_variables) == rownames(future_variables))/length(future_variables$PC1)
  
  # predicting current and future "genomic space" based on the model
  predicted.now <- predict(rda_model, newdata = current_variables[c("PC1","PC2","PC3")], type = "lc")
  predicted.new <- predict(rda_model, newdata = future_variables[c("PC1","PC2","PC3")], type = "lc")
  
  # Changing columns names
  df.now <- predicted.now %>% as.data.frame() %>%
    mutate(
      now_at1 = RDA1 * weights[1],
      now_at2 = RDA2 * weights[2],
      now_at3 = RDA3 * weights[3]
    ) %>% dplyr::select(-RDA1, -RDA2, -RDA3)
  
  df.new <- predicted.new %>% as.data.frame() %>%
    mutate(
      new_at1 = RDA1 * weights[1],
      new_at2 = RDA2 * weights[2],
      new_at3 = RDA3 * weights[3]
    ) %>% dplyr::select(-RDA1, -RDA2, -RDA3)
  
  #dout <- df.env
  dout <- cbind(current_variables, df.now)
  dout <- cbind(dout, df.new)
  
  get_distance = function(x, output){
    A = x[c("now_at1","now_at2","now_at3")]
    B = x[c("new_at1","new_at2","new_at3")]
    ed <- pdist(A,B)@dist
    return(ed)
  }
  
  dout$offset_rda <- apply(dout, 1, get_distance)
  
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

run.offset <- calculate_genetic_offset(run.inputs[["meta"]], run.inputs[['snps']], trait, garden, group_name)
run.offset$offset_rda
head(run.offset)



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
      
      run.offset <- try(run.offset <- calculate_genetic_offset(run.inputs[["meta"]], run.inputs[['snps']], trait, garden, group_name), silent=TRUE)
      if (inherits(run.offset, "try-error")) {
        message("Skipped element due to error")
        next
        }
     
    }
  }
}


