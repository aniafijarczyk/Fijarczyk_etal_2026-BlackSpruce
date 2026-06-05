rm(list=ls())
setwd("C:/Users/aniaf/Projects/BlackSpruce/959B_run_rda_var5")

library(gradientForest)
library(tidyverse)
library(cowplot)
library(ggrepel)
library(readxl)
library(RColorBrewer)
library(dplyr)
library(pdist)
library(vegan)
library(LEA)




#=================#
#      DATA       #
#=================#

dclim <- read.csv("../31_climate/02_garden_climate_with_gardens.tab", sep="\t", header=T)
head(dclim)
dim(dclim)

# climate
clim.vars <- c('dPP', 'CMI', 'fallMT')


dclim_garden <- dclim %>% filter(KeyID %in% c("Chibougamau","MontLaurier","Acadia","PeaceRiver")) %>% 
  dplyr::select(all_of(c("KeyID", clim.vars)))
dclim_garden
dclim_garden$SITE_ID <- c("AC","CH","ML","PR")
dclim_garden

# Three climate PCs for 70 provenances (n=70)
dclim_pop <- dclim %>% filter(dataset == "Provenances")
dclim_pop$POP <- as.integer(dclim_pop$KeyID)
dclim_pop <- dclim_pop %>% dplyr::select(all_of(c("POP", clim.vars)))
dclim_pop %>% head()
dclim_pop$POP

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

### Removing groups ME and WI !!!!!!!!!!!
###dm_TRAIN <- dm_TRAIN %>% filter(!group %in% c("ME","WI"))

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
meta_current <- dm_TEST %>% dplyr::select(all_of(c("Trait_name", "sample", "SITE_ID", "group", clim.vars)))
meta_future <- dfut_TEST %>% dplyr::select(all_of(c("Trait_name", "sample", "SITE_ID", "group", clim.vars)))
dim(meta_current)
dim(meta_future)

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


prepare_inputs <- function(trait, filepath, garden, suffix) {
  
  filename_snps <-  filepath
  
  # Subset of TRAIN samples from the garden, trait; not in the group
  df.snps <- read.csv(filename_snps, sep="\t", header=TRUE, row.names=1)
  
  dt.sub <- dm_TRAIN %>% 
    filter(sample %in% rownames(df.snps)) %>% 
    filter(SITE_ID == garden) %>%
    filter(Trait_name == trait)
  samples_train <- dt.sub %>% pull(sample)
  
  # Filtering SNP data
  dff.snps <- df.snps %>% filter(rownames(df.snps) %in% samples_train)
  print(dim(dff.snps))
  
  # MAF filter 0.05
  snps.freqs <- colSums(dff.snps)/(nrow(dff.snps))
  snps.maf <- c(which(snps.freqs < 0.05), which(snps.freqs > (1-0.05)))
  dff.filtered <- dff.snps %>% dplyr::select(-all_of(names(snps.maf)))
  
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
  
  write.table(output[["meta"]], file = paste0("./samples/01_run_gradient_forest_TRAIN_subset_",garden,"_",suffix,"_",trait,".tsv"),
              sep="\t", row.names = F, col.names = T, quote=F, append=F)
  
  
  return(output)
}





#---------------------#
#   Genetic offset    #
#---------------------#


calculate_genetic_offset <- function(datasets_meta, datasets_snps, trait, garden, suffix) {
  
  
  # calculate the model using current climate and allele frequencies
  rda_model <- rda(datasets_snps ~ dPP + CMI + fallMT, data = datasets_meta)
  weights <- rda_model$CCA$eig/sum(rda_model$CCA$eig)
  
  # check for overfitting (if unconstrained variance == 0, the model is overfitted)
  unconstrained_var <- rda_model$CA$tot.chi
  
  # getting current and future variables
  current_variables <- meta_current[((meta_current$Trait_name == trait) & (meta_current$SITE_ID == garden)),]
  future_variables <- meta_future[((meta_future$Trait_name == trait) & (meta_future$SITE_ID == garden)),]
  rownames(current_variables) <- current_variables$sample
  rownames(future_variables) <- future_variables$sample
  sum(rownames(current_variables) == rownames(future_variables))/length(future_variables$PC1)
  
  # predicting current and future "genomic space" based on the model
  predicted.now <- predict(rda_model, newdata = current_variables[clim.vars], type = "lc")
  predicted.new <- predict(rda_model, newdata = future_variables[clim.vars], type = "lc")
  
  # Changing columns names
  df.now <- predicted.now %>% 
    as.data.frame() %>%
    mutate(
      # 'across' loops through any column starting with "RDA"
      across(
        starts_with("RDA"), 
        # .x is the column data, cur_column() gets its name (e.g., "RDA3")
        # We extract the number from the column name to use as the weight index
        ~ .x * weights[as.numeric(gsub("RDA", "", cur_column()))],
        # This renames the columns from "RDA1" to "now_at1" on the fly
        .names = "now_at{gsub('RDA', '', .col)}"
      )
    ) %>% 
    # Drop the original RDA columns dynamically
    dplyr::select(-starts_with("RDA"))
  
  df.new <- predicted.new %>% 
    as.data.frame() %>%
    mutate(
      across(
        starts_with("RDA"), 
        ~ .x * weights[as.numeric(gsub("RDA", "", cur_column()))],
        .names = "new_at{gsub('RDA', '', .col)}"
      )
    ) %>% 
    dplyr::select(-starts_with("RDA"))
  
  #dout <- df.env
  dout <- cbind(current_variables, df.now)
  dout <- cbind(dout, df.new)
  
  get_distance = function(x, output){
    A = x[startsWith(names(x), "now")]
    B = x[startsWith(names(x), "new")]
    ed <- pdist(A,B)@dist
    return(ed)
  }
  
  if (unconstrained_var > 0) {
    dout$offset_rda <- apply(dout, 1, get_distance)
  } else {
    dout$offset_rda <- NA
  }
  
  # Saving table
  write.table(dout, file = paste0("./results/01_run_gradient_forest_",garden,"_",suffix,"_",trait,".tsv"),
              sep="\t", row.names = F, col.names = T, quote=F, append=F)
  
  return(dout)
}




#--------------------------#
#      ONE EXAMPLE         #
#--------------------------#

trait = "Rr"
garden = "AC"
path = "../945_GEA_lfmm/"
filecore = "05_outlier_genotypes_"
dset = "1000"
suffix = "1000"
filename = paste0(path, filecore, dset, "_", trait, ".tsv")
print(filename)


run.inputs <- prepare_inputs(trait, filename, garden, suffix)
dim(run.inputs[["meta"]])
head(run.inputs[["meta"]])
run.offset <- calculate_genetic_offset(run.inputs[["meta"]], run.inputs[['snps']], trait, garden, suffix)
run.offset %>% head()



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
  
  for (dset in c("1000")) {
    print(dset)
    suffix <- paste0(dset, "var5")
    print(suffix)
    
    for (garden in c("AC","ML","CH","PR")) {
      print(garden)
      
      path = "../945_GEA_lfmm/"
      filecore = "05_outlier_genotypes_"
      filename = paste0(path, filecore, dset, "_", trait, ".tsv")
      print(filename)

      run.inputs <- prepare_inputs(trait, filename, garden, dset)
      
      run.gradient <- try(run.gradient <- calculate_genetic_offset(run.inputs[["meta"]], run.inputs[['snps']], trait, garden, suffix), silent=TRUE)
      if (inherits(run.gradient, "try-error")) {
        message("Skipped element due to error")
        next
      }
    }
  }
}



### END
#################################################################################


