rm(list=ls())
setwd("C:/Users/aniaf/Projects/BlackSpruce/957B_garden_offset_run_var5")

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
#   Gradient forest   #
#---------------------#


run_gradient_forest <- function(datasets_meta, datasets_snps, trait, garden, suffix) {
  df.env <- datasets_meta %>% dplyr::select(all_of(clim.vars))
  
  ### 100
  nSites <- dim(datasets_snps)[1]
  nSpecs <- dim(datasets_snps)[2]
  lev <- floor(log2(nSites * 0.368/2))
  gf <- gradientForest(data = cbind(df.env, datasets_snps), predictor.vars = clim.vars,
                           response.vars = colnames(datasets_snps), ntree = 500, transform = NULL,
                           corr.threshold = 0.5, compact = T, nbin = 201, maxLevel = lev)
  save(gf, file = paste0("./data/01_run_gradient_forest_",garden,"_",suffix,"_",trait,".RData"))
  return(gf)
}




#---------------------#
#   Genetic offset    #
#---------------------#


calculate_genetic_offset <- function(gradient_forest_output, trait, garden, suffix) {
  
  current_variables <- meta_current[((meta_current$Trait_name == trait) & (meta_current$SITE_ID == garden)),]
  future_variables <- meta_future[((meta_future$Trait_name == trait) & (meta_future$SITE_ID == garden)),]
  rownames(current_variables) <- current_variables$sample
  rownames(future_variables) <- future_variables$sample
  sum(rownames(current_variables) == rownames(future_variables))/length(future_variables$PC1)
  
  # Prediction results - allelic turnover
  predicted.now <- predict(gradient_forest_output, current_variables[clim.vars])
  predicted.new <- predict(gradient_forest_output, future_variables[clim.vars])
  
  # Changing columns names
  
  now.clim.vars <- c()
  for (i in 1:length(clim.vars)) { 
    j = paste0("now_",clim.vars[i])
    now.clim.vars[i] <- j
  }
  
  new.clim.vars <- c()
  for (i in 1:length(clim.vars)) { 
    j = paste0("new_",clim.vars[i])
    new.clim.vars[i] <- j
  }

  colnames(predicted.now) <- now.clim.vars
  colnames(predicted.new) <- new.clim.vars
  
  #dout <- df.env
  dout <- cbind(current_variables, predicted.now)
  dout <- cbind(dout, predicted.new)
  
  get_distance = function(x, output){
    A = x[now.clim.vars]
    B = x[new.clim.vars]
    ed <- pdist(A,B)@dist
    return(ed)
  }
  
  dout$offset <- apply(dout, 1, get_distance)
  
  # Saving table
  write.table(dout, file = paste0("01_run_gradient_forest_",garden,"_",suffix,"_",trait,".tsv"),
              sep="\t", row.names = F, col.names = T, quote=F, append=F)
  
  return(dout)
}



#--------------------------#
#      ONE EXAMPLE         #
#--------------------------#

trait = "Biomass_Increment"
garden = "ML"
path = "../945_GEA_lfmm/"
filecore = "05_outlier_genotypes_"
dset = "1000"
suffix = "1000"
filename = paste0(path, filecore, dset, "_", trait, ".tsv")
print(filename)


run.inputs <- prepare_inputs(trait, filename, garden, suffix)
run.inputs$meta %>% group_by(SITE_ID) %>% dplyr::summarise(n = n())
gradient_forest_output <- run_gradient_forest(run.inputs[["meta"]], run.inputs[['snps']], trait, garden, suffix)
gradient_forest_output

run.offset <- calculate_genetic_offset(gradient_forest_output, trait, garden, suffix)
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

      run.gradient <- try(run.gradient <- run_gradient_forest(run.inputs[["meta"]], run.inputs[['snps']], trait, garden, suffix), silent=TRUE)
      if (inherits(run.gradient, "try-error")) {
        message("Skipped element due to error")
        next
      }
      
      run.offset <- calculate_genetic_offset(run.gradient, trait, garden, suffix)           
      
    }
  }
}



### END
#################################################################################


