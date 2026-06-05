rm(list=ls())
setwd("C:/Users/aniaf/Projects/BlackSpruce/974_provenances_offset_run")

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

### Environment
dmeta <- read.csv("../57_climate_clusters/04_all_climate_PCA_Kmeans.tsv", sep="\t", header=TRUE)
head(dmeta)
# 14 climate variables + dPP
#clim.vars <- c('CMI', 'mayMinT', 'TP', 'MAT', 'MMinT',
#               'MMaxT', 'sumTP', 'WSTP', 'winMT', 'sprMT', 'sumMT', 'fallMT',
#               'MCMT', 'MWMT', 'dPP')
clim.vars <- c("PC1","PC2","PC3")
#clim.vars <- c('dPP', 'CMI', 'TP', 'WSTP', 'fallMT')

dmeta$sample <- dmeta$KeyID
dmeta$sample
head(dmeta)
dim(dmeta)

# adding group
dgroup <- read.csv("../DATA_intermediate/23_filter_EPN_indiv_metrics.tsv", sep="\t", header=T)
dgroup$KeyID <- dgroup$POP
dgroup <- dgroup %>% dplyr::select(KeyID, group) %>% distinct()
head(dgroup)

dmeta <- merge(dmeta, dgroup, by = "KeyID", sort=F, all.x = T)
dim(dmeta)
tail(dmeta)
#dmeta %>% filter(!group %in% c("West","Central", "East", "ME","WI")) %>% pull(sample) %>% unique()
dmeta$sample

#==================#
#     FUNCTIONS    #
#==================#



#---------------------#
#     Input preps     #
#---------------------#

prepare_inputs <- function(filepath) {
  
  filename_snps <-  filepath

  ### SNP data
  df.snps <- read.csv(filename_snps, sep=",", header=TRUE, row.names=1)
  print(dim(df.snps))

  dt.sub <- dmeta %>% 
    filter(sample %in% rownames(df.snps)) %>% 
    filter(Year == 1960) %>%
    filter(Scenario == 1990)
  
  samples_train <- dt.sub %>% pull(sample)
  
  # Filtering SNP data
  df.snps$sample <- rownames(df.snps)
  dff.snps <- df.snps %>% filter(sample %in% samples_train)
  row.names(dff.snps) <- dff.snps$sample
  dff.snps <- dff.snps %>% dplyr::select(-sample)
  print(dim(dff.snps))
  
  # MAF filter 0.05
  snps.freqs <- colSums(dff.snps)/(nrow(dff.snps))
  snps.maf <- c(which(snps.freqs < 0.05), which(snps.freqs > (1-0.05)))
  dff.filtered <- dff.snps %>% dplyr::select(-all_of(names(snps.maf)))
  
  # Subset random 1000
  all_snps <- colnames(dff.filtered)
  set.seed(778343)
  random_snp <- sample(all_snps, 1000)
  dff.random <- dff.filtered[,random_snp]
  
  # sorting
  df.sort <- data.frame("sample" = rownames(dff.random))
  df.sort.meta <- merge(df.sort, dt.sub, by = "sample", sort=F)
  
  # Check if order matches
  print("Checking order of metadata")
  print(sum(rownames(dff.random) == df.sort.meta$sample)/length(df.sort.meta$sample))
  
  output <- list()
  output[["snps"]] <- dff.random
  output[["meta"]] <- df.sort.meta
  return(output)
}





#---------------------#
#   Gradient forest   #
#---------------------#


run_gradient_forest <- function(datasets_meta, datasets_snps) {
  df.env <- datasets_meta %>% dplyr::select(all_of(clim.vars))
  
  ### 
  nSites <- dim(datasets_snps)[1]
  nSpecs <- dim(datasets_snps)[2]
  lev <- floor(log2(nSites * 0.368/2))
  gf <- gradientForest(data = cbind(df.env, datasets_snps), predictor.vars = clim.vars,
                       response.vars = colnames(datasets_snps), ntree = 500, transform = NULL,
                       corr.threshold = 0.5, compact = T, nbin = 201, maxLevel = lev)
  save(gf, file = paste0("./data/02_run_gradient_forest.RData"))
  return(gf)
}



#---------------------#
#   Genetic offset    #
#---------------------#


calculate_genetic_offset <- function(gradient_forest_output, year, scenario) {
  
  current_variables <- dmeta %>% 
    filter(Year == 1960) %>%
    filter(Scenario == 1990)

  future_variables <- dmeta %>% 
    filter(Year == year) %>%
    filter(Scenario == scenario)  
  
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
  write.table(dout, file = paste0("02_run_gradient_forest_",year,"_",scenario,".tsv"),
              sep="\t", row.names = F, col.names = T, quote=F, append=F)
  
  return(dout)
  }







#=================#
#   ONE EXAMPLE   #
#=================#

# 
# path = "../58_provenance_genetic_offset_run/"
# 
# inputs <- prepare_inputs(paste0(path, "01_frequencies_POP.tsv"))
# inputs
# dim(inputs$snps)
# dim(inputs$meta)
# head(inputs$meta)
# 
# gradient <- run_gradient_forest(inputs[["meta"]], inputs[['snps']])
# gradient
# 
# offset <- calculate_genetic_offset(gradient, "2031", "245")
# head(offset)
# 

#=================#
#       LOOP      #
#=================#


path = "../926_MoransI/"
inputs <- prepare_inputs(paste0(path, "01_moran_Pop_AFs.csv"))
dim(inputs[["snps"]])
inputs[["meta"]] %>% head()
gradient <- run_gradient_forest(inputs[["meta"]], inputs[['snps']])
gradient
  
for (year in c(2031, 2050, 2070, 2090)) {
  for (scenario in c("245")) {
    print(paste(year, scenario))
    offset <- calculate_genetic_offset(gradient, year, scenario)
  }
}











#=================#
#       END       #
#=================#









