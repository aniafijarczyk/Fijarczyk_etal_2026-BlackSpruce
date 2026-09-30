rm(list=ls())
setwd("C:/Users/aniaf/Projects/BlackSpruceEA/43_blups")


library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(cowplot)
library(lmerTest)



#===========#
#   DATA    #
#===========#

# All phenotyped trees n = ~2200
dpheno = read.csv("../../BlackSpruce/33_phenotypes/08_combine_phenotypes_table.tsv", sep="\t", header=T)
head(dpheno)
dpheno.sub <- dpheno %>% dplyr::select(genotype, POP_ID, SITE_ID, Block_for_analyses, Average_Ring_Density, Biomass_Increment,
                                       Biomass_Increment_1980, Biomass_Increment_1985, Biomass_Increment_1990, Biomass_Increment_1995,
                                       Biomass_Increment_2000, Biomass_Increment_2005, Biomass_Increment_2010, Biomass_Increment_2015)

dpheno.sub$POP_ID <- as.character(dpheno.sub$POP_ID)
#dpheno.sub$POP_SITE <- paste0(dpheno.sub$POP_ID,"_",dpheno.sub$SITE_ID)
head(dpheno.sub)
dim(dpheno.sub)
colnames(dpheno.sub)

# Any rows with missing data
dpheno.sub[!complete.cases(dpheno.sub), ]
dim(dpheno.sub)

# Gather 
dg <- dpheno.sub %>% gather(key = "Trait", value = "value", -genotype, -POP_ID, -SITE_ID, -Block_for_analyses)
head(dg)
dim(dg)

# Removing rows with missing data
dset <- dg[complete.cases(dg),]
dim(dset)

# Checking how much data there is
colnames(dset)
head(dset)
dblocks <- dset %>% dplyr::select(genotype, POP_ID, SITE_ID, Block_for_analyses) %>% distinct() %>% group_by(SITE_ID, Block_for_analyses) %>% summarize(n_pops = n_distinct(POP_ID), n_ind = n_distinct(genotype))
dblocks

dset$SITE_ID_Block <- paste0(dset$SITE_ID, "_", dset$Block_for_analyses)
dblocks$SITE_ID_Block <- paste0(dblocks$SITE_ID, "_", dblocks$Block_for_analyses)
dblocks

# Filtering out samples with low number of prvenances per block (and low number of ind)
blocks_to_remove <- dblocks %>% filter(n_ind<10) %>% pull(SITE_ID_Block)
blocks_to_remove

dset %>% filter(SITE_ID_Block %in% blocks_to_remove)
dset <- dset %>% filter(!SITE_ID_Block %in% blocks_to_remove)
dset$Block_for_analyses <- droplevels(factor(dset$Block_for_analyses))
dset$Block_for_analyses
head(dset)
dim(dset)

# Logs of values
head(dset)
dset$value
dset$log_value <- log(dset$value)
head(dset)


# Split into gardens
dset_ch <- dset %>% filter(SITE_ID == "CH")
dset_ac <- dset %>% filter(SITE_ID == "AC")
dset_ml <- dset %>% filter(SITE_ID == "ML")
dset_pr <- dset %>% filter(SITE_ID == "PR")
head(dset_pr)



#==================#
#     FUNCTIONS    #
#==================#


dset_ch$Height
head(dset_ch)
dset_ch$Block_for_analyses
dset_ch$POP_ID

head(dset_ch)

head(dset_ch[c('POP_ID','value')])


run_model_full_v1B <- function(dataset, trait, trait_name) {
  
  # 1. clean data first: complete rows, block and provenance as factors, no empty levels
  d <- dataset[complete.cases(dataset[, c(trait, "Block_for_analyses", "POP_ID")]) &
                 is.finite(dataset[[trait]]), ]
  d$Block_for_analyses <- droplevels(factor(d$Block_for_analyses))
  d$POP_ID        <- droplevels(factor(d$POP_ID))
  
  # 2. fit model on the cleaned data
  f   <- as.formula(paste(trait, "~ Block_for_analyses + (1|POP_ID) + (1|Block_for_analyses:POP_ID)"))
  mod <- lmer(f, data = d, control = lmerControl(optimizer = "bobyqa"))
  stopifnot(nrow(d) == nobs(mod))
  
  # 3. block effect from the model's own design matrix (always conformable)
  X         <- getME(mod, "X")
  block_eff <- as.vector(X %*% fixef(mod)) - fixef(mod)[["(Intercept)"]]
  
  # 4. plot (block x provenance) BLUP for each tree
  re_plot  <- ranef(mod)$`Block_for_analyses:POP_ID`
  plot_eff <- re_plot[paste(d$Block_for_analyses, d$POP_ID, sep = ":"), 1]
  
  # 5. adjusted value - use this for downstream
  d$value_corrected <- d[[trait]] - block_eff - plot_eff
  
  # B. provenance blup = intercept + average block effect + provenance BLUP - one value per provenance
  # baseline = intercept + average block effect (trait scale of an "average" block)
  b        <- fixef(mod)
  mean_blk <- mean(c(0, b[-1]))          # block 1 = 0, others relative to it
  baseline <- b[["(Intercept)"]] + mean_blk
  
  # provenance BLUP for each tree
  re_pop   <- ranef(mod)$`POP_ID`
  pop_blup <- re_pop[as.character(d$POP_ID), 1]
  stopifnot(!anyNA(pop_blup))
  
  # provenance value (same for all trees of a provenance)
  d$value_blup <- baseline + pop_blup
  
  # C. raw value + provenance BLUP - do not use it
  d$value_old <- d[[trait]] + pop_blup
  
  df_combined <- d[, c("genotype","POP_ID","Trait","value_corrected","value_blup","value_old")]

  return(df_combined)
}


#==================#
#   MODEL FITTING  #
#==================#

### Chibougamau
traits <- c("Average_Ring_Density", "Biomass_Increment", "Biomass_Increment_1980","Biomass_Increment_1985",
            "Biomass_Increment_1990", "Biomass_Increment_1995", "Biomass_Increment_2000",
            "Biomass_Increment_2005", "Biomass_Increment_2010", "Biomass_Increment_2015")


outputs <- list()
for (trait in c("Biomass_Increment", "Biomass_Increment_1980","Biomass_Increment_1985",
                "Biomass_Increment_1990", "Biomass_Increment_1995", "Biomass_Increment_2000",
                "Biomass_Increment_2005", "Biomass_Increment_2010", "Biomass_Increment_2015")) {
  trait.mod <- run_model_full_v1B(dset_ch[dset_ch$Trait == trait,], "log_value", trait)
  outputs[[trait]] <- trait.mod
}

for (trait in c("Average_Ring_Density")) {
  trait.mod <- run_model_full_v1B(dset_ch[dset_ch$Trait == trait,], "value", trait)
  outputs[[trait]] <- trait.mod
}

outputs[["Biomass_Increment"]] %>% head()
outputs[["Average_Ring_Density"]] %>% head()

combined_df <- do.call(rbind, outputs)
rownames(combined_df) <- NULL
head(combined_df)

dcomb_ch <- merge(dset_ch, combined_df, by = c("genotype","POP_ID","Trait"), sort=F)
head(dcomb_ch)



### Mont Laurier

outputs <- list()
for (trait in c("Biomass_Increment", "Biomass_Increment_1980","Biomass_Increment_1985",
                "Biomass_Increment_1990", "Biomass_Increment_1995", "Biomass_Increment_2000",
                "Biomass_Increment_2005", "Biomass_Increment_2010", "Biomass_Increment_2015")) {
  trait.mod <- run_model_full_v1B(dset_ml[dset_ml$Trait == trait,], "log_value", trait)
  outputs[[trait]] <- trait.mod
}

for (trait in c("Average_Ring_Density")) {
  trait.mod <- run_model_full_v1B(dset_ml[dset_ml$Trait == trait,], "value", trait)
  outputs[[trait]] <- trait.mod
}


combined_df <- do.call(rbind, outputs)
rownames(combined_df) <- NULL
head(combined_df)

dcomb_ml <- merge(dset_ml, combined_df, by = c("genotype","POP_ID","Trait"), sort=F)
head(dcomb_ml)




### Acadia

outputs <- list()
for (trait in c("Biomass_Increment", "Biomass_Increment_1980","Biomass_Increment_1985",
                "Biomass_Increment_1990", "Biomass_Increment_1995", "Biomass_Increment_2000",
                "Biomass_Increment_2005", "Biomass_Increment_2010", "Biomass_Increment_2015")) {
  trait.mod <- run_model_full_v1B(dset_ac[dset_ac$Trait == trait,], "log_value", trait)
  outputs[[trait]] <- trait.mod
}

for (trait in c("Average_Ring_Density")) {
  trait.mod <- run_model_full_v1B(dset_ac[dset_ac$Trait == trait,], "value", trait)
  outputs[[trait]] <- trait.mod
}


combined_df <- do.call(rbind, outputs)
rownames(combined_df) <- NULL
head(combined_df)

dcomb_ac <- merge(dset_ac, combined_df, by = c("genotype","POP_ID","Trait"), sort=F)
head(dcomb_ac)
dim(dcomb_ac)



### Peace River

outputs <- list()
for (trait in c("Biomass_Increment", "Biomass_Increment_1980","Biomass_Increment_1985",
                "Biomass_Increment_1990", "Biomass_Increment_1995", "Biomass_Increment_2000",
                "Biomass_Increment_2005", "Biomass_Increment_2010","Biomass_Increment_2015")) {
  print(trait)
  trait.mod <- run_model_full_v1B(dset_pr[dset_pr$Trait == trait,], "log_value", trait)
  outputs[[trait]] <- trait.mod
}


for (trait in c("Average_Ring_Density")) {
  trait.mod <- run_model_full_v1B(dset_pr[dset_pr$Trait == trait,], "value", trait)
  outputs[[trait]] <- trait.mod
}

combined_df <- do.call(rbind, outputs)
rownames(combined_df) <- NULL
head(combined_df)

dcomb_pr <- merge(dset_pr, combined_df, by = c("genotype","POP_ID","Trait"), sort=F)
head(dcomb_pr)
unique(dcomb_pr$Trait)

###############


### Combine datasets

dcombo <- rbind(dcomb_ch, dcomb_ml, dcomb_ac, dcomb_pr)
head(dcombo)

write.table(dcombo, "10_correct_AB_phenotypes_corrected.tsv", sep="\t", col.names=T, row.names = F, quote=F, append=F)




