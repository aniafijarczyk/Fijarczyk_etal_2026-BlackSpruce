rm(list=ls())
setwd("C:/Users/aniaf/Projects/BlackSpruce/943X_climate_transfer_distance")


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
dpheno = read.csv("../../BlackSpruce/33_phenotypes/06_overview_phenotypes_Patrick_data.tsv", sep="\t", header=T)
dpheno$POP_ID <- as.character(dpheno$POP_ID)
dpheno$POP_SITE <- paste0(dpheno$POP_ID,"_",dpheno$SITE_ID)
head(dpheno)
dim(dpheno)
colnames(dpheno)

# Any rows with missing data
dpheno[!complete.cases(dpheno), ]
dim(dpheno)

# Removing rows with missing data
dset <- dpheno[complete.cases(dpheno),]
dim(dset)

# Checking how much data there is
colnames(dset)
dblocks <- dset %>% group_by(SITE_ID, Block_Patrick) %>% summarize(n_pops = n_distinct(POP_ID), n_ind = n())
dblocks

dset$SITE_ID_Block <- paste0(dset$SITE_ID, "_", dblocks$Block_Patrick)
dblocks$SITE_ID_Block <- paste0(dblocks$SITE_ID, "_", dblocks$Block_Patrick)
dblocks

# Filtering out samples with low number of prvenances per block (and low number of ind)
blocks_to_remove <- dblocks %>% filter(n_ind<10) %>% pull(SITE_ID_Block)

dset <- dset %>% filter(!SITE_ID_Block %in% blocks_to_remove)
dset$Block_Patrick <- droplevels(factor(dset$Block_Patrick))
head(dset)
dim(dset)

# Split into gardens
dset_ch <- dset %>% filter(SITE_ID == "CH")
dset_ac <- dset %>% filter(SITE_ID == "AC")
dset_ml <- dset %>% filter(SITE_ID == "ML")
dset_pr <- dset %>% filter(SITE_ID == "PR")
head(dset_pr)




#==================#
#     EXAMPLE      #
#==================# 


dset_ch$Height
head(dset_ch)
dset_ch$Block_Patrick
dset_ch$POP_ID
trait <- "Height"
dataset <- dset_ch

# 1. clean data first: complete rows, block and provenance as factors, no empty levels
d <- dataset[complete.cases(dataset[, c(trait, "Block_Patrick", "POP_ID")]), ]
d$Block_Patrick <- droplevels(factor(d$Block_Patrick))
d$POP_ID        <- droplevels(factor(d$POP_ID))

# 2. fit model on the cleaned data
f   <- as.formula(paste(trait, "~ Block_Patrick + (1|POP_ID) + (1|Block_Patrick:POP_ID)"))
mod <- lmer(f, data = d, control = lmerControl(optimizer = "bobyqa"))
stopifnot(nrow(d) == nobs(mod))
summary(mod)
VarCorr(mod)
fixef(mod)
fixef(mod)[["(Intercept)"]]
fixef(mod)[-1]

# 3. block effect from the model's own design matrix (always conformable)
X         <- getME(mod, "X")
X
# X %*% fixef(mod) goes through each tree's row, multiplies each 0/1 by the matching number
# intercept + my block's coefficient - intercept
# same as: predict(mod, re.form = NA) - fixef(mod)[["(Intercept)"]]
block_eff <- as.vector(X %*% fixef(mod)) - fixef(mod)[["(Intercept)"]]
block_eff

# 4. plot (block x provenance) BLUP for each tree
ranef(mod)$`POP_ID`
re_plot  <- ranef(mod)$`Block_Patrick:POP_ID`
re_plot
paste(d$Block_Patrick, d$POP_ID, sep = ":")
re_plot[paste(d$Block_Patrick, d$POP_ID, sep = ":"),1]
plot_eff <- re_plot[paste(d$Block_Patrick, d$POP_ID, sep = ":"), 1]
length(plot_eff)
plot_eff[is.na(plot_eff)] <- 0     # safety: should not happen
plot_eff

# 5. adjusted value
d$value_corrected <- d[[trait]] - block_eff - plot_eff

out <- d[, c("genotype", "POP_ID", "value_corrected")]
#colnames(out)[3] <- paste0(trait, "_corrected")
head(out)



#==================#
#     FUNCTIONS    #
#==================#



run_model_full_v1 <- function(dataset, trait) {
  
  # 1. clean data first: complete rows, block and provenance as factors, no empty levels
  d <- dataset[complete.cases(dataset[, c(trait, "Block_Patrick", "POP_ID")]) &
                 is.finite(dataset[[trait]]), ]
  d$Block_Patrick <- droplevels(factor(d$Block_Patrick))
  d$POP_ID        <- droplevels(factor(d$POP_ID))
  
  # 2. fit model on the cleaned data
  f   <- as.formula(paste(trait, "~ Block_Patrick + (1|POP_ID) + (1|Block_Patrick:POP_ID)"))
  mod <- lmer(f, data = d, control = lmerControl(optimizer = "bobyqa"))
  stopifnot(nrow(d) == nobs(mod))
  
  # 3. block effect from the model's own design matrix (always conformable)
  X         <- getME(mod, "X")
  block_eff <- as.vector(X %*% fixef(mod)) - fixef(mod)[["(Intercept)"]]
  
  #or
  #block_eff  <- as.vector(X %*% fixef(mod)) - fixef(mod)[["(Intercept)"]]   # per tree, block 1 = 0
  #lvl_eff    <- c(0, fixef(mod)[-1])                                          # per block level
  #block_eff  <- block_eff - mean(lvl_eff)                                    # centre on average block
  
  # 4. plot (block x provenance) BLUP for each tree
  re_plot  <- ranef(mod)$`Block_Patrick:POP_ID`
  plot_eff <- re_plot[paste(d$Block_Patrick, d$POP_ID, sep = ":"), 1]
  
  # 5. adjusted value - value used downstream
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
  
  # C. raw value + provenance BLUP - dont use this one
  d$value_old <- d[[trait]] + pop_blup
  
  

  df_combined <- d[, c("genotype","POP_ID","value_corrected","value_blup","value_old")]
  colnames(df_combined) <- c("genotype","POP_ID",
                             paste(trait, "corrected", sep="_"),
                             paste(trait, "blup", sep="_"),
                             paste(trait, "old", sep="_"))
  return(df_combined)
}




#==================#
#   MODEL FITTING  #
#==================#

height.ch <- run_model_full_v1(dset_ch, "Height")
head(height.ch)
dbh.ch <- run_model_full_v1(dset_ch, "DBH")
head(dbh.ch)
dm1 <- merge(dset_ch, height.ch, by = c("genotype", "POP_ID"), sort=F)
dm2 <- merge(dm1, dbh.ch, by = c("genotype", "POP_ID"), sort=F)
head(dm2)


height.ml <- run_model_full_v1(dset_ml, "Height")
dbh.ml <- run_model_full_v1(dset_ml, "DBH")
dm3 <- merge(dset_ml, height.ml, by = c("genotype", "POP_ID"), sort=F)
dm4 <- merge(dm3, dbh.ml, by = c("genotype", "POP_ID"), sort=F)


height.ac <- run_model_full_v1(dset_ac, "Height")
dbh.ac <- run_model_full_v1(dset_ac, "DBH")
dm5 <- merge(dset_ac, height.ac, by = c("genotype", "POP_ID"), sort=F)
dm6 <- merge(dm5, dbh.ac, by = c("genotype", "POP_ID"), sort=F)


height.pr <- run_model_full_v1(dset_pr, "Height")
dbh.pr <- run_model_full_v1(dset_pr, "DBH")
dm7 <- merge(dset_pr, height.pr, by = c("genotype", "POP_ID"), sort=F)
dm8 <- merge(dm7, dbh.pr, by = c("genotype", "POP_ID"), sort=F)



dsets <- rbind(dm2, dm4, dm6, dm8)
head(dsets)

write.table(dsets, "08_correct_Height_phenotypes_corrected.tsv", sep="\t", col.names=T, row.names = F, quote=F, append=F)





dsets



### checks
head(dsets)

p1 <- ggplot(dsets) + aes(x = Height, y = Height_corrected) +
  geom_point() +
  facet_grid(~SITE_ID)

p2 <- ggplot(dsets) + aes(x = Height, y = Height_blup) +
  geom_point() +
  facet_grid(~SITE_ID)

p3 <- ggplot(dsets) + aes(x = Height, y = Height_old) +
  geom_point() +
  facet_grid(~SITE_ID)


plot_grid(p1, p2, p3)



ggplot(dsets) + aes(x = Height_old, y = Height_corrected) +
  geom_point() +
  facet_wrap(~SITE_ID, ncol=2)


cor.test(dsets$Height_corrected, dsets$Height_old, method = "pearson")

