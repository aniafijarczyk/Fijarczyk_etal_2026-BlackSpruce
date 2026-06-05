rm(list=ls())
setwd("C:/Users/aniaf/Projects/BlackSpruce/943_climate_transfer_distance")


library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape)
library(cowplot)
library(grid)
library(RColorBrewer)



#===========#
#   DATA    #
#===========#


######################
# Phenotypes - Patrick
dpheno = read.csv("../43_climate_transfer_distance/08_correct_Patrick_phenotypes_corrected.tsv", sep="\t", header=T)
head(dpheno)
dpheno$POP_ID <- as.character(dpheno$POP_ID)
dpheno$POP_SITE <- paste0(dpheno$POP_ID,"_",dpheno$SITE_ID)
head(dpheno)
colnames(dpheno)

ggplot(dpheno) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(x = Height, y = Height_corrected))

ggplot(dpheno) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(x = DBH, y = DBH_corrected))




### Mean of height
dmean <- dpheno %>% group_by(POP_ID, POP_SITE, SITE_ID) %>% summarise(Height = mean(Height_corrected), DBH = mean(DBH_corrected))
dg_mean <- dmean %>% gather(key = "Trait", value = "mean", c("Height","DBH"))
head(dg_mean)

dsd <- dpheno %>% group_by(POP_ID, POP_SITE, SITE_ID) %>% summarise(Height = sd(Height_corrected), DBH = sd(DBH_corrected))
dg_sd <- dsd %>% gather(key = "Trait", value = "sd", c("Height","DBH"))
head(dg_sd)

dn <- dpheno %>% group_by(POP_ID, POP_SITE, SITE_ID) %>% summarise(Height = sum(!is.na(Height_corrected)), DBH = sum(!is.na(DBH_corrected)))
dg_n <- dn %>% gather(key = "Trait", value = "n", c("Height","DBH"))
head(dg_n)

dcomb_ <- merge(dg_mean, dg_sd, by = c('POP_ID','POP_SITE','SITE_ID','Trait'))
dcomb <- merge(dcomb_, dg_n, by = c('POP_ID','POP_SITE','SITE_ID','Trait'), all.x=T, sort=F)
head(dcomb)



######################
# Phenotypes - Etienne
dpheno2 <- read.csv("../43_climate_transfer_distance/10_correct_Etienne_phenotypes_corrected.tsv", sep="\t", header=T)
dpheno2$POP_ID <- as.character(dpheno2$POP_ID)
dpheno2$POP_SITE <- paste0(dpheno2$POP_ID,"_",dpheno2$SITE_ID)
head(dpheno2)
unique(dpheno2$Trait)

ggplot(dpheno2[dpheno2$Trait != "Average_Ring_Density",]) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(x = log_value , y = value_corrected)) +
  facet_grid(Trait~SITE_ID, scales="free")

ggplot(dpheno2[dpheno2$Trait == "Average_Ring_Density",]) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(x = value , y = value_corrected)) +
  facet_grid(Trait~SITE_ID, scales="free")


dmean2 <- dpheno2 %>% group_by(POP_ID, POP_SITE, SITE_ID, Trait) %>% summarise(mean = mean(value_corrected))
dsd2 <- dpheno2 %>% group_by(POP_ID, POP_SITE, SITE_ID, Trait) %>% summarise(sd = sd(value_corrected))
dn2 <- dpheno2 %>% group_by(POP_ID, POP_SITE, SITE_ID, Trait) %>% summarise(n = sum(!is.na(value_corrected)))

dcomb2_ <- merge(dmean2, dsd2, by = c('POP_ID','POP_SITE','SITE_ID','Trait'))
dcomb2 <- merge(dcomb2_, dn2, by = c('POP_ID','POP_SITE','SITE_ID','Trait'), all.x=T, sort=F)
head(dcomb2)


######################
# Phenotypes - survival
dpheno3 = read.csv("../33_phenotypes/07_overview_phenotypes_Patrick_survival_data.tsv", sep="\t", header=T)
dim(dpheno3)
dpheno3 <- dpheno3[!is.na(dpheno3$POP_ID), ]
gardens3 = data.frame('E353B3' = 'CH','E353B1' = 'ML',
                      'E60A' = 'AC','G348M' = 'PR')
dgard <- as.data.frame(t(gardens3))
dgard$garden <- rownames(dgard)
colnames(dgard) <- c("SITE_ID", "garden")
head(dgard)
mpheno3 <- merge(dpheno3, dgard, by = "garden", sort=F)
mpheno3$POP_SITE <- paste0(mpheno3$POP_ID,"_",mpheno3$SITE_ID)
mpheno3$Trait <- "Survival"
mpheno3$mean <- mpheno3$Survival
mpheno3$sd <- NA
mpheno3$n <- mpheno3$Planted
head(mpheno3)

dmean3 <- mpheno3 %>% dplyr::select(POP_ID, POP_SITE, SITE_ID, Trait, mean, sd, n)
head(dmean3)


######################
# Phenotypes - Extreme
dpheno4 <- read.csv("../43_climate_transfer_distance/12_correct_Extreme_phenotypes_corrected.tsv", sep="\t", header=T)
dpheno4$POP_ID <- as.character(dpheno4$POP_ID)
dpheno4$POP_SITE <- paste0(dpheno4$POP_ID,"_",dpheno4$SITE_ID)
head(dpheno4)
unique(dpheno4$Trait)

dmean4 <- dpheno4 %>% group_by(POP_ID, POP_SITE, SITE_ID, Trait) %>% summarise(mean = mean(value_corrected))
dsd4 <- dpheno4 %>% group_by(POP_ID, POP_SITE, SITE_ID, Trait) %>% summarise(sd = sd(value_corrected))
dn4 <- dpheno4 %>% group_by(POP_ID, POP_SITE, SITE_ID, Trait) %>% summarise(n = sum(!is.na(value_corrected)))

dcomb4_ <- merge(dmean4, dsd4, by = c('POP_ID','POP_SITE','SITE_ID','Trait'))
dcomb4 <- merge(dcomb4_, dn4, by = c('POP_ID','POP_SITE','SITE_ID','Trait'), all.x=T, sort=F)
head(dcomb4)


########################

all.phenos <- rbind(dcomb, dcomb2, dmean3, dcomb4)
head(all.phenos)

#########################
# Adding strata
head(dpheno)
dstrata <- dpheno %>% dplyr::select(POP_ID, POP_GR) %>% distinct()
dstrata$POP_ID

gall.phenos <- merge(all.phenos, dstrata, by = "POP_ID", all.x=T, sort=F)
head(gall.phenos)


######################
# Climate
dclim <- read.csv("../34_climate_transfer_distance/03_climate_distance.tsv", sep="\t", header=T)
head(dclim)
colnames(dclim) <- c("POP_ID","AC","CH","ML","PR")
gclim <- dclim %>% gather(key="SITE_ID", value = "EuclDist", -POP_ID)
gclim$POP_SITE <- paste0(gclim$POP_ID,"_",gclim$SITE_ID)
head(gclim)

#dclim2 <- read.csv("../34_climate_transfer_distance/03_climate_distance_All_vars.tsv", sep="\t", header=T)
#head(dclim2)
#unique(dclim2$CommonGarden)
#dim(dclim)
#dclim <- dclim %>% dplyr::rename(POP_ID = Provenance)
#dclim$POP_SITE <- paste0(dclim$POP_ID,"_",dclim$SITE_ID)
#head(dclim)
#length(unique(dclim$POP_SITE))

################
# Combining with climate
clim.phenos <- merge(gall.phenos, gclim, by = c("POP_ID", "SITE_ID", "POP_SITE"), sort=F, all.x=T)
dim(clim.phenos)
head(clim.phenos)
length(unique(clim.phenos$POP_SITE))

#####################
# Filtering for number of ind per group ???
flim.phenos <- clim.phenos %>% filter(n>4)
dim(flim.phenos)

##################
# Saving dataset

write.table(clim.phenos, file = "13_plot_climate_distance_input.tsv", sep="\t", row.names = F, col.names=T, quote=F, append=F)



