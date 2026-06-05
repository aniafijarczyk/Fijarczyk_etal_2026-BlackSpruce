rm(list=ls())
setwd("C:/Users/aniaf/Projects/BlackSpruce/932_gradient_forest_env")


library(gradientForest)
library(tidyverse)
library(cowplot)
library(ggrepel)
library(readxl)
library(RColorBrewer)
library(dplyr)
library(ggnewscale)


####################
#---  DATA  -------#
####################

load("gradient_02_run_results_random.RData")


###########################
#---  IMPORTANCES  -------#
###########################

# A function to extract importances from multiple files and average them
# Importances in a list are sorted alphabetically


getMeanImportances <- function(gf.list) {
  # Loop to combine importances
  imps <- as.data.frame(gf.list[[1]]$overall.imp)
  colnames(imps) <- "mcmc_1"
  for (i in 2:length(gf.list)) {
    imp <- as.data.frame(gf.list[[i]]$overall.imp)
    colnames(imp) <- paste0("mcmc_",i)
    imps <- cbind(imps, imp)
  }
  
  # Taking mean
  imps$variable <- rownames(imps)
  dg_imps <- imps %>% gather(key = "mcmc", value = "Importance",-variable) %>% 
    group_by(variable) %>% summarise(mean = mean(Importance))
  dg_imps <- as.data.frame(dg_imps) 
  rownames(dg_imps) <- dg_imps$variable
  dg_imps <- dg_imps %>% dplyr::select(-variable)
  return(dg_imps)
  
}

# A function to get a dataframe with importances 
getImportances <- function(gf.object) {
  # Loop to combine importances
  imp <- as.data.frame(gf.object$overall.imp)
  colnames(imp) <- "Importance"
  return(imp)
}

imp.random <- getImportances(gf)
head(imp.random)

dt.sorted <- imp.random %>% dplyr::arrange(-Importance)
dt.sorted
write.table(dt.sorted, "gradient_03_result_importances.tsv", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE, append=FALSE)


############### Plot most important - bars

dt <- imp.random
dt$variable <- rownames(dt)
head(dt)
subset <- dt %>% arrange(-Importance) %>% head(10)
subset

#Importance variable
#MEM1    0.0013876938     MEM1
#dPP     0.0005693769      dPP
#CMI     0.0004573133      CMI
#TP      0.0004336062       TP
#fallFP4 0.0003782788  fallFP4
#MEM3    0.0002972816     MEM3

ordered_vars <- subset %>% arrange(Importance) %>% pull(variable) %>% as.vector()
ordered_vars

#subset$variable <- factor(subset$variable, levels = ordered_vars)
subset$variable <- factor(subset$variable, levels = ordered_vars, labels=rev(c("MEM1",
                                                                               "MEM2",
                                                                               "Photoperiod >10h",
                                                                               "Climate moisture index",
                                                                               "Precipitation",
                                                                               "Fall frost prob. at -4C >10%",
                                                                               "Fall frost prob. at -2C >10%",
                                                                               "Fall frost prob. at 0C >10%",
                                                                               "MEM3",
                                                                               "Solar radiation")))


subset

p0 <- ggplot(subset) + aes(x = Importance, y = variable) +
  geom_col() + 
  labs(x = "Importance", y = "") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=18, angle=45, hjust=1),
        axis.text.y = element_text(size=22),
        axis.title = element_text(size=26),
        )
p0


png("gradient_03_results_importances_barplot_10.png", w=2000,h=1400,res=300)
p0
dev.off()


############### Plot all  - bars
dt <- imp.random
dt$variable <- rownames(dt)
head(dt)
ordered_vars <- dt %>% arrange(Importance) %>% pull(variable) %>% as.vector()
ordered_vars

#subset$variable <- factor(subset$variable, levels = ordered_vars)
dt$variable <- factor(dt$variable, levels = ordered_vars)



p0 <- ggplot(dt) + aes(x = Importance, y = variable) +
  geom_col() + 
  labs(x = "Importance", y = "") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14, angle=45, hjust=1))
p0

png("gradient_03_results_importances_barplot_all.png", w=800,h=1200,res=150)
p0
dev.off()




