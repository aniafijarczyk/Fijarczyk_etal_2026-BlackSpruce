rm(list=ls())
setwd("C:/Users/aniaf/Projects/BlackSpruce/950_garden_offset_plot")

library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
library(RColorBrewer)



#=================#
#      DATA       #
#=================#


trait <- "Biomass_Increment"

dclim <- read.csv("../34_climate_transfer_distance/03_climate_distance.tsv", sep="\t", header=T)
colnames(dclim) <- c("POP_ID","AC","CH","ML","PR")
gclim <- dclim %>% gather(key="SITE_ID", value = "EuclDist", -POP_ID)
gclim$POP_SITE <- paste0(gclim$POP_ID,"_",gclim$SITE_ID)
head(gclim)
gclim <- gclim %>% dplyr::select(POP_SITE, EuclDist)
head(gclim)
length(unique(gclim$POP_SITE))


### phenotypes
dpheno <- read.csv("../44_train_and_test/00_combine_samples_means.tsv", sep="\t", header=T)
dpheno <- dpheno %>% filter(SET == 'TEST')
dim(dpheno)
dpheno$sample <- dpheno$POP_SITE
head(dpheno)
mpheno <- merge(dpheno, gclim, by = "POP_SITE", sort=F)
head(mpheno)

dpheno_sub <- mpheno %>% filter((Trait_name == "Biomass_Increment") & (SITE_ID == "PR"))
dim(dpheno_sub)
dpheno_sub %>% head()

### genetic offsets
df <- read.csv("../948_garden_offset_run_separate/results/01_run_gradient_forest_PR_1000_Biomass_Increment.tsv", header=T, sep="\t")
dim(df)
df %>% head()

### combining
dm <- merge(df, dpheno_sub, by = c("Trait_name","SITE_ID","sample"), sort=F, all.x = T)
dim(dm)
dm %>% head()


### Filtering small clusters
unique(dm$group)
dm <- dm %>% filter(group != "East")

############ Plot ctd

dm$group <- factor(dm$group, levels = c("West","Central","East"))
p1 <- ggplot(dm) + aes(x = EuclDist, y = mean) +
  geom_point(aes(fill = group), pch=21, size=3) +
  scale_fill_manual(values = c("#DFC27D","#C7EAE5", "black"), name="Cluster") +
  geom_smooth(method = "lm", se = FALSE, aes(color=group, linetype = group)) +
  geom_smooth(method = "lm", se = FALSE, aes(color = "Combined", linetype = "Combined")) +
  scale_color_manual(values = c("#DFC27D","#C7EAE5", "black"), name="Cluster") +
  scale_linetype_manual(name = "Cluster", values = c("solid","solid","dashed")) +
  labs(x = "Climate transfer\ndistance", y = "log ABtot") +
  
  geom_smooth(method = "lm", se = FALSE, color="black", linetype="dashed") +
  theme(
    panel.background = element_rect(fill=NA, color="grey"),
    panel.grid = element_blank(),
    axis.text = element_text(size=12),
    axis.title = element_text(size=16)
  )
p1


############ Plot offset


p2 <- ggplot(dm) + aes(x = offset, y = mean) +
  geom_point(aes(fill = group), pch=21, size=3) +
  scale_fill_manual(values = c("#DFC27D","#C7EAE5", "black"), name="Cluster") +
  geom_smooth(method = "lm", se = FALSE, aes(color=group, linetype = group)) +
  geom_smooth(method = "lm", se = FALSE, aes(color = "Combined", linetype = "Combined")) +
  scale_color_manual(values = c("#DFC27D","#C7EAE5", "black"), name="Cluster") +
  scale_linetype_manual(name = "Cluster", values = c("solid","solid","dashed")) +
  labs(x = "Genomic offset (GF)\n(1000 SNPs)", y = "log ABtot") +
  
  geom_smooth(method = "lm", se = FALSE, color="black", linetype="dashed") +
  theme(
    panel.background = element_rect(fill=NA, color="grey"),
    panel.grid = element_blank(),
    axis.text = element_text(size=12),
    axis.title = element_text(size=16)
  )
p2


############ Plot offset for lfmm


### genetic offsets
df3 <- read.csv("../948_garden_offset_run_separate/results/01_run_gradient_forest_PR_lfmm_Biomass_Increment.tsv", header=T, sep="\t")
dim(df3)
df3 %>% head()
dm3 <- merge(df3, dpheno_sub, by = c("Trait_name","SITE_ID","sample"), sort=F, all.x = T)
dim(dm3)
dm3 %>% head()
unique(dm3$group)
dm3$group <- factor(dm3$group, levels = c("West","Central","East"))
dm3 <- dm3 %>% filter(group != "East")


p3 <- ggplot(dm3) + aes(x = offset, y = mean) +
  geom_point(aes(fill = group), pch=21, size=3) +
  scale_fill_manual(values = c("#DFC27D","#C7EAE5", "black"), name="Cluster") +
  geom_smooth(method = "lm", se = FALSE, aes(color=group, linetype = group)) +
  geom_smooth(method = "lm", se = FALSE, aes(color = "Combined", linetype = "Combined")) +
  scale_color_manual(values = c("#DFC27D","#C7EAE5", "black"), name="Cluster") +
  scale_linetype_manual(name = "Cluster", values = c("solid","solid","dashed")) +
  labs(x = "Genomic offset (GF)\n(LFMM outliers)", y = "log ABtot") +
  
  geom_smooth(method = "lm", se = FALSE, color="black", linetype="dashed") +
  theme(
    panel.background = element_rect(fill=NA, color="grey"),
    panel.grid = element_blank(),
    axis.text = element_text(size=12),
    axis.title = element_text(size=16)
  )
p3


############ Plot offset for rda


### genetic offsets
df4 <- read.csv("../991_run_rda_separate/results/01_run_gradient_forest_PR_lfmm_Biomass_Increment.tsv", header=T, sep="\t")
dim(df4)
df4 %>% head()
dm4 <- merge(df4, dpheno_sub, by = c("Trait_name","SITE_ID","sample"), sort=F, all.x = T)
dim(dm4)
dm4 %>% head()
unique(dm4$group)
dm4$group <- factor(dm4$group, levels = c("West","Central","East"))
dm4 <- dm4 %>% filter(group != "East")
head(dm4)

p4 <- ggplot(dm4) + aes(x = offset_rda, y = mean) +
  geom_point(aes(fill = group), pch=21, size=3) +
  scale_fill_manual(values = c("#DFC27D","#C7EAE5", "black"), name="Cluster") +
  geom_smooth(method = "lm", se = FALSE, aes(color=group, linetype = group)) +
  geom_smooth(method = "lm", se = FALSE, aes(color = "Combined", linetype = "Combined")) +
  scale_color_manual(values = c("#DFC27D","#C7EAE5", "black"), name="Cluster") +
  scale_linetype_manual(name = "Cluster", values = c("solid","solid","dashed")) +
  
  labs(x = "Genomic offset (RDA)\n(LFMM outliers)", y = "log ABtot") +
  
  #geom_smooth(method = "lm", se = FALSE, color="black", linetype="dashed") +

  theme(
    panel.background = element_rect(fill=NA, color="grey"),
    panel.grid = element_blank(),
    axis.text = element_text(size=12),
    axis.title = element_text(size=16)
  )
p4



plot_grid(p1, p2, p3, p4, ncol=2)

png("19_plot_selected.png", w = 2500, h=2000, res=300)
plot_grid(p1, p2, p3, p4, ncol=2)
dev.off()




