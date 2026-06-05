rm(list=ls())
setwd("C:/Users/aniaf/Projects/BlackSpruce/951_garden_plot")

library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
library(RColorBrewer)




#=================#
#      DATA       #
#=================#

df1 <- read.csv("../951_garden_plot/02_bootstrap.tsv", sep="\t", header=T)
dim(df1)
df1 <- df1 %>% drop_na()
df1$test<- "GF"
df1$group<- "combined"
dim(df1)
head(df1)

df2 <- read.csv("../951_garden_plot/02B_bootstrap_rda.tsv", sep="\t", header=T)
dim(df2)
df2 <- df2 %>% drop_na()
df2$test <- "RDA"
df2$group<- "combined"
dim(df2)
head(df2)
df12 <- rbind(df1, df2)
df12 <- df12 %>% dplyr::select("SET","marker","trait","SITE_ID","group","rho","CI_low","CI_high","n_samples","n_boot","test")
head(df12)


df5 <- read.csv("../951_garden_plot/03_bootstrap_clusters.tsv", sep="\t", header=T)
dim(df5)
df5 <- df5 %>% drop_na()
df5$test<- "GF"
dim(df5)
head(df5)

df6 <- read.csv("../951_garden_plot/03B_bootstrap_clusters_rda.tsv", sep="\t", header=T)
dim(df6)
df6 <- df6 %>% drop_na()
df6$test <- "RDA"
dim(df6)
head(df6)
df56 <- rbind(df5, df6)
df_tot <- rbind(df12, df56)
df_tot$lab <- paste0(df_tot$group,"-",df_tot$test)
head(df_tot)


df_climate1 <- read.csv("../950_garden_offset_plot/12_get_ctd.tsv", sep="\t", header=T)
df_climate1$test <- "climate"
df_climate1$group <- "combined"
df_climate1$lab <- paste0(df_climate1$group,"-",df_climate1$test)
df_climate1$rho <- df_climate1$spearman_rho
df_climate1$CI_low <- NA
df_climate1$CI_high <- NA
df_climate1$n_boot <- NA
df_climate1$marker <- "1000"
df_climate1 <- df_climate1 %>% dplyr::select("SET","marker","trait","SITE_ID","group",
                                           "rho","CI_low","CI_high","n_samples","n_boot","test","lab")
head(df_climate1)

df_climate2 <- read.csv("../950_garden_offset_plot/12_get_ctd_clusters.tsv", sep="\t", header=T)
df_climate2 <- df_climate2 %>% filter(!group %in% c("ME","WI")) %>% filter(n_samples >=5)
df_climate2$test <- "climate"
df_climate2$lab <- paste0(df_climate2$group,"-",df_climate2$test)
df_climate2$rho <- df_climate2$spearman_rho
df_climate2$CI_low <- NA
df_climate2$CI_high <- NA
df_climate2$n_boot <- NA
df_climate2$marker <- "1000"
df_climate2 <- df_climate2 %>% dplyr::select("SET","marker","trait","SITE_ID","group",
                                             "rho","CI_low","CI_high","n_samples","n_boot","test","lab")
head(df_climate2)




df_tot2 <- rbind(df_tot, df_climate1, df_climate2)
head(df_tot2)

unique(df_tot2$lab)

#======================#
#        PLOTS         #
#======================#



##########################   SUBSET


df <- df_tot2 %>% dplyr::filter(marker == "1000") %>% filter(SET == "TEST")


traits <- c("Height", "Biomass_Increment","Biomass_Increment_1985","Biomass_Increment_2010", "DBH", "Average_Ring_Density", "Rs")
trait.labels <- c("Height", "ABtot", "AB5-1985","AB5-2010","DBH","Density", "Rs")

ds <- df %>% filter(trait %in% traits)
head(ds)

ds$trait_labels <- factor(ds$trait, levels = traits, labels = trait.labels)
ds$SITE_ID <- factor(ds$SITE_ID, levels = c("PR","ML","CH","AC"))
ds$model <- factor(ds$lab, levels = c("combined-climate","combined-GF","combined-RDA",
                                      "West-climate","West-GF","West-RDA",
                                      "Central-climate","Central-GF","Central-RDA",
                                      "East-climate","East-GF","East-RDA"))
unique(ds$model)
head(ds)

p2 <- ggplot(ds) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey") +
  
  geom_errorbar(data=ds, aes(x = trait_labels, ymin = CI_low, ymax = CI_high, group = model, color = model),
                position = position_dodge(width=0.7),
                width=0.01) +
  
  geom_point(data=ds, aes(x = trait_labels, y = rho, 
                          fill = model, color = model, shape=model),
             position = position_dodge(width=0.7),
             size=2.5) +
  
  scale_fill_manual(values = c("black","white","black","#DFC27D","white","#DFC27D",
                               "#C7EAE5","white","#C7EAE5","#35978F","white","#35978F"), name="Model") +
  scale_color_manual(values = c("black","black","black","#DFC27D","#DFC27D","#DFC27D",
                                "#C7EAE5","#C7EAE5","#C7EAE5","#35978F","#35978F","#35978F"), name="Model") +
  scale_shape_manual(values = c(24,21,21,24,21,21,
                                24,21,21,24,21,21), name="Model") +
  
  facet_grid(SITE_ID~.) +
  
  labs(x = "Trait", y = "Spearman rho") + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5,size=16),
        panel.background = element_rect(fill=NA, color="grey"),
        panel.grid = element_blank(),
        axis.text.y = element_text(size=10),
        axis.title = element_text(size=18),
        strip.text = element_text(size=16),
        legend.title = element_text(size=16),
        legend.text = element_text(size=12),
        legend.position = "right")
p2


png("16_plot_rho_basic_II.png", w = 3600, h = 1600, res=300)
p2
dev.off()



##########################   ALL PART I


df <- df_tot2 %>% dplyr::filter(marker == "1000") %>% filter(SET == "TEST")

traits <- c("Height", "Biomass_Increment", "DBH", "Average_Ring_Density", "Biomass_Increment_1980",
            "Biomass_Increment_1985","Biomass_Increment_1990","Biomass_Increment_1995")
trait.labels <- c("Height", "ABtot", "DBH","Density", "AB5-1980","AB5-1985","AB5-1990","AB5-1995")


ds <- df %>% filter(trait %in% traits)
head(ds)

ds$trait_labels <- factor(ds$trait, levels = traits, labels = trait.labels)
ds$SITE_ID <- factor(ds$SITE_ID, levels = c("PR","ML","CH","AC"))
ds$model <- factor(ds$lab, levels = c("combined-climate","combined-GF","combined-RDA",
                                      "West-climate","West-GF","West-RDA",
                                      "Central-climate","Central-GF","Central-RDA",
                                      "East-climate","East-GF","East-RDA"))
unique(ds$model)
head(ds)

p3 <- ggplot(ds) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey") +
  
  geom_errorbar(data=ds, aes(x = trait_labels, ymin = CI_low, ymax = CI_high, group = model, color = model),
                position = position_dodge(width=0.7),
                width=0.01) +
  
  geom_point(data=ds, aes(x = trait_labels, y = rho, 
                          fill = model, color = model, shape=model),
             position = position_dodge(width=0.7),
             size=2.5) +
  
  scale_fill_manual(values = c("black","white","black","#DFC27D","white","#DFC27D",
                               "#C7EAE5","white","#C7EAE5","#35978F","white","#35978F"), name="Model") +
  scale_color_manual(values = c("black","black","black","#DFC27D","#DFC27D","#DFC27D",
                                "#C7EAE5","#C7EAE5","#C7EAE5","#35978F","#35978F","#35978F"), name="Model") +
  scale_shape_manual(values = c(24,21,21,24,21,21,
                                24,21,21,24,21,21), name="Model") +
  
  facet_grid(SITE_ID~.) +
  
  labs(x = "Trait", y = "Spearman rho") + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5,size=16),
        panel.background = element_rect(fill=NA, color="grey"),
        panel.grid = element_blank(),
        axis.text.y = element_text(size=10),
        axis.title = element_text(size=18),
        strip.text = element_text(size=16),
        legend.title = element_text(size=16),
        legend.text = element_text(size=12),
        legend.position = "right")
p3





##########################   ALL PART II


df <- df_tot2 %>% dplyr::filter(marker == "1000") %>% filter(SET == "TEST")


traits <- c("Biomass_Increment_2000",
            "Biomass_Increment_2005", "Biomass_Increment_2010", "Biomass_Increment_2015",
            "Rc", "Rr", "Rl", "Rs")
trait.labels <- c("AB5-2000","AB5-2005", "AB5-2010", "AB5-2015", "Rc", "Rr", "Rs", "Rt")

ds <- df %>% filter(trait %in% traits)
head(ds)

ds$trait_labels <- factor(ds$trait, levels = traits, labels = trait.labels)
ds$SITE_ID <- factor(ds$SITE_ID, levels = c("PR","ML","CH","AC"))
ds$model <- factor(ds$lab, levels = c("combined-climate","combined-GF","combined-RDA",
                                      "West-climate","West-GF","West-RDA",
                                      "Central-climate","Central-GF","Central-RDA",
                                      "East-climate","East-GF","East-RDA"))
unique(ds$model)
head(ds)

p4 <- ggplot(ds) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey") +
  
  geom_errorbar(data=ds, aes(x = trait_labels, ymin = CI_low, ymax = CI_high, group = model, color = model),
                position = position_dodge(width=0.7),
                width=0.01) +
  
  geom_point(data=ds, aes(x = trait_labels, y = rho, 
                          fill = model, color = model, shape=model),
             position = position_dodge(width=0.7),
             size=2.5) +
  
  scale_fill_manual(values = c("black","white","black","#DFC27D","white","#DFC27D",
                               "#C7EAE5","white","#C7EAE5","#35978F","white","#35978F"), name="Model") +
  scale_color_manual(values = c("black","black","black","#DFC27D","#DFC27D","#DFC27D",
                                "#C7EAE5","#C7EAE5","#C7EAE5","#35978F","#35978F","#35978F"), name="Model") +
  scale_shape_manual(values = c(24,21,21,24,21,21,
                                24,21,21,24,21,21), name="Model") +
  
  facet_grid(SITE_ID~.) +
  
  labs(x = "Trait", y = "Spearman rho") + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5,size=16),
        panel.background = element_rect(fill=NA, color="grey"),
        panel.grid = element_blank(),
        axis.text.y = element_text(size=10),
        axis.title = element_text(size=18),
        strip.text = element_text(size=16),
        legend.title = element_text(size=16),
        legend.text = element_text(size=12),
        legend.position = "right")
p4





png("16_plot_rho_basic_all.png", h = 3600, w = 3400, res=300)
plot_grid(p3, p4, ncol = 1)
dev.off()

