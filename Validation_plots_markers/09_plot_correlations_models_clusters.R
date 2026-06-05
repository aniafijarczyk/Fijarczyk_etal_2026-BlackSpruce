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

df1 <- read.csv("../951_garden_plot/03_bootstrap_clusters.tsv", sep="\t", header=T)
dim(df1)
df1 <- df1 %>% drop_na()
df1$test<- "GF"
dim(df1)
head(df1)

df2 <- read.csv("../951_garden_plot/03B_bootstrap_clusters_rda.tsv", sep="\t", header=T)
dim(df2)
df2 <- df2 %>% drop_na()
df2$test <- "RDA"
dim(df2)
head(df2)
df12 <- rbind(df1, df2)


df3 <- read.csv("../950_garden_offset_plot/12_get_ctd_clusters.tsv", sep="\t", header=T)
df3 <- df3 %>% filter(!group %in% c("ME","WI")) %>% filter(n_samples >=5)
df3$test <- "GF"
df4 <- df3
df4$test <- "RDA"
df34 <- rbind(df3, df4)
head(df34)




#======================#
#     ALTERNATIVE      #
#======================#

df <- rbind(df12) %>% filter(SET == "TEST")

traits <- c("Height", "Biomass_Increment")
trait.labels <- c("Height", "ABtot")
ds <- df %>% filter(trait %in% traits) %>% filter(marker %in% c("100","100LF","1000","1000LF","all","lfmm","RDA","RDAcorrected"))
head(ds)

ds$test <- factor(ds$test, levels = c("GF","RDA"))
ds$trait_labels <- factor(ds$trait, levels = traits, labels = trait.labels)
#ds$SET <- factor(ds$SET, levels = c("TEST","TRAIN"))
ds$SITE_ID <- factor(ds$SITE_ID, levels = c("PR","ML","CH","AC"))
ds$group <- factor(ds$group, levels = c("West","Central","East"))
ds$marker_labels <- factor(ds$marker, levels = c("100","100LF","1000","1000LF","all","lfmm","RDA","RDAcorrected"),
                           labels = c("0.1 K", "0.1 K-lf", "1 K", "1 K-lf", "~29 K", "LFMM", "RDA", "RDA-struct"))

df3_sub <- df34 %>% filter(SET == "TEST") %>% filter(trait %in% traits)
df3_sub$test <- factor(df3_sub$test, levels = c("GF","RDA"))
df3_sub$trait_labels <- factor(df3_sub$trait, levels = traits, labels = trait.labels)
#df3_sub$SET <- factor(df3_sub$SET, levels = c("TEST","TRAIN"))
df3_sub$SITE_ID <- factor(df3_sub$SITE_ID, levels = c("PR","ML","CH","AC"))
df3_sub$group <- factor(df3_sub$group, levels = c("West","Central","East"))

df.a <- ds %>% filter(test == "GF")
df3.a <- df3_sub %>% filter(test == "GF")
head(df.a)
head(df3.a)

df3.a <- df3.a %>%
  mutate(x_start = as.numeric(factor(group)) - 0.5,   # start a bit to the left
         x_end   = as.numeric(factor(group)) + 0.5)   # end a bit to the right)
head(df3.a)


p1 <- ggplot(df3.a) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey") +
  
  geom_errorbar(data=df.a, aes(x = group, ymin = CI_low, ymax = CI_high, group = marker_labels,
                             color = marker_labels),
                position = position_dodge(width=0.9),
                width=0.01) +
  
  geom_segment(data=df3.a,  aes(x = x_start, xend = x_end, y = spearman_rho, yend = spearman_rho),
               linetype="dotted", linewidth=1,color="black") +
  
  geom_point(data=df.a, aes(x = group, y = rho, shape = marker_labels,
                          fill = marker_labels, color = marker_labels),
             position = position_dodge(width=0.9),
             size=3) +
  
  scale_shape_manual(values = c(22,22,22,22,22,10,8,12), name="Markers") +
  scale_fill_manual(values = c("#0072b2","white","#e69f00","white","black","black","black","black"), name="Markers") +
  scale_color_manual(values = c("#0072b2","#0072b2","#e69f00","#e69f00","black","black","black","black"), name="Markers") +
  
  facet_grid(trait_labels~SITE_ID) +
  
  labs(x = "Cluster", y = "Spearman rho") + 
  ggtitle("Gradient Forest") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        panel.background = element_rect(fill=NA, color="grey"),
        panel.grid = element_blank(),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        strip.text = element_text(size=16),
        legend.title = element_text(size=16),
        legend.text = element_text(size=12),
        legend.position = "right")
p1

png("09_plot_correlations_models_clusters_GF.png", w = 3600, h = 1600, res=300)
p1
dev.off()



###################################################################################################

df.a <- ds %>% filter(test == "RDA")
df3.a <- df3_sub %>% filter(test == "RDA")
head(df.a)
head(df3.a)

df3.a <- df3.a %>%
  mutate(x_start = as.numeric(factor(group)) - 0.5,   # start a bit to the left
         x_end   = as.numeric(factor(group)) + 0.5)   # end a bit to the right)
head(df3.a)


p2 <- ggplot(df3.a) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey") +
  
  geom_errorbar(data=df.a, aes(x = group, ymin = CI_low, ymax = CI_high, group = marker_labels,
                               color = marker_labels),
                position = position_dodge(width=0.9),
                width=0.01) +
  
  geom_segment(data=df3.a,  aes(x = x_start, xend = x_end, y = spearman_rho, yend = spearman_rho),
               linetype="dotted", linewidth=1,color="black") +
  
  geom_point(data=df.a, aes(x = group, y = rho, shape = marker_labels,
                            fill = marker_labels, color = marker_labels),
             position = position_dodge(width=0.9),
             size=3) +
  
  
  
  scale_shape_manual(values = c(22,22,22,22,22,10,8,12), name="Markers") +
  scale_fill_manual(values = c("#0072b2","white","#e69f00","white","black","black","black","black"), name="Markers") +
  scale_color_manual(values = c("#0072b2","#0072b2","#e69f00","#e69f00","black","black","black","black"), name="Markers") +
  
  facet_grid(trait_labels~SITE_ID) +
  
  labs(x = "Cluster", y = "Spearman rho") + 
  ggtitle("RDA") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        panel.background = element_rect(fill=NA, color="grey"),
        panel.grid = element_blank(),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        strip.text = element_text(size=16),
        legend.title = element_text(size=16),
        legend.text = element_text(size=12),
        legend.position = "right")
p2

png("09_plot_correlations_models_clusters_RDA.png", w = 3600, h = 1600, res=300)
p2
dev.off()


###################################################################################################

##########################   ALL part I - GF


df <- rbind(df12) %>% filter(SET == "TEST")

traits <- c("Height", "Biomass_Increment", "DBH", "Average_Ring_Density", "Biomass_Increment_1980",
            "Biomass_Increment_1985","Biomass_Increment_1990","Biomass_Increment_1995")
trait.labels <- c("Height", "ABtot", "DBH","Density", "AB5-1980","AB5-1985","AB5-1990","AB5-1995")

ds <- df %>% filter(trait %in% traits) %>% filter(marker %in% c("100","100LF","1000","1000LF","all","lfmm","RDA","RDAcorrected"))
head(ds)

ds$test <- factor(ds$test, levels = c("GF","RDA"))
ds$trait_labels <- factor(ds$trait, levels = traits, labels = trait.labels)
#ds$SET <- factor(ds$SET, levels = c("TEST","TRAIN"))
ds$SITE_ID <- factor(ds$SITE_ID, levels = c("PR","ML","CH","AC"))
ds$group <- factor(ds$group, levels = c("West","Central","East"))
ds$marker_labels <- factor(ds$marker, levels = c("100","100LF","1000","1000LF","all","lfmm","RDA","RDAcorrected"),
                           labels = c("0.1 K", "0.1 K-lf", "1 K", "1 K-lf", "~29 K", "LFMM", "RDA", "RDA-struct"))

df3_sub <- df34 %>% filter(SET == "TEST") %>% filter(trait %in% traits)
df3_sub$test <- factor(df3_sub$test, levels = c("GF","RDA"))
df3_sub$trait_labels <- factor(df3_sub$trait, levels = traits, labels = trait.labels)
#df3_sub$SET <- factor(df3_sub$SET, levels = c("TEST","TRAIN"))
df3_sub$SITE_ID <- factor(df3_sub$SITE_ID, levels = c("PR","ML","CH","AC"))
df3_sub$group <- factor(df3_sub$group, levels = c("West","Central","East"))

df.a <- ds %>% filter(test == "GF")
df3.a <- df3_sub %>% filter(test == "GF")
head(df.a)
head(df3.a)

df3.a <- df3.a %>%
  mutate(x_start = as.numeric(factor(group)) - 0.5,   # start a bit to the left
         x_end   = as.numeric(factor(group)) + 0.5)   # end a bit to the right)
head(df3.a)


p3 <- ggplot(df3.a) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey") +
  
  geom_errorbar(data=df.a, aes(x = group, ymin = CI_low, ymax = CI_high, group = marker_labels,
                             color = marker_labels),
                position = position_dodge(width=0.9),
                width=0.01) +
  
  geom_segment(data=df3.a,  aes(x = x_start, xend = x_end, y = spearman_rho, yend = spearman_rho),
               linetype="solid", linewidth=0.5,color="black") +
  
  geom_point(data=df.a, aes(x = group, y = rho, shape = marker_labels,
                          fill = marker_labels, color = marker_labels),
             position = position_dodge(width=0.9),
             size=2) +
  
  scale_shape_manual(values = c(22,22,22,22,22,10,8,12), name="Markers") +
  scale_fill_manual(values = c("#0072b2","white","#e69f00","white","black","black","black","black"), name="Markers") +
  scale_color_manual(values = c("#0072b2","#0072b2","#e69f00","#e69f00","black","black","black","black"), name="Markers") +
  
  facet_grid(SITE_ID~trait_labels) +
  
  labs(x = "Cluster", y = "Spearman rho") + 
  ggtitle("Gradient Forest") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        panel.background = element_rect(fill=NA, color="grey"),
        panel.grid = element_blank(),
        axis.text = element_text(size=12),
        axis.title = element_text(size=16),
        strip.text = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.position = "bottom")
p3



##########################   ALL PART I RDA

df.a <- ds %>% filter(test == "RDA")
df3.a <- df3_sub %>% filter(test == "RDA")
head(df.a)
head(df3.a)

df3.a <- df3.a %>%
  mutate(x_start = as.numeric(factor(group)) - 0.5,   # start a bit to the left
         x_end   = as.numeric(factor(group)) + 0.5)   # end a bit to the right)
head(df3.a)


p4 <- ggplot(df3.a) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey") +
  
  geom_errorbar(data=df.a, aes(x = group, ymin = CI_low, ymax = CI_high, group = marker_labels,
                               color = marker_labels),
                position = position_dodge(width=0.9),
                width=0.01) +
  
  geom_segment(data=df3.a,  aes(x = x_start, xend = x_end, y = spearman_rho, yend = spearman_rho),
               linetype="solid", linewidth=0.5,color="black") +
  
  geom_point(data=df.a, aes(x = group, y = rho, shape = marker_labels,
                            fill = marker_labels, color = marker_labels),
             position = position_dodge(width=0.9),
             size=2) +
  
  scale_shape_manual(values = c(22,22,22,22,22,10,8,12), name="Markers") +
  scale_fill_manual(values = c("#0072b2","white","#e69f00","white","black","black","black","black"), name="Markers") +
  scale_color_manual(values = c("#0072b2","#0072b2","#e69f00","#e69f00","black","black","black","black"), name="Markers") +
  
  facet_grid(SITE_ID~trait_labels) +
  
  labs(x = "Cluster", y = "Spearman rho") + 
  ggtitle("RDA") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        panel.background = element_rect(fill=NA, color="grey"),
        panel.grid = element_blank(),
        axis.text = element_text(size=12),
        axis.title = element_text(size=16),
        strip.text = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.position = "bottom")
p4



png("09_plot_correlations_models_clusters_partI.png", h = 3600, w = 3400, res=300)
plot_grid(p3, p4, ncol = 1)
dev.off()








###################################################################################################

##########################   ALL part II - GF


df <- rbind(df12) %>% filter(SET == "TEST")

traits <- c("Biomass_Increment_2000",
             "Biomass_Increment_2005", "Biomass_Increment_2010", "Biomass_Increment_2015",
             "Rc", "Rr", "Rl", "Rs")
trait.labels <- c("AB5-2000",
                   "AB5-2005", "AB5-2010", "AB5-2015", "Rc", "Rr", "Rs", "Rt")

ds <- df %>% filter(trait %in% traits) %>% filter(marker %in% c("100","100LF","1000","1000LF","all","lfmm","RDA","RDAcorrected"))
head(ds)

ds$test <- factor(ds$test, levels = c("GF","RDA"))
ds$trait_labels <- factor(ds$trait, levels = traits, labels = trait.labels)
#ds$SET <- factor(ds$SET, levels = c("TEST","TRAIN"))
ds$SITE_ID <- factor(ds$SITE_ID, levels = c("PR","ML","CH","AC"))
ds$group <- factor(ds$group, levels = c("West","Central","East"))
ds$marker_labels <- factor(ds$marker, levels = c("100","100LF","1000","1000LF","all","lfmm","RDA","RDAcorrected"),
                           labels = c("0.1 K", "0.1 K-lf", "1 K", "1 K-lf", "~29 K", "LFMM", "RDA", "RDA-struct"))

df3_sub <- df34 %>% filter(SET == "TEST") %>% filter(trait %in% traits)
df3_sub$test <- factor(df3_sub$test, levels = c("GF","RDA"))
df3_sub$trait_labels <- factor(df3_sub$trait, levels = traits, labels = trait.labels)
#df3_sub$SET <- factor(df3_sub$SET, levels = c("TEST","TRAIN"))
df3_sub$SITE_ID <- factor(df3_sub$SITE_ID, levels = c("PR","ML","CH","AC"))
df3_sub$group <- factor(df3_sub$group, levels = c("West","Central","East"))

df.a <- ds %>% filter(test == "GF")
df3.a <- df3_sub %>% filter(test == "GF")
head(df.a)
head(df3.a)

df3.a <- df3.a %>%
  mutate(x_start = as.numeric(factor(group)) - 0.5,   # start a bit to the left
         x_end   = as.numeric(factor(group)) + 0.5)   # end a bit to the right)
head(df3.a)


p5 <- ggplot(df3.a) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey") +
  
  geom_errorbar(data=df.a, aes(x = group, ymin = CI_low, ymax = CI_high, group = marker_labels,
                               color = marker_labels),
                position = position_dodge(width=0.9),
                width=0.01) +
  
  geom_segment(data=df3.a,  aes(x = x_start, xend = x_end, y = spearman_rho, yend = spearman_rho),
               linetype="solid", linewidth=0.5,color="black") +
  
  geom_point(data=df.a, aes(x = group, y = rho, shape = marker_labels,
                            fill = marker_labels, color = marker_labels),
             position = position_dodge(width=0.9),
             size=2) +
  
  scale_shape_manual(values = c(22,22,22,22,22,10,8,12), name="Markers") +
  scale_fill_manual(values = c("#0072b2","white","#e69f00","white","black","black","black","black"), name="Markers") +
  scale_color_manual(values = c("#0072b2","#0072b2","#e69f00","#e69f00","black","black","black","black"), name="Markers") +
  
  facet_grid(SITE_ID~trait_labels) +
  
  labs(x = "Cluster", y = "Spearman rho") + 
  ggtitle("Gradient Forest") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        panel.background = element_rect(fill=NA, color="grey"),
        panel.grid = element_blank(),
        axis.text = element_text(size=12),
        axis.title = element_text(size=16),
        strip.text = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.position = "bottom")
p5



##########################   ALL PART II RDA

df.a <- ds %>% filter(test == "RDA")
df3.a <- df3_sub %>% filter(test == "RDA")
head(df.a)
head(df3.a)

df3.a <- df3.a %>%
  mutate(x_start = as.numeric(factor(group)) - 0.5,   # start a bit to the left
         x_end   = as.numeric(factor(group)) + 0.5)   # end a bit to the right)
head(df3.a)


p6 <- ggplot(df3.a) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey") +
  
  geom_errorbar(data=df.a, aes(x = group, ymin = CI_low, ymax = CI_high, group = marker_labels,
                               color = marker_labels),
                position = position_dodge(width=0.9),
                width=0.01) +
  
  geom_segment(data=df3.a,  aes(x = x_start, xend = x_end, y = spearman_rho, yend = spearman_rho),
               linetype="solid", linewidth=0.5,color="black") +
  
  geom_point(data=df.a, aes(x = group, y = rho, shape = marker_labels,
                            fill = marker_labels, color = marker_labels),
             position = position_dodge(width=0.9),
             size=2) +
  
  scale_shape_manual(values = c(22,22,22,22,22,10,8,12), name="Markers") +
  scale_fill_manual(values = c("#0072b2","white","#e69f00","white","black","black","black","black"), name="Markers") +
  scale_color_manual(values = c("#0072b2","#0072b2","#e69f00","#e69f00","black","black","black","black"), name="Markers") +
  
  facet_grid(SITE_ID~trait_labels) +
  
  labs(x = "Cluster", y = "Spearman rho") + 
  ggtitle("RDA") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        panel.background = element_rect(fill=NA, color="grey"),
        panel.grid = element_blank(),
        axis.text = element_text(size=12),
        axis.title = element_text(size=16),
        strip.text = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.position = "bottom")
p6



png("09_plot_correlations_models_clusters_partII.png", h = 3600, w = 3400, res=300)
plot_grid(p5, p6, ncol = 1)
dev.off()


