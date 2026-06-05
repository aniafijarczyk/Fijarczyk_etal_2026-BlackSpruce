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
dim(df1)
head(df1)

df2 <- read.csv("../951_garden_plot/02B_bootstrap_rda.tsv", sep="\t", header=T)
dim(df2)
df2 <- df2 %>% drop_na()
df2$test <- "RDA"
dim(df2)
head(df2)
df12 <- rbind(df1, df2)


df3 <- read.csv("../950_garden_offset_plot/12_get_ctd.tsv", sep="\t", header=T)
df3$test <- "GF"
df4 <- df3
df4$test <- "RDA"
df34 <- rbind(df3, df4)
head(df34)




#======================#
#     ALTERNATIVE      #
#======================#


df <- df12 %>% filter(SET == "TEST")
df$marker %>% unique()

traits <- c("Height", "Biomass_Increment")
trait.labels <- c("Height", "ABtot")
ds <- df %>% filter(trait %in% traits) %>% filter(marker %in% c("100","100LF","1000","1000LF","all","lfmm","RDA","RDAcorrected"))
head(ds)

ds$test <- factor(ds$test, levels = c("GF","RDA"))
ds$trait_labels <- factor(ds$trait, levels = traits, labels = trait.labels)
#ds$SET <- factor(ds$SET, levels = c("TRAIN","TEST"))
ds$SITE_ID <- factor(ds$SITE_ID, levels = c("PR","ML","CH","AC"))
ds$marker_labels <- factor(ds$marker, levels = c("100","100LF","1000","1000LF","all","lfmm","RDA","RDAcorrected"),
                           labels = c("0.1 K", "0.1 K-lf", "1 K", "1 K-lf", "~29 K", "LFMM", "RDA", "RDA-struct"))


df3_sub <- df34 %>% filter(SET == "TEST") %>% filter(trait %in% traits)
df3_sub$trait_labels <- factor(df3_sub$trait, levels = traits, labels = trait.labels)
#df3_sub$SET <- factor(df3_sub$SET, levels = c("TRAIN","TEST"))
df3_sub$SITE_ID <- factor(df3_sub$SITE_ID, levels = c("PR","ML","CH","AC"))
df3_sub$test <- factor(df3_sub$test, levels = c("GF","RDA"))


##################################################################################

df3_sub <- df3_sub %>%
  mutate(x_start = as.numeric(factor(test)) - 0.5,   # start a bit to the left
         x_end   = as.numeric(factor(test)) + 0.5)   # end a bit to the right)



p1 <- ggplot(df3_sub) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey") +
  
  geom_errorbar(data=ds, aes(x = test, ymin = CI_low, ymax = CI_high, group = marker_labels,
                             color = marker_labels),
                position = position_dodge(width=0.7),
                width=0.01) +
  
  geom_segment(data=df3_sub,  aes(x = x_start, xend = x_end, y = spearman_rho, yend = spearman_rho),
               linetype="dotted", linewidth=1,color="black") +
  
  geom_point(data=ds, aes(x = test, y = rho, shape = marker_labels,
                           fill = marker_labels, color = marker_labels),
             position = position_dodge(width=0.7),
             size=3) +
  
  
  
  scale_shape_manual(values = c(22,22,22,22,22,10,8,12), name="Markers") +
  scale_fill_manual(values = c("#0072b2","white","#e69f00","white","black","black","black","black"), name="Markers") +
  scale_color_manual(values = c("#0072b2","#0072b2","#e69f00","#e69f00","black","black","black","black"), name="Markers") +
  
  facet_grid(trait_labels~SITE_ID) +
  
  labs(x = "Dataset", y = "Spearman rho") + 
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

png("08_plot_correlations_model.png", w = 3600, h = 1600, res=300)
p1
dev.off()





##########################   ALL part I


df <- df12 %>% filter(SET == "TEST")
df$marker %>% unique()

traits <- c("Height", "Biomass_Increment", "DBH", "Average_Ring_Density", "Biomass_Increment_1980",
            "Biomass_Increment_1985","Biomass_Increment_1990","Biomass_Increment_1995")
trait.labels <- c("Height", "ABtot", "DBH","Density", "AB5-1980","AB5-1985","AB5-1990","AB5-1995")

ds <- df %>% filter(trait %in% traits) %>% filter(marker %in% c("100","100LF","1000","1000LF","all","lfmm","RDA","RDAcorrected"))
head(ds)

ds$test <- factor(ds$test, levels = c("GF","RDA"))
ds$trait_labels <- factor(ds$trait, levels = traits, labels = trait.labels)
#ds$SET <- factor(ds$SET, levels = c("TRAIN","TEST"))
ds$SITE_ID <- factor(ds$SITE_ID, levels = c("PR","ML","CH","AC"))
ds$marker_labels <- factor(ds$marker, levels = c("100","100LF","1000","1000LF","all","lfmm","RDA","RDAcorrected"),
                           labels = c("0.1 K", "0.1 K-lf", "1 K", "1 K-lf", "~29 K", "LFMM", "RDA", "RDA-struct"))


df3_sub <- df34 %>% filter(SET == "TEST") %>% filter(trait %in% traits)
df3_sub$trait_labels <- factor(df3_sub$trait, levels = traits, labels = trait.labels)
#df3_sub$SET <- factor(df3_sub$SET, levels = c("TRAIN","TEST"))
df3_sub$SITE_ID <- factor(df3_sub$SITE_ID, levels = c("PR","ML","CH","AC"))
df3_sub$test <- factor(df3_sub$test, levels = c("GF","RDA"))


df3_sub <- df3_sub %>%
  mutate(x_start = as.numeric(factor(test)) - 0.5,   # start a bit to the left
         x_end   = as.numeric(factor(test)) + 0.5)   # end a bit to the right)
head(df3_sub)





p2 <- ggplot(df3_sub) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey") +
  
  geom_errorbar(data=ds, aes(x = test, ymin = CI_low, ymax = CI_high, group = marker_labels,
                               color = marker_labels),
                position = position_dodge(width=0.9),
                width=0.01) +
  
  geom_segment(data=df3_sub,  aes(x = x_start, xend = x_end, y = spearman_rho, yend = spearman_rho),
               linetype="solid", linewidth=0.5,color="black") +
  
  geom_point(data=ds, aes(x = test, y = rho, shape = marker_labels,
                            fill = marker_labels, color = marker_labels),
             position = position_dodge(width=0.9),
             size=2) +
  
  scale_shape_manual(values = c(22,22,22,22,22,10,8,12), name="Markers") +
  scale_fill_manual(values = c("#0072b2","white","#e69f00","white","black","black","black","black"), name="Markers") +
  scale_color_manual(values = c("#0072b2","#0072b2","#e69f00","#e69f00","black","black","black","black"), name="Markers") +
  
  facet_grid(SITE_ID~trait_labels) +
  
  labs(x = "Cluster", y = "Spearman rho") + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        panel.background = element_rect(fill=NA, color="grey"),
        panel.grid = element_blank(),
        axis.text = element_text(size=12),
        axis.title = element_text(size=16),
        strip.text = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.position = "bottom")
p2


##########################   ALL PART II



df <- df12 %>% filter(SET == "TEST")
df$marker %>% unique()

traits <- c("Biomass_Increment_2000",
             "Biomass_Increment_2005", "Biomass_Increment_2010", "Biomass_Increment_2015",
             "Rc", "Rr", "Rl", "Rs")
trait.labels <- c("AB5-2000",
                   "AB5-2005", "AB5-2010", "AB5-2015", "Rc", "Rr", "Rs", "Rt")

ds <- df %>% filter(trait %in% traits) %>% filter(marker %in% c("100","100LF","1000","1000LF","all","lfmm","RDA","RDAcorrected"))
head(ds)

ds$test <- factor(ds$test, levels = c("GF","RDA"))
ds$trait_labels <- factor(ds$trait, levels = traits, labels = trait.labels)
#ds$SET <- factor(ds$SET, levels = c("TRAIN","TEST"))
ds$SITE_ID <- factor(ds$SITE_ID, levels = c("PR","ML","CH","AC"))
ds$marker_labels <- factor(ds$marker, levels = c("100","100LF","1000","1000LF","all","lfmm","RDA","RDAcorrected"),
                           labels = c("0.1 K", "0.1 K-lf", "1 K", "1 K-lf", "~29 K", "LFMM", "RDA", "RDA-struct"))


df3_sub <- df34 %>% filter(SET == "TEST") %>% filter(trait %in% traits)
df3_sub$trait_labels <- factor(df3_sub$trait, levels = traits, labels = trait.labels)
#df3_sub$SET <- factor(df3_sub$SET, levels = c("TRAIN","TEST"))
df3_sub$SITE_ID <- factor(df3_sub$SITE_ID, levels = c("PR","ML","CH","AC"))
df3_sub$test <- factor(df3_sub$test, levels = c("GF","RDA"))


df3_sub <- df3_sub %>%
  mutate(x_start = as.numeric(factor(test)) - 0.5,   # start a bit to the left
         x_end   = as.numeric(factor(test)) + 0.5)   # end a bit to the right)
head(df3_sub)




p3 <- ggplot(df3_sub) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey") +
  
  geom_errorbar(data=ds, aes(x = test, ymin = CI_low, ymax = CI_high, group = marker_labels,
                               color = marker_labels),
                position = position_dodge(width=0.9),
                width=0.01) +
  
  geom_segment(data=df3_sub,  aes(x = x_start, xend = x_end, y = spearman_rho, yend = spearman_rho),
               linetype="solid", linewidth=0.5,color="black") +
  
  geom_point(data=ds, aes(x = test, y = rho, shape = marker_labels,
                            fill = marker_labels, color = marker_labels),
             position = position_dodge(width=0.9),
             size=2) +
  
  
  
  scale_shape_manual(values = c(22,22,22,22,22,10,8,12), name="Markers") +
  scale_fill_manual(values = c("#0072b2","white","#e69f00","white","black","black","black","black"), name="Markers") +
  scale_color_manual(values = c("#0072b2","#0072b2","#e69f00","#e69f00","black","black","black","black"), name="Markers") +
  
  facet_grid(SITE_ID~trait_labels) +
  
  labs(x = "Cluster", y = "Spearman rho") + 
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


png("08_plot_correlations_model_all.png", h = 3600, w = 3400, res=300)
plot_grid(p2, p3, ncol = 1)
dev.off()

