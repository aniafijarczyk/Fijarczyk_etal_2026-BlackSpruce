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

df <- read.csv("12_compare_offsets.tsv", sep="\t", header=T)
dim(df)
df <- df %>% drop_na()
dim(df)
head(df)

df %>% filter(n < 10)

df$marker1 <- factor(df$marker1, levels = unique(df$marker1))
df$marker2 <- factor(df$marker2, levels = unique(df$marker2))
head(df)



#======================#
#      TEST PLOT       #
#======================#

df_sub <- df %>% filter(SITE_ID == "ML") %>% filter(trait == "Height")
head(df_sub)


p1 <- ggplot(df_sub, aes(x = marker1, y = marker2, fill = r)) +
  geom_tile(color = "white") +
  geom_text(aes(label=round(r,2)), color="black", size=8) +
  scale_fill_distiller(palette = "Greys", limits = c(0.5, 1),
                       name = "Pearson's r", direction=-1) +
  coord_fixed() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size=12),
    axis.text.y = element_text(size=12),
    panel.grid = element_blank()
  ) +
  labs(x = "", y = "")
p1

################################################################################
traits <- c("Height", "Biomass_Increment", "DBH", "Average_Ring_Density", "Biomass_Increment_1980",
            "Biomass_Increment_1985","Biomass_Increment_1990","Biomass_Increment_1995")
trait.labels <- c("Height", "ABtot", "DBH","Density", "AB5-1980","AB5-1985","AB5-1990","AB5-1995")

df.sub <- df %>% filter(trait %in% traits)
df.sub$trait_labels <- factor(df.sub$trait, levels = traits, labels = trait.labels)
df.sub$SITE_ID <- factor(df.sub$SITE_ID, levels = c("PR","ML","CH","AC"))
head(df.sub)


p2 <- ggplot(df.sub, aes(x = marker1, y = marker2, fill = r)) +
  geom_tile(color = "white") +
  geom_text(aes(label=round(r,2)), color="black", size=3) +
  scale_fill_distiller(palette = "RdYlBu", limits = c(-1, 1),
                       name = "Pearson r", direction=-1) +
  #coord_fixed() +
  facet_grid(trait_labels~SITE_ID) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size=16),
    axis.text.y = element_text(size=16),
    panel.grid = element_blank(),
    strip.text = element_text(size=18),
    legend.title = element_text(size=16),
    legend.text = element_text(size=16)
  ) +
  labs(x = "", y = "")
p2



png("13_plot_offset_corrs_I.png", w = 3400, h = 4000, res=300)
p2
dev.off()


###############################

traits2 <- c("Biomass_Increment_2000",
             "Biomass_Increment_2005", "Biomass_Increment_2010", "Biomass_Increment_2015",
             "Rc", "Rr", "Rl", "Rs")
trait.labels2 <- c("AB5-2000","AB5-2005", "AB5-2010", "AB5-2015", "Rc", "Rr", "Rs", "Rt")

df.sub <- df %>% filter(trait %in% traits2)
df.sub$trait_labels <- factor(df.sub$trait, levels = traits2, labels = trait.labels2)
df.sub$SITE_ID <- factor(df.sub$SITE_ID, levels = c("PR","ML","CH","AC"))
head(df.sub)


p3 <- ggplot(df.sub, aes(x = marker1, y = marker2, fill = r)) +
  geom_tile(color = "white") +
  geom_text(aes(label=round(r,2)), color="black", size=3) +
  scale_fill_distiller(palette = "RdYlBu", limits = c(-1, 1),
                       name = "Pearson r", direction=-1) +
  #coord_fixed() +
  facet_grid(trait_labels~SITE_ID) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size=16),
    axis.text.y = element_text(size=16),
    panel.grid = element_blank(),
    strip.text = element_text(size=18),
    legend.title = element_text(size=16),
    legend.text = element_text(size=16)
  ) +
  labs(x = "", y = "")
p3



png("13_plot_offset_corrs_II.png", w = 3400, h = 4000, res=300)
p3
dev.off()



