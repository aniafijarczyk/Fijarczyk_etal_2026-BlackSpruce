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

df <- read.csv("15_compute_correlations.tsv", sep="\t", header=T)
head(df)

### only correlations
dg <- df %>% dplyr::select(Trait, SITE_ID, rho, p.value, n_trees, adj.p) %>% distinct()
dim(dg)
head(dg)

### small sample sizes
dg %>% filter(n_trees < 10)

#============#
#   PLOTS    #
#============#


unique(dg$Trait)

traits <- c("Height", "Biomass_Increment", "Survival", "DBH", "Average_Ring_Density", "Biomass_Increment_1980","Biomass_Increment_1985","Biomass_Increment_1990","Biomass_Increment_1995","Biomass_Increment_2000",
            "Biomass_Increment_2005", "Biomass_Increment_2010", "Biomass_Increment_2015",
            "Rc", "Rr", "Rl", "Rs")

#trait.labels <- c("Height", "Biomass inc.", "Survival", "DBH","Density", "BI-1980","BI-1985","BI-1990","BI-1995","BI-2000",
#            "BI-2005", "BI-2010", "BI-2015", "Rc", "Rr", "Rl", "Rs")
trait.labels <- c("Height", "ABtot", "Survival", "DBH","Density", "AB5-1980","AB5-1985","AB5-1990","AB5-1995","AB5-2000",
                  "AB5-2005", "AB5-2010", "AB5-2015", "Rc", "Rr", "Rs", "Rt")


dg$Trait <- factor(dg$Trait, levels = traits, labels = trait.labels)
dg$SITE_ID <- factor(dg$SITE_ID, levels = c("PR","ML","CH", "AC"))

p1 <- ggplot(dg) + aes(y = SITE_ID, x = Trait) +
  geom_tile(data=dg, aes(fill = rho), color="white") +
  geom_text(data=dg, aes(label = round(rho, 2)), color = "black") +
  geom_text(data=dg[dg$adj.p<0.05,], aes(label = round(rho, 2)), fontface="bold", color = "black") +
  scale_fill_distiller(palette = "Spectral", direction = -1, name = "Spearman\nrho") +
  coord_fixed() +
  labs(y = "Common garden", x = "Trait") +
  theme(panel.background = element_blank(),
        axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5, size=14),
        axis.text.y = element_text(size=14),
        axis.title = element_text(size=18),
        legend.title = element_text(size=14),
        axis.ticks = element_blank())
p1



png("16_plot_heatmaps.png", h=1000, w=3400, res=300)
p1
dev.off()
