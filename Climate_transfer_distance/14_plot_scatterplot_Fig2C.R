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

dclim <- read.csv("./13_plot_climate_distance_input.tsv", sep="\t", header=T)
dclim %>% head()
class(dclim$POP_ID)






#########

traits <- c("Height", "Biomass_Increment", "Survival")
dsub <- dclim %>% dplyr::filter(Trait %in% traits)
head(dsub)

dsub$label <- dsub$Trait
dsub$label <- factor(dsub$label, levels = c("Height", "Biomass_Increment", "Survival"), 
                     labels = c("Height\n[cm]","Log ABtot\n[kg]","Survival\n[%]"))
dsub$SITE_ID <- factor(dsub$SITE_ID, levels = c("PR","ML","CH", "AC"))
head(dsub)

p2 <- ggplot(dsub) + 
  aes(x = EuclDist, y = mean) +
  #geom_smooth(method = "lm", color="grey", fill = "lightgrey") +
  geom_errorbar(aes(ymin = mean - sd,
                    ymax = mean + sd),
                width = 0.0001) +
  geom_point(aes(fill=POP_GR), pch=21, size=2) +
  scale_fill_manual(values = brewer.pal(9, "BrBG"), name="Region") +
  labs(x = "Climate transfer distance", y = "") +
  
  #facet_wrap(Trait~SITE_ID, scales = "free", ncol=4) +
  facet_grid(label ~ SITE_ID, scales = "free_y") +
  theme(panel.background = element_rect(color="black", fill=NA),
        panel.grid = element_blank(),
        strip.text.x = element_text(size=14),
        strip.text.y = element_text(size=14),
        axis.title = element_text(size=18),
        legend.position = "right",
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_text(size=14),
        legend.text = element_text(size=14),
  ) +
  guides(fill = guide_legend(override.aes = list(size = 5)))
p2




png("14_plot_scatterplot.png", h=1200, w=3400, res=300)
p2
dev.off()
