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


##################
# Saving dataset

clim.phenos <- read.csv("13_plot_climate_distance_input.tsv", sep="\t", header=T)
head(clim.phenos)
dim(clim.phenos)

#####################
# Filtering for number of ind per group ???
#flim.phenos <- clim.phenos %>% filter(n>4)
#dim(flim.phenos)





############################################################################################################
######################### adding cluster info


### Metadata with groups (for 66 pops)
dmeta <- read.csv("../DATA_intermediate/23_filter_EPN_indiv_metrics.tsv", sep="\t", header=T)
dmeta$POP_ID <- dmeta$POP
# All filtered provenances (n=66) 
dmeta <- dmeta %>% dplyr::select(POP_ID, group) %>% distinct()
dim(dmeta)
head(dmeta)

dsub.gr <- merge(clim.phenos, dmeta, by = "POP_ID", all.x = T)
dsub.gr %>% filter(!group %in% c("West","Central","East","WI","ME"))
dsub.gr <- dsub.gr %>% filter(group %in% c("West","Central","East"))
head(dsub.gr)

small_sample_sizes.a <- dsub.gr %>% filter(SITE_ID == "ML" & group == "West") %>% pull(POP_SITE) %>% unique()
small_sample_sizes.b <- dsub.gr %>% filter(SITE_ID == "PR" & group == "East") %>% pull(POP_SITE) %>% unique()
small_sample_sizes <- c(small_sample_sizes.b, small_sample_sizes.a)
small_sample_sizes

dsub.gr <- dsub.gr %>% filter(!POP_SITE %in% small_sample_sizes)
head(dsub.gr)




############################################################################################################


traits <- c("Height", "Biomass_Increment", "DBH", "Survival", "Average_Ring_Density", "Biomass_Increment_1980",
            "Biomass_Increment_1985","Biomass_Increment_1990","Biomass_Increment_1995")
trait.labels <- c("Height", "Biomass inc.", "DBH", "Survival","Density", "BI-1980","BI-1985","BI-1990","BI-1995")

dsub <- dsub.gr %>% dplyr::filter(Trait %in% traits)
dsub$group <- factor(dsub$group, levels = c("West","Central","East"))
dsub$SITE_ID <- factor(dsub$SITE_ID, levels = c("PR","ML","CH","AC"))
dsub$label <- factor(dsub$Trait, levels = traits, labels = trait.labels)
dim(dsub)
head(dsub)


p3a <- ggplot(dsub) + 
  aes(x = EuclDist, y = mean) +
  #geom_smooth(method = "lm", color="grey", fill = "lightgrey") +
  geom_errorbar(aes(ymin = mean - sd,
                    ymax = mean + sd),
                width = 0.0001) +
  geom_point(aes(fill=group), pch=21, size=2) +
  scale_fill_manual(values = c("#DFC27D","#C7EAE5", "#35978F"), name="Cluster") +
  labs(x = "Climate transfer distance", y = "") +
  geom_smooth(aes(color=group), method='lm', se = FALSE) +
  scale_color_manual(values = c("#DFC27D","#C7EAE5", "#35978F"), name="Cluster") +
  #facet_wrap(Trait~SITE_ID, scales = "free", ncol=4) +
  facet_grid(label~SITE_ID, scales = "free_y") +
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
p3a

##################

traits2 <- c("Biomass_Increment_2000",
             "Biomass_Increment_2005", "Biomass_Increment_2010", "Biomass_Increment_2015",
             "Rc", "Rr", "Rl", "Rs")
#trait.labels2 <- c("BI-2000",
#                   "BI-2005", "BI-2010", "BI-2015", "Rc", "Rr", "Rl", "Rs")
trait.labels2 <- c("BI-2000",
                   "BI-2005", "BI-2010", "BI-2015", "Rc", "Rr", "Rs", "Rt")

dsub <- dsub.gr %>% dplyr::filter(Trait %in% traits2)
dsub$group <- factor(dsub$group, levels = c("West","Central","East"))
dsub$SITE_ID <- factor(dsub$SITE_ID, levels = c("PR","ML","CH","AC"))
dsub$label <- factor(dsub$Trait, levels = traits2, labels = trait.labels2)
dim(dsub)
head(dsub)

p3b <- ggplot(dsub) + 
  aes(x = EuclDist, y = mean) +
  #geom_smooth(method = "lm", color="grey", fill = "lightgrey") +
  geom_errorbar(aes(ymin = mean - sd,
                    ymax = mean + sd),
                width = 0.0001) +
  geom_point(aes(fill=group), pch=21, size=2) +
  scale_fill_manual(values = c("#DFC27D","#C7EAE5", "#35978F"), name="Cluster") +
  labs(x = "Climate transfer distance", y = "") +
  geom_smooth(aes(color=group), method='lm', se = FALSE) +
  scale_color_manual(values = c("#DFC27D","#C7EAE5", "#35978F"), name="Cluster") +
  #facet_wrap(Trait~SITE_ID, scales = "free", ncol=4) +
  facet_grid(label~SITE_ID, scales = "free_y") +
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
p3b



png("17_plot_climate_distance_scatterplots_clustersA.png", h=3500, w=2400, res=300)
p3a
dev.off()

png("17_plot_climate_distance_scatterplots_clustersB.png", h=3500, w=2400, res=300)
p3b
dev.off()
