rm(list=ls())
setwd("C:/Users/aniaf/Projects/BlackSpruce/947_GEA_RDA")

library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
library(RColorBrewer)

library(sf)
library(raster)
library(spData)
library(tmap)




#------------#
#    DATA    #
#------------#

df <- read.csv("03_heterozygosity_calculate_MAF.tsv", sep="\t", header=T)
head(df)
dim(df)

samples <- df$POP_SITE
length(samples)


# Metadata
dmeta <- read.csv("../DATA_intermediate/23_filter_EPN_indiv_metrics.tsv", sep="\t")
dmeta <- dmeta %>% dplyr::select(POP, STRATA, SITE_ID, lat, lon) %>% distinct()
dmeta$POP_SITE <- paste0(dmeta$POP, "_", dmeta$SITE_ID)
head(dmeta)
dim(dmeta)


dm <- merge(df, dmeta, by="POP_SITE", sort=F)
dim(dm)
head(dm)

write.table(dm, "04_plot_heterozygosity.tsv", sep="\t", col.names = T, row.names = F, quote=F, append=F)



#########################################################
#-                       Plots                         -#
#########################################################


#####################
head(dm)

dh <- dm %>% dplyr::select(POP_SITE, Hobs_all, Hobs_RDA, Hobs_RDAcor, Hobs_LFMM, Trait) %>% gather(key = "Set", value = "Hobs", -POP_SITE, -Trait)
head(dh)

traits <- c("Height", "Biomass_Increment", "DBH", "Average_Ring_Density", "Biomass_Increment_1980",
            "Biomass_Increment_1985","Biomass_Increment_1990","Biomass_Increment_1995",
            "Biomass_Increment_2000","Biomass_Increment_2005", "Biomass_Increment_2010", "Biomass_Increment_2015",
            "Rc", "Rr", "Rl", "Rs")
trait.labels <- c("Height", "ABtot", "DBH","Density", "AB5-1980","AB5-1985","AB5-1990","AB5-1995",
                  "AB5-2000","AB5-2005", "AB5-2010", "AB5-2015", "Rc", "Rr", "Rs", "Rt")


dh$trait_labels <- factor(dh$Trait, levels = traits, labels = trait.labels)
dh$Set <- factor(dh$Set, levels = c("Hobs_all","Hobs_LFMM","Hobs_RDA","Hobs_RDAcor"),
                 labels = c("All","LFMM","RDA","RDA-struct"))
head(dh)

p2 <- ggplot(dh) +
  aes(x = Set, y = Hobs) +
  geom_jitter() +
  geom_boxplot(outlier.shape = NA, alpha=0.5) +
  #scale_fill_manual(values = c(brewer.pal(9, "BrBG"),"firebrick1"), name="Region") +
  labs(y = "Observed heterozygosity", x = "Dataset") +
  facet_wrap(.~trait_labels) +
  theme(panel.background = element_rect(fill=NA, color="grey20"),
        panel.grid = element_blank(),
        axis.title = element_text(size=20),
        axis.text.x = element_text(angle=90, size=14, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=14),
        strip.text = element_text(size=14),
        legend.position = "none",
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20),
  )
p2



png("04_plot_heterozygosity_Hobs.png", w = 2600, h = 2200, res=300)
p2
dev.off()


##################################

dp <- dm %>% dplyr::select(POP_SITE, Poly_all, Poly_RDA, Poly_RDAcor, Poly_LFMM, Trait) %>% gather(key = "Set", value = "Poly", -POP_SITE, -Trait)
head(dp)

p3 <- ggplot(dp) +
  aes(x = Set, y = Poly) +
  geom_jitter() +
  geom_boxplot(outlier.shape = NA, alpha=0.5) +
  #scale_fill_manual(values = c(brewer.pal(9, "BrBG"),"firebrick1"), name="Region") +
  labs(y = "Polymorphic sites", x = "Dataset") +
  facet_wrap(.~Trait) +
  theme(panel.background = element_rect(fill=NA, color="grey20"),
        panel.grid = element_blank(),
        axis.title = element_text(size=20),
        axis.text.x = element_text(angle=90),
        legend.position = "none",
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20),
  )
p3

png("04_plot_heterozygosity_Poly.png", w = 2600, h = 2200, res=300)
p3
dev.off()









