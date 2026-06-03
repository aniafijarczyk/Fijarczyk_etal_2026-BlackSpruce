rm(list=ls())
setwd("C:/Users/aniaf/Projects/BlackSpruce/929_heterozygosity")

library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
library(RColorBrewer)

library(sf)
library(raster)
library(spData)
library(tmap)




#===================#
#       DATA        #
#===================#


df.all <- read.csv("02_calculate_heterozygosity.tsv", sep="\t", header=T)
head(df.all)

# Mean over iterations
df.all %>% head()
df <- df.all %>% group_by(POP) %>% summarise(Ho_mean = mean(Ho), Ho_sd = sd(Ho),
                                             Ho_err = sd(Ho) / sqrt(sum(!is.na(Ho))),
                                             He_mean = mean(He),
                                             He_sd = sd(He),
                                             He_err = sd(He) / sqrt(sum(!is.na(He))),
                                             N_poly_mean = mean(N_poly),
                                             N_poly_sd = sd(N_poly),
                                             N_poly_err = sd(N_poly) / sqrt(sum(!is.na(N_poly))))

dim(df)
df
samples <- df$POP
length(samples)



#===================#
#     Metadata      #
#===================#

# Adding geographic coordinates

dmeta <- read.csv("../DATA_intermediate/23_filter_EPN_indiv_metrics.tsv", sep="\t")
head(dmeta)
dmeta <- dmeta %>% dplyr::select(POP, STRATA, lat, lon) %>% distinct()
head(dmeta)
dim(dmeta)

# adding admixture level and call rate

dadm <- read.csv("../03_admixture_EPR_EPN/admixture_04_combine_output_assigned_clusters.tsv", sep="\t")
dadm <- dadm %>% dplyr::select(id, POP, BestK6, BestK6_Q, BATCH_ID, call_rate)
head(dadm)
dadm$BestK6_Q
dadm$class <- ifelse(dadm$BestK6_Q > 0.75, 0, 1)
head(dadm)
pop_adm <- dadm %>% group_by(POP) %>% dplyr::summarise(prop_mixed = sum(class)/n(),
                                            mean_Q = mean(BestK6_Q),
                                            mean_call_rate = mean(call_rate))
pop_adm
  

# Merging table

dm0 <- merge(df, dmeta, by="POP", sort=F)
dim(dm0)
head(dm0)
dm <- merge(dm0, pop_adm, by = "POP", sort=F)
head(dm)
write.table(dm, "03_plot_heterozygosity.tsv", sep="\t", col.names = T, row.names = F, quote=F, append=F)



#===================#
#      Plots        #
#===================#

# Polymorphic sites vs Ho, colored by strata

head(dm)

p2 <- ggplot(dm) +
  aes(x = N_poly_mean, y = Ho_mean) +
  geom_errorbar(aes(xmin = N_poly_mean - N_poly_sd,
                    xmax = N_poly_mean + N_poly_sd)) +
  geom_errorbar(aes(ymin = Ho_mean - Ho_sd,
                    ymax = Ho_mean + Ho_sd)) +
  geom_point(pch = 21, aes(fill = STRATA), size=4) +
  scale_fill_manual(values = c(brewer.pal(9, "BrBG"),"firebrick1"), name="Region") +
  labs(y = "Observed heterozygosity", x = "Polymorphic sites") +
  theme(panel.background = element_rect(fill=NA, color="grey20"),
        panel.grid = element_blank(),
        axis.title = element_text(size=20),
        legend.position = "none",
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20),
  )
p2


png("03_plot_heterozygosity_Ho.png", w = 1600, h = 1200, res=300)
p2
dev.off()




# Polymorphic sites vs Ho, colored by admixture level

p3 <- ggplot(dm) +
  aes(x = N_poly_mean, y = Ho_mean) +
  geom_errorbar(aes(xmin = N_poly_mean - N_poly_sd,
                    xmax = N_poly_mean + N_poly_sd)) +
  geom_errorbar(aes(ymin = Ho_mean - Ho_sd,
                    ymax = Ho_mean + Ho_sd)) +
  geom_point(pch = 21, aes(fill = prop_mixed), size=4) +
  scale_fill_continuous(name="Proportion\nof\nadmixed\ntrees") +
  labs(y = "Observed heterozygosity", x = "Polymorphic sites") +
  theme(panel.background = element_rect(fill=NA, color="grey20"),
        panel.grid = element_blank(),
        axis.title = element_text(size=20),
        #legend.position = "none",
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_text(size=16),
        legend.text = element_text(size=20),
  )
p3


# Polymorphic sites vs Ho, colored by call rate

p5 <- ggplot(dm) +
  aes(x = N_poly_mean, y = Ho_mean) +
  geom_errorbar(aes(xmin = N_poly_mean - N_poly_sd,
                    xmax = N_poly_mean + N_poly_sd)) +
  geom_errorbar(aes(ymin = Ho_mean - Ho_sd,
                    ymax = Ho_mean + Ho_sd)) +
  geom_point(pch = 21, aes(fill = mean_call_rate), size=4) +
  scale_fill_continuous(name="Mean_call_rate") +
  labs(y = "Observed heterozygosity", x = "Polymorphic sites") +
  theme(panel.background = element_rect(fill=NA, color="grey20"),
        panel.grid = element_blank(),
        axis.title = element_text(size=20),
        #legend.position = "none",
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_text(size=16),
        legend.text = element_text(size=20),
  )
p5







# Polymorphic sites vs He, colored by strata

p4 <- ggplot(dm) +
  aes(x = N_poly_mean, y = He_mean) +
  geom_errorbar(aes(xmin = N_poly_mean - N_poly_sd,
                    xmax = N_poly_mean + N_poly_sd)) +
  geom_errorbar(aes(ymin = He_mean - He_sd,
                    ymax = He_mean + He_sd)) +
  geom_point(pch = 21, aes(fill = STRATA), size=4) +
  scale_fill_manual(values = c(brewer.pal(9, "BrBG"),"firebrick1"), name="Region") +
  labs(y = "Expected heterozygosity", x = "Polymorphic sites") +
  theme(panel.background = element_rect(fill=NA, color="grey20"),
        panel.grid = element_blank(),
        axis.title = element_text(size=20),
        legend.position = "none",
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20),
  )
p4


png("03_plot_heterozygosity_He.png", w = 1600, h = 1200, res=300)
p4
dev.off()
