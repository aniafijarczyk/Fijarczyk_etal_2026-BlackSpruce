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


##################################################################################

df <- df12
df$marker %>% unique()

traits <- c("Height", "Biomass_Increment", "DBH", "Average_Ring_Density", "Biomass_Increment_1980",
            "Biomass_Increment_1985","Biomass_Increment_1990","Biomass_Increment_1995",
            "Biomass_Increment_2000",
            "Biomass_Increment_2005", "Biomass_Increment_2010", "Biomass_Increment_2015",
            "Rc", "Rr", "Rl", "Rs")
trait.labels <- c("Height", "Biomass inc.", "DBH","Density", "BI-1980","BI-1985","BI-1990","BI-1995",
                  "BI-2000",
                  "BI-2005", "BI-2010", "BI-2015", "Rc", "Rr", "Rs", "Rt")
ds <- df %>% filter(trait %in% traits) %>% filter(marker %in% c("100","100LF","1000","1000LF","all","lfmm","RDA","RDAcorrected"))
head(ds)

ds$trait_labels <- factor(ds$trait, levels = traits, labels = trait.labels)
ds$SET <- factor(ds$SET, levels = c("TRAIN","TEST"))
ds$SITE_ID <- factor(ds$SITE_ID, levels = c("PR","ML","CH","AC"))
ds$marker_labels <- factor(ds$marker, levels = c("100","100LF","1000","1000LF","all","lfmm","RDA","RDAcorrected"),
                           labels = c("0.1 K", "0.1 K-lf", "1 K", "1 K-lf", "~29 K", "LFMM", "RDA", "RDA-struct"))


#ds <- ds %>% filter(trait %in% c("Height"))

ds_wide <- ds %>% dplyr::select(-n_samples, -n_boot) %>%
  pivot_wider(
    names_from = SET,
    values_from = c(rho, CI_low, CI_high),
    names_glue = "{SET}_{.value}"
  )

ds_wide$sign_neg <- ifelse(((ds_wide$TRAIN_CI_low<0) & (ds_wide$TRAIN_CI_high<0) & (ds_wide$TEST_CI_low<0) & (ds_wide$TEST_CI_high<0)),1,0)
ds_wide$sign_pos <- ifelse(((ds_wide$TRAIN_CI_low>0) & (ds_wide$TRAIN_CI_high>0) & (ds_wide$TEST_CI_low>0) & (ds_wide$TEST_CI_high>0)),1,0)
ds_wide$sign <- rowSums(ds_wide[, c("sign_neg", "sign_pos")], na.rm = TRUE)
as.data.frame(ds_wide) %>% head()


p1 <- ggplot(ds_wide) +
  aes(x = TRAIN_rho, y = TEST_rho) +
  geom_hline(yintercept = 0, linetype="solid", color = "grey80", linewidth=1) +
  geom_vline(xintercept = 0, linetype="solid", color = "grey80", linewidth=1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +

  #geom_errorbar(data=ds_wide, aes(x = TRAIN_rho, xmin = TRAIN_CI_low, xmax = TRAIN_CI_high,
  #                                group = marker_labels, color = SITE_ID), width=0.001) +
  #geom_errorbar(data=ds_wide, aes(y = TEST_rho, ymin = TEST_CI_low, ymax = TEST_CI_high,
  #                                group = marker_labels, color = SITE_ID), width=0.001) +
  
  geom_point(data=ds_wide, aes(shape = marker_labels, color=SITE_ID, alpha = as.factor(sign))) +
    
  scale_shape_manual(values = c(1, 2, 3, 4, 8, 15, 17, 23), name = "Marker") +
  scale_color_manual(values = c("#0072b2","#009e73","#e69f00",'black'), name="Garden") +
  scale_alpha_manual(values = c(0.25, 1), labels = c("No","Yes"), name="Different\nfrom 0") +
  facet_wrap(trait_labels~test, ncol = 4) +
  #ggtitle("Gradient Forest") +
  theme(panel.background = element_rect(fill=NA, color="grey"),
        panel.grid = element_blank(),
        axis.text = element_text(size=10),
        axis.title = element_text(size=14),
        strip.text = element_text(size=10),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = "right")
p1


png("10_plot_correlations_train_all.png", h = 3600, w = 3400, res=300)
p1
dev.off()



##################################################################################
