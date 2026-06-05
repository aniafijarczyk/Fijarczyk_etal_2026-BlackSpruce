rm(list=ls())
setwd("C:/Users/aniaf/Projects/BlackSpruce/956_garden_offset_plot")

library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
library(RColorBrewer)


#=================#
#      DATA       #
#=================#



df1 <- read.csv("./02_garden_pops.tsv", sep="\t", header=T)
df1$Dataset <- "Populations"
df1$Model <- "GF"
df1 <- df1 %>% filter(n_samples >=5)
head(df1)

df2 <- read.csv("./03_garden_pops_rda.tsv", sep="\t", header=T)
df2$Dataset <- "Populations"
df2$Model <- "RDA"
df2 <- df2 %>% filter(n_samples >=5)
head(df2)

dm <- rbind(df1, df2)
dm <- dm %>% filter(SET == "TEST")
head(dm)


################################################################################


traits <- c("Height", "Biomass_Increment", "DBH", "Average_Ring_Density", "Biomass_Increment_1980",
            "Biomass_Increment_1985","Biomass_Increment_1990","Biomass_Increment_1995",
            "Biomass_Increment_2000",
            "Biomass_Increment_2005", "Biomass_Increment_2010", "Biomass_Increment_2015",
            "Rc", "Rr", "Rl", "Rs")
trait.labels <- c("Height", "ABtot", "DBH","Density", "AB5-1980","AB5-1985","AB5-1990","AB5-1995",
                  "AB5-2000","AB5-2005", "AB5-2010", "AB5-2015", "Rc", "Rr", "Rs", "Rt")


#traits <- c("Height")
#trait.labels <- c("Height")

df.sub <- dm %>% filter(trait %in% traits)

df.sub$trait_labels <- factor(df.sub$trait, levels = traits, labels = trait.labels)
unique(df.sub$n_samples)
df.sub$SITE_ID <- factor(df.sub$SITE_ID, levels = c("PR","ML","CH","AC"))
head(df.sub)


p1 <- ggplot(df.sub) + aes(x = n_samples, y = spearman_rho, color=SITE_ID, shape = Model) +
  geom_point(size=2, alpha=0.1) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2),
              se = TRUE, linewidth=1,
              aes(fill = SITE_ID, linetype = Model, alpha = Model)) +
  
  facet_wrap(.~trait_labels, ncol=8) +
  scale_color_manual(values = c("#0072b2","#009e73","#e69f00",'black'), name="Garden") +
  scale_fill_manual(values = c("#0072b2","#009e73","#e69f00",'black'), name="Garden") +
  scale_shape_manual(values = c(19,21)) +
  scale_linetype_manual(values = c("solid","dashed")) +
  scale_alpha_manual(values = c(0.5,0.15)) +
  scale_x_continuous(trans = "log2") +
  
  labs(x = "Number of populations in train set", y = "Spearman rho") + 
  
  theme(#axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=11),
    axis.text.x = element_text(size=11),
    axis.text.y = element_text(size=12),
    axis.title=element_text(size=16),
    strip.text = element_text(size=13),
    panel.background = element_rect(color="grey", fill=NA),
    panel.grid = element_blank())
p1



#png("04_plot_correlations_pops_all.png", h = 2900, w = 3200, res=300)
#p1
#dev.off()



################################################################################


df.means <- df.sub %>%
  group_by(n_samples, SITE_ID, trait_labels, Model) %>%
  summarize(mean_spearman = mean(spearman_rho, na.rm = TRUE), .groups = "drop")
df.means


p2 <- ggplot(df.sub) +
  
  geom_point(data=df.sub,
             aes(x = n_samples, y = spearman_rho, color=SITE_ID, shape = Model),
             size=1,
             alpha=0.2) +
  
  geom_line(data = df.means, 
            aes(x = n_samples, y = mean_spearman, color=SITE_ID, linetype=Model,
                group = interaction(SITE_ID, Model)),
            linewidth = 1) +
  
  facet_wrap(.~trait_labels, ncol=8) +
  scale_color_manual(values = c("#0072b2","#009e73","#e69f00",'black'), name="Garden") +
  scale_fill_manual(values = c("#0072b2","#009e73","#e69f00",'black'), name="Garden") +
  scale_shape_manual(values = c(19,21)) +
  scale_linetype_manual(values = c("solid","dashed")) +
  #scale_alpha_manual(values = c(0.5,0.15)) +
  scale_x_continuous(trans = "log2") +
  
  labs(x = "Maximum number of populations in the training dataset", y = "Spearman rho") + 
  
  theme(#axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=11),
    axis.text.x = element_text(size=11),
    axis.text.y = element_text(size=12),
    axis.title=element_text(size=16),
    strip.text = element_text(size=13),
    panel.background = element_rect(color="grey", fill=NA),
    panel.grid = element_blank())
p2

png("04_plot_correlations_pops_all_II.png", h = 2900, w = 3200, res=300)
p2
dev.off()




