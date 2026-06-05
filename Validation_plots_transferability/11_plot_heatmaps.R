rm(list=ls())
setwd("C:/Users/aniaf/Projects/BlackSpruce/996_garden_plot_clusters")

library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
library(RColorBrewer)


#=================#
#      DATA       #
#=================#

df1 <- read.csv("09_bootstrap.tsv", sep="\t", header=T)
dim(df1)
df1 <- df1 %>% filter(n_samples >=5) %>% drop_na()
head(df1)
dim(df1)

df2 <- read.csv("10_bootstrap_clusters.tsv", sep="\t", header=T)
dim(df2)
df2 <- df2 %>% filter(n_samples >=5) %>% drop_na()
dim(df2)
head(df2)

df3 <- read.csv("../951_garden_plot/02_bootstrap.tsv", sep="\t", header=T)
df3 <- df3 %>% filter(n_samples >=5) %>% drop_na() %>% filter(marker == "1000")
dim(df3)
head(df3)
df3$train_group <- "combined"
df3$group <- "combined"
df3 <- df3[c('SET', 'marker', 'trait', 'SITE_ID', 'group', 'train_group','rho', 'CI_low','CI_high',
             'n_samples', 'n_boot')]
head(df3)


df4 <- read.csv("../951_garden_plot/03_bootstrap_clusters.tsv", sep="\t", header=T)
df4 <- df4 %>% filter(n_samples >=5) %>% drop_na() %>% filter(marker == "1000")
head(df4)
df4$train_group <- "combined"
df4 <- df4[c('SET', 'marker', 'trait', 'SITE_ID', 'group', 'train_group','rho', 'CI_low','CI_high',
             'n_samples', 'n_boot')]


### Combining
df <- rbind(df1, df2, df3, df4)
head(df)

### Selecting test only
df <- df %>% filter(SET == "TEST")


head(df)
df$nsign <- ifelse((df$CI_low<0 & df$CI_high<0), -1, 0)
head(df)
df$sign <- ifelse((df$CI_low>0 & df$CI_high>0), 1, df$nsign)
head(df)

###
#df5 <- read.csv("../50_garden_offset_plot/12_get_ctd_clusters.tsv", sep="\t", header=T)
#df5 <- df5 %>% filter(!group %in% c("ME","WI")) %>% filter(n_samples >=5)
#head(df5)


#======================#
#      TEST PLOT       #
#======================#



traits <- c("Height", "Biomass_Increment", "DBH", "Average_Ring_Density", "Biomass_Increment_1980","Biomass_Increment_1985","Biomass_Increment_1990","Biomass_Increment_1995","Biomass_Increment_2000",
            "Biomass_Increment_2005", "Biomass_Increment_2010", "Biomass_Increment_2015",
            "Rc", "Rr", "Rl", "Rs")

trait.labels <- c("Height", "ABtot", "DBH","Density", "AB5-1980","AB5-1985","AB5-1990","AB5-1995",
                  "AB5-2000","AB5-2005", "AB5-2010", "AB5-2015", "Rc", "Rr", "Rs", "Rt")

df.sub <- df
head(df.sub)
df.sub$trait_labels <- factor(df.sub$trait, levels = traits, labels = trait.labels)
df.sub$group <- factor(df.sub$group, levels = c("West","Central","East","combined"))
df.sub$train_group <- factor(df.sub$train_group, levels = c("West","Central","East","combined"))
df.sub$SITE_ID <- factor(df.sub$SITE_ID, levels = c("PR","ML","CH","AC"))


head(df.sub)
df.sub[df.sub$sign != 0,]

p1 <- ggplot(df.sub) +

  geom_tile(data = df.sub, aes(x = train_group, y = group, fill = rho)) +
  geom_point(data = df.sub[df.sub$sign != 0,], aes(x = train_group, y = group), color="black") +
  
  facet_grid(SITE_ID~trait_labels) +
  #facet_grid(trait_labels~SITE_ID) +
  scale_fill_distiller(palette = "RdYlBu", name="Spearman\nrho", limits = c(-1, 1)) +

  labs(x = "Cluster trained on", y = "Cluster tested") + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        panel.background = element_rect(color="grey", fill=NA),
        panel.grid = element_blank(),
        legend.position = "right")
p1


#png("09_plot_correlations_clusters_ALL.png", h = 1500, w = 5000, res=300)
#p1
#dev.off()

####################################################### PART I

traits <- c("Height", "Biomass_Increment", "DBH", "Average_Ring_Density", "Biomass_Increment_1980",
            "Biomass_Increment_1985","Biomass_Increment_1990","Biomass_Increment_1995")
trait.labels <- c("Height", "ABtot", "DBH","Density", "AB5-1980","AB5-1985","AB5-1990","AB5-1995")

df.sub <- df %>% filter(trait %in% traits)
df.sub$trait_labels <- factor(df.sub$trait, levels = traits, labels = trait.labels)
df.sub$group <- factor(df.sub$group, levels = c("West","Central","East","combined"))
df.sub$train_group <- factor(df.sub$train_group, levels = c("West","Central","East","combined"))
df.sub$SITE_ID <- factor(df.sub$SITE_ID, levels = c("PR","ML","CH","AC"))



p1a <- ggplot(df.sub) +
  
  geom_tile(data = df.sub, aes(x = train_group, y = group, fill = rho)) +
  geom_point(data = df.sub[df.sub$sign != 0,], aes(x = train_group, y = group), color="black") +
  
  facet_grid(SITE_ID~trait_labels) +
  #facet_grid(trait_labels~SITE_ID) +
  scale_fill_distiller(palette = "RdYlBu", name="Spearman\nrho", limits = c(-1, 1)) +
  
  labs(x = "Cluster trained on", y = "Cluster tested") + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        panel.background = element_rect(color="grey", fill=NA),
        panel.grid = element_blank(),
        legend.position = "right")
p1a



##################################################################################

traits2 <- c("Biomass_Increment_2000",
             "Biomass_Increment_2005", "Biomass_Increment_2010", "Biomass_Increment_2015",
             "Rc", "Rr", "Rl", "Rs")
trait.labels2 <- c("AB5-2000","AB5-2005", "AB5-2010", "AB5-2015", "Rc", "Rr", "Rs", "Rt")

df.sub <- df %>% filter(trait %in% traits2)
head(df.sub)
df.sub$trait_labels <- factor(df.sub$trait, levels = traits2, labels = trait.labels2)
df.sub$group <- factor(df.sub$group, levels = c("West","Central","East","combined"))
df.sub$train_group <- factor(df.sub$train_group, levels = c("West","Central","East","combined"))
df.sub$SITE_ID <- factor(df.sub$SITE_ID, levels = c("PR","ML","CH","AC"))

p1b <- ggplot(df.sub) +
  
  geom_tile(data = df.sub, aes(x = train_group, y = group, fill = rho)) +
  geom_point(data = df.sub[df.sub$sign != 0,], aes(x = train_group, y = group), color="black") +
  
  facet_grid(SITE_ID~trait_labels) +
  scale_fill_distiller(palette = "RdYlBu", name="Spearman\nrho", limits = c(-1, 1)) +
  
  labs(x = "Cluster trained on", y = "Cluster tested") + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        panel.background = element_rect(color="grey", fill=NA),
        panel.grid = element_blank(),
        legend.position = "right")
p1b


png("11_plot_heatmaps.png", h = 3000, w = 3000, res=300)
plot_grid(p1a, p1b, ncol = 1)
dev.off()



#======================#
#     PRETTY PLOT      #
#======================#


traits <- c("Height", "Biomass_Increment")
trait.labels <- c("Height", "ABtot")


df.sub <- df %>% filter(trait %in% traits)
head(df.sub)
df.sub$trait_labels <- factor(df.sub$trait, levels = traits, labels = trait.labels)
df.sub$group <- factor(df.sub$group, levels = c("West","Central","East","combined"))
df.sub$train_group <- factor(df.sub$train_group, levels = c("West","Central","East","combined"))
df.sub$SITE_ID <- factor(df.sub$SITE_ID, levels = c("PR","ML","CH","AC"))

head(df.sub)

p2 <- ggplot(df.sub) +
  
  geom_tile(data = df.sub, aes(x = train_group, y = group, fill = rho)) +
  geom_text(aes(x = train_group, y = group, label = round(rho,2)), color = "black", size = 3) + # Add labels
  
  # bold labels for sign != 0
  geom_text(
    data = subset(df.sub, sign != 0),
    aes(x = train_group, y = group, label = round(rho, 2)),
    fontface = "bold", color = "black", size = 3
  ) +
  
  facet_grid(trait_labels~SITE_ID) +
  scale_fill_distiller(palette = "RdYlBu", name="Spearman\nrho", limits = c(-1, 1)) +
  
  labs(x = "Cluster trained on", y = "Cluster tested") + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        panel.background = element_rect(color="grey", fill=NA),
        panel.grid = element_blank(),
        legend.position = "right")
p2



png("11_plot_heatmaps_selected.png", w = 2200, h = 1000, res=300)
p2
dev.off()


###################################
