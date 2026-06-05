rm(list=ls())
setwd("C:/Users/aniaf/Projects/BlackSpruce/980_provenances_offset_plot")

library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape)
library(cowplot)
library(grid)
library(RColorBrewer)

library(sf)
library(raster)
library(spData)
library(tmap)

#===========#
#   DATA    #
#===========#


dpheno1 <- read.csv("01_garden_compare_PCs.tsv", sep="\t", header=T)
dpheno2 <- read.csv("01_garden_compare_var3.tsv", sep="\t", header=T)


unique(dpheno1$year)
unique(dpheno1$scenario)
unique(dpheno1$cluster)
unique(dpheno1$set)
class(dpheno1$year)

dmeta <- read.csv("../METADATA/populations_genetic_groups.tsv", sep="\t", header=T)
head(dmeta)
dmeta$POP_ID
dim(dmeta)
tail(dmeta)

dpheno1 <- dpheno1 %>% dplyr::select(-group)
dpheno1 <- merge(dpheno1, dmeta, by.x = "KeyID", by.y = "POP_ID", all.x=T, sort=F)
head(dpheno1)

dpheno2 <- dpheno2 %>% dplyr::select(-group)
dpheno2 <- merge(dpheno2, dmeta, by.x = "KeyID", by.y = "POP_ID", all.x=T, sort=F)


### Subset populations in


#=====================#
#   GENETIC OFFSET    #
#=====================#






####################################################################################

#=====================#
#        MAPS         #
#=====================#


world_NA = world[world$continent == "North America", ]
world_NA$geom
head(world_NA)
my_sf <- read_sf("C:/Users/aniaf/Projects/BlackSpruce/METADATA/maps/BlackSpruce/data/commondata/data0/picemari.shp")
species_distribution_outline <- st_union(my_sf)
ne_sf <- st_transform(species_distribution_outline, crs = st_crs(world))
head(ne_sf)

head(dpheno1)
dpheno1$cluster %>% unique()
dpheno1$set %>% unique()
dpheno1$scenario %>% unique()
dpheno1$year %>% unique()
dpheno1$KeyID %>% unique()
dpheno1$group %>% unique()

dm <- dpheno1 %>% filter(cluster == "combined") %>% filter(set == "PC") %>% filter(scenario == "245") %>% 
  filter(dataset == "Provenances")
dm$group %>% unique()
head(dm)

dm.CE <- dm %>% filter(group %in% c("Central","East")) %>% filter(Longitude > -125)
#dm.CE <- dm %>% filter(group != "West") %>% filter(Longitude > -125)

p2 <- ggplot(dm.CE) + aes(x = Longitude, y = Latitude, fill = offset) +
  geom_point(pch=21, size=5) +
  facet_grid(scenario~year) +
  scale_fill_viridis_b(option="magma", direction = -1)
p2


# scale_color_viridis_b(option="magma", direction = -1, name="Genetic offset",
#                       breaks = seq(min(dataset$offset), max(dataset$offset), by = 0.005),
#                       labels = function(x) {
#                         labs <- sprintf("%.2f", x)
#                         labs[seq_along(labs) %% 2 == 0] <- ""  # blank every 2nd label
#                         labs})


################################################################################


plot_map2 <- function(dataset, min_long=-170, max_long=-52) {
  
  dataset$year <- factor(dataset$year, levels = c(2031, 2050, 2070, 2090),
                         labels = c("2022-40", "2041-60", "2061-80", "2081-2100"))
  dataset$scenario <- factor(dataset$scenario, levels = c(245, 370, 585),
                             labels = c("SSP2-4.5", "SSP3-7.0", "SSP5-8.5"))
  
  m1 <- ggplot() + 
    geom_sf(data = ne_sf, fill = "white") +
    
    geom_point(data = dataset, aes(x = Longitude, y = Latitude, color = offset), pch=19, size=2) +
    scale_color_viridis_b(option="magma", direction = -1, name="Genomic offset",
                          breaks = seq(min(dataset$offset), max(dataset$offset), by = 0.005),
                          labels = function(x) {
                            labs <- sprintf("%.2f", x)
                            labs[seq_along(labs) %% 4 != 1] <- ""  # blank every 5nd label
                            labs
                          }
    ) +
    #scale_x_continuous(limits = c(-170, -52)) +
    scale_x_continuous(limits = c(min_long, max_long)) +
    scale_y_continuous(limits = c(40, 70)) +
    labs(x = "Longitude", y = "Latitude") +
    facet_wrap(year~., ncol=1) +
    theme(
      panel.background = element_rect(fill="white", color=NA),
      panel.grid = element_blank(),
      #panel.background = element_rect(fill="aliceblue", color="grey20"),
      legend.position = "bottom",
      legend.background = element_blank(),
      legend.key = element_blank(),
      #legend.title = element_blank(),
      strip.background = element_rect(fill=NA, color=NA),
      strip.text.x = element_text(size=12),
      #strip.text.y.left = element_text(size=16, angle = 0),  # rotate facet labe
      strip.text.y = element_blank(),
      legend.text = element_text(size=10),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    ) +
    guides(
      color = guide_colorbar(
        barwidth = 10,    # increase width (in legend units)
        barheight = 1,     # adjust height if needed
        title.position = "top",   # move title above the bar
        title.hjust = 0.5
      ))
  print(m1)
}

plot_map2(dm.CE, min_long=-110, max_long=-52)

####################################

#dm1 <- dpheno1 %>% filter(cluster == "combined")
#dm2 <- dpheno1 %>% filter(cluster == "West")
#dm3 <- dpheno1 %>% filter(cluster == "Central")
#dm4 <- dpheno1 %>% filter(cluster == "East")

#m1 <- plot_map(dm1)
#m2 <- plot_map(dm2)
#m3 <- plot_map(dm3)
#m4 <- plot_map(dm4)
#plot_grid(m1, m2, m3, m4, ncol=1)


#### Model PC3 - Central/East

dm <- dpheno1 %>% filter(cluster == "combined") %>% filter(set == "PC") %>% filter(scenario == "245") %>% 
  filter(dataset == "Provenances")
dm$group %>% unique()
head(dm)
dm.CE <- dm %>% filter(group %in% c("Central","East")) %>% filter(Longitude > -125)
max(dm.CE$offset)
min(dm.CE$offset)

plot.ce <- plot_map2(dm.CE, min_long=-110, max_long=-52)



#### Model var3 - West
dm.2 <- dpheno2 %>% filter(cluster == "West") %>% filter(set == "var3") %>% filter(scenario == "245") %>% 
  filter(dataset == "Provenances")
dm.2$group %>% unique()
head(dm.2)
dm.W <- dm.2 %>% filter(group %in% c("West"))
max(dm.W$offset)
min(dm.W$offset)



plot.ce <- plot_map2(dm.CE, min_long=-110, max_long=-52)
plot.w <- plot_map2(dm.W, max_long=-90)



png("03_plot_offset_separate.png", res=300, w=2500,h=2500)
plot_grid(plot.w, plot.ce, ncol=2, labels = c("A","B"))
dev.off()


#png("03_plot_offset_separate.png", res=300, w=3400,h=2000)
#plot_grid(plot.w, plot.ce, ncol=1)
#dev.off()



#

# combine results from two models


##################################################################################


plot_map3 <- function(dataset, min_long=-170, max_long=-52) {
  
  dataset$year <- factor(dataset$year, levels = c(2031, 2050, 2070, 2090),
                         labels = c("2022-40", "2041-60", "2061-80", "2081-2100"))
  dataset$scenario <- factor(dataset$scenario, levels = c(245, 370, 585),
                             labels = c("SSP2-4.5", "SSP3-7.0", "SSP5-8.5"))
  
  m1 <- ggplot() + 
    geom_sf(data = ne_sf, fill = "white") +
    
    geom_point(data = dataset, aes(x = Longitude, y = Latitude, color = offset), pch=19, size=2) +
    scale_color_viridis_b(option="magma", direction = -1, name="Genetic offset",
                          breaks = seq(min(dataset$offset), max(dataset$offset), by = 0.005),
                          labels = function(x) {
                            labs <- sprintf("%.2f", x)
                            labs[seq_along(labs) %% 4 != 1] <- ""  # blank every 5nd label
                            labs
                          }
    ) +
    #scale_x_continuous(limits = c(-170, -52)) +
    scale_x_continuous(limits = c(min_long, max_long)) +
    scale_y_continuous(limits = c(40, 70)) +
    labs(x = "Longitude", y = "Latitude") +
    facet_wrap(year~., ncol=2) +
    theme(
      panel.background = element_rect(fill="white", color=NA),
      panel.grid = element_blank(),
      #panel.background = element_rect(fill="aliceblue", color="grey20"),
      legend.position = "right",
      legend.background = element_blank(),
      legend.key = element_blank(),
      #legend.title = element_blank(),
      strip.background = element_rect(fill=NA, color=NA),
      strip.text.x = element_text(size=12),
      #strip.text.y.left = element_text(size=16, angle = 0),  # rotate facet labe
      strip.text.y = element_blank(),
      legend.text = element_text(size=10),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    ) +
    guides(
      color = guide_colorbar(
        barwidth = 1,    # increase width (in legend units)
        barheight = 10,     # adjust height if needed
        title.position = "top",   # move title above the bar
        title.hjust = 0.5
      ))
  print(m1)
}

plot_map3(dm.CE, min_long=-110, max_long=-52)

####################################


#### Model PC3 - Central/East

dm <- dpheno1 %>% filter(cluster == "combined") %>% filter(set == "PC") %>% filter(scenario == "245") %>% 
  filter(dataset == "Provenances")
dm$group %>% unique()
head(dm)
dm.CE <- dm %>% filter(group %in% c("Central","East")) %>% filter(Longitude > -125)
max(dm.CE$offset)
min(dm.CE$offset)


#### Model var5 - West
dm.2 <- dpheno2 %>% filter(cluster == "West") %>% filter(set == "var5") %>% filter(scenario == "245") %>% 
  filter(dataset == "Provenances")
dm.2$group %>% unique()
head(dm.2)
dm.W <- dm.2 %>% filter(group %in% c("West"))
max(dm.W$offset)
min(dm.W$offset)



plot.ce <- plot_map3(dm.CE, min_long=-110, max_long=-52)
plot.w <- plot_map3(dm.W, max_long=-90)



png("03_plot_offset_separate_2.png", res=300, w=2500,h=3500)
plot_grid(plot.w, plot.ce, ncol=1, labels = c("A","B"))
dev.off()


png("03_plot_offset_separate_1.png", res=300, w=2500,h=1500)
plot.ce
dev.off()

png("03_plot_offset_separate_2.png", res=300, w=2500,h=1500)
plot.w
dev.off()



