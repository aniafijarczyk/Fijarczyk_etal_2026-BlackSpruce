rm(list=ls())
setwd("C:/Users/aniaf/Projects/BlackSpruce/04_maps")


library(ggplot2)
library(tidyr)
library(dplyr)
library(dartR)
library(ggrepel)

library(sf)
library(raster)
library(spData)
library(tmap)


library(scatterpie)
library(RColorBrewer)
library(readxl)



#######################
#------- DATA---------#
#######################

library(sf)
my_sf <- read_sf("C:/Users/aniaf/Projects/BlackSpruce/METADATA/maps/BlackSpruce/data/commondata/data0/picemari.shp")
my_sf

head(my_sf)

species_distribution_outline <- st_union(my_sf)


plants <- data.frame("Peace_River" = c(56.00, -116.39, "PR"),
                     "Mont_Laurier" = c(46.60, -75.80, "ML"),
                     "Chibougamau" = c(50.02, -74.21,"CHI"),
                     "Acadia" = c(46.00, -66.38,"AC"),
                     "Chapleau" = c(47.963757,-83.433724, "CHA"),
                     "Petawawa" = c(46.00, -77.50, "PE"),
                     "Roddickton" = c(50.90,-56.10, "RO"))

dplants <- as.data.frame(t(plants))
colnames(dplants) <- c("lat", "lon", "code")
head(dplants)
dplants$set <- "gardens"
dplants$lat <- as.numeric(dplants$lat)
dplants$lon <- as.numeric(dplants$lon)
dplants

### provenances
dp <- read.csv("../METADATA/SampleSheet_DARTseq_BILAN_08-2022.csv", sep=",", header=TRUE)
head(dp)
dpop <- dp %>% dplyr::select(POP_ID, lat, lon, POP_GR) %>% distinct()
head(dpop)
dim(dpop)
dpop$code <- dpop$POP_ID
dpop$set <- "provenances"
head(dpop)
row.names(dpop) <- dpop$POP_ID

### North America -  datasets in the spData package
world_NA = world[world$continent == "North America", ]
world_NA$geom
head(world_NA)


df_coords <- st_as_sf(dpop, coords = c("lat", "lon"), crs = 4326)
head(df_coords)

#####################
#----- PLOTS -------#
#####################


###################################################
### Plotting provenances for all sequenced data

head(dpop)
head(dplants)

p1 <- ggplot() + 
  geom_sf(data = world_NA, fill = "floralwhite") +
  geom_sf(data = species_distribution_outline, fill = "darkseagreen3") +
  geom_point(data = dpop, aes(x=lon, y=lat), shape = 19, size = 2, color="grey10") +
  geom_point(data = dplants, aes(x=lon, y=lat), shape = 24, size = 2, fill = "red", color="white") +

  #geom_label(data = dplants, aes(x=lon, y=lat, label = code)) +
  geom_label_repel(data = dplants, aes(x=lon, y=lat, label = code),
                   min.segment.length = unit(0, 'lines'), segment.size=0.5) +

  scale_x_continuous(limits = c(-170, -52)) +
  scale_y_continuous(limits = c(40, 70)) +
  labs(x = "Longitude", y = "Latitude") +
  theme(panel.background = element_rect(fill="aliceblue", color="grey20"),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        legend.position = c(0.15, 0.1),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=20),
  )
p1




png("maps_03_gardens_provs.png", w=1400,h=1000, res=150,type="cairo")
p1
dev.off()



################################################################################


dc <- read.csv("../23_filter/01_remove_mixed_admixture_metadata_clean.tsv", sep="\t", header=T)
dc <- dc %>% filter(SPECIES_ID == "EPN") %>% dplyr::select(POP, STRATA, lat, lon) %>% distinct()
head(dc)
dim(dc)


head(dplants)
dplants.sub <- dplants %>% filter(code %in% c("PR","ML","AC","CHI"))
dplants.sub$code.2 <- dplants.sub$code
dplants.sub$code.2[dplants.sub$code.2 == "CHI"] <- "CH"
dplants.sub
dplants.sub$x <- c(-121, -80, -72, -60)
dplants.sub$y <- c(54, 43, 54, 43)

head(dpop)

dpop.sub <- dpop %>% filter(POP_ID %in% dc$POP) %>% dplyr::select(lat, lon, POP_GR, set, POP_ID)
#dpop.sub <- dpop.sub %>% dplyr::rename(code = POP_GR, code.2 = POP_ID)
head(dpop.sub)
dim(dpop.sub)

p2 <- ggplot() + 
  geom_sf(data = world_NA, fill = "floralwhite") +
  geom_sf(data = species_distribution_outline, fill = "white") +
  geom_point(data = dpop.sub, aes(x=lon, y=lat, fill = POP_GR), shape = 21, size = 3, color="grey10") +
  scale_fill_manual(values = c(brewer.pal(9, "BrBG"),"firebrick1"), name="Region") +
  
  geom_segment(data = dplants.sub,
               aes(x = lon, xend = x, y = lat, yend = y), linewidth=0.5) +
  
  geom_point(data = dplants.sub, aes(x=lon, y=lat), shape = 24, size = 4, fill = "black", color="white") +

  geom_label(data = dplants.sub,
            aes(x = x, y = y, label = code.2), color = "black", size=4) +
  
  scale_x_continuous(limits = c(-170, -52)) +
  scale_y_continuous(limits = c(40, 70)) +
  labs(x = "Longitude", y = "Latitude") +
  theme(panel.background = element_rect(fill="aliceblue", color="grey20"),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size=12),
        legend.position = "right",
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_text(size=14),
        legend.text = element_text(size=14),
  ) +
  guides(fill = guide_legend(override.aes = list(size = 5)))
p2




png("maps_03_gardens_groups_4gardens.png", w=2800,h=2000, res=300,type="cairo")
p2
dev.off()

