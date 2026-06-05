rm(list=ls())
setwd("//wsl.localhost/Ubuntu/home/BlackSpruce/03_admixture_EPR_EPN")


library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape)
library(gridExtra)
library(grid)
library(RColorBrewer)
display.brewer.all(colorblindFriendly=TRUE)
display.brewer.all()
library(pals)


library(sf)
library(raster)
library(spData)
library(tmap)
library(scatterpie)





#######################
#------- DATA---------#
#######################


# Reading cluster assignement table
adm <- read.csv("admixture_04_combine_output_assigned_clusters.tsv", sep="\t", header=TRUE)
dim(adm)
adm <- adm %>% filter(SPECIES_ID != "UNK")
adm <- adm %>% filter(SPECIES_ID != "EPR") %>% filter(call_rate > 0.85)

dim(adm)
head(adm)
adm$BestK6 %>% unique()

dim(adm)



# BS distribution
my_sf <- read_sf("//wsl.localhost/Ubuntu/home/BlackSpruce/METADATA/maps/BlackSpruce/data/commondata/data0/picemari.shp")
my_sf



target_crs <- 4326
shapefile_transformed <- st_transform(my_sf, crs = target_crs)
shapefile_transformed

df.dist <- as.data.frame(shapefile_transformed)
head(df.dist)


spec_trans <- st_union(shapefile_transformed) %>%
  fortify(region = "id")
spec_trans

spec_df <- as.data.frame(st_coordinates(spec_trans))
head(spec_df)
spec_df$group <- 1
spec_df$order <- c(1:length(spec_df$X))
spec_df$region <- "black_spruce"
spec_df$id <- "black_spruce"
colnames(spec_df)[1:2] <- c("x","y")

#spec <- st_sf(geometry = shapefile_transformed)
#spec_df <- st_as_sf(st_cast(spec, "MULTIPOLYGON"))
#head(spec_df)

worldmap <- map_data("world", region = c("USA","Canada","Mexico"))
head(worldmap)
dim(worldmap)


#######################
#------- PLOT---------#
#######################


plotPie <- function(df_with_meta, colors_vector) {
  
  df_cluster_pie <- df_with_meta %>% dplyr::select(lat, lon, POP, Clusters) %>%
    group_by(lat, lon, POP, Clusters) %>% summarize(n = n()) %>% spread(key = Clusters, value = n) %>% replace(is.na(.), 0)
  
  df_clusts <- df_cluster_pie[c(4:dim(df_cluster_pie)[2])]
  selected_cols <- colnames(df_clusts)
  df_cluster_pie$n <- rowSums(df_clusts)
  df_cluster_pie$radius <- df_cluster_pie$n/max(df_cluster_pie$n) + 0.5
  
  worldmap<- map_data ("world", region = c("USA","Canada","Mexico"))
  
  p1 <- ggplot(worldmap, aes(long, lat)) +
    geom_map(map=worldmap, aes(map_id=region), fill = "floralwhite",color="darkgrey",linewidth=0.5)
    #geom_map(map=spec_df, aes(map_id=region), fill = "honeydew2",color="darkgrey",linewidth=0.5)
    #geom_polygon(data = spec_df, aes(x=x,y=y), fill = "honeydew2", color = "darkgrey", linewidth = 0.1)
  
  p <- p1 + #geom_map(map=worldmap, aes(map_id=region), fill = "white",color="black") +
    geom_scatterpie(data = df_cluster_pie, aes(x=lon, y=lat, group=POP, r = radius),
                    cols=selected_cols, color="NA", alpha=1.) +
    scale_fill_manual(values = colors_vector) +
    scale_x_continuous(limits = c(-152, -52)) +
    scale_y_continuous(limits = c(30, 70)) +
    coord_fixed() +
    labs(x = "",y="") +
    theme(panel.background = element_rect(fill=NA,color="black"),
          legend.position = "none")
  print(p)
  
}

###########
#---K=6---#
###########


head(adm)
adm_K6 <- adm %>% dplyr::select(id, POP, lat, lon, BestK6)
names(adm_K6)[names(adm_K6) == 'BestK6'] <- 'Clusters'
adm_K6$Clusters <- factor(adm_K6$Clusters, levels = c("K6.West","K6.Central","K6.East","K6.WI","K6.ME","K6.RedSpruce"))
cols_K6 <- c("#DFC27D","#C7EAE5", "#35978F", "black","goldenrod1","firebrick1")
plot.K6 <- plotPie(adm_K6, cols_K6)



png("admixture_06_plot_pie_map_K6_max.png", w=1400,h=1500,res=150)
plot.K6
dev.off()



###########
#---K=5---#
###########

adm_K5 <- adm %>% dplyr::select(id, POP, lat, lon, BestK5)
names(adm_K5)[names(adm_K5) == 'BestK5'] <- 'Clusters'
adm_K5$Clusters <- factor(adm_K5$Clusters, levels = c("K5.West","K5.Central","K5.East","K5.WI","K5.RedSpruce"))
cols_K5 <- c("#DFC27D","#C7EAE5", "#35978F", "black","firebrick1")
plot.K5 <- plotPie(adm_K5, cols_K5)


png("admixture_06_plot_pie_map_K5.png", w=2000,h=800,res=150)
plot.K5
dev.off()



###########
#---K=4---#
###########

adm_K4 <- adm %>% dplyr::select(id, POP, lat, lon, BestK4)
names(adm_K4)[names(adm_K4) == 'BestK4'] <- 'Clusters'
adm_K4$Clusters <- factor(adm_K4$Clusters, levels = c("K4.West","K4.Central","K4.East","K4.RedSpruce"))
cols_K4 <- c("#DFC27D","#C7EAE5", "#35978F", "firebrick1")
plot.K4 <- plotPie(adm_K4, cols_K4)


png("admixture_06_plot_pie_map_K4.png", w=2000,h=800,res=150)
plot.K4
dev.off()







##################
#-- call rates---#
##################

head(adm)


adm.cr <- adm %>% dplyr::select(id, POP, lat, lon, BestK6, call_rate)
head(adm.cr)
adm.cr %>% group_by(BestK6) %>% summarize(mean = mean(call_rate))

#BestK6        mean
#<chr>        <dbl>
#  1 K6.Central   0.805
#2 K6.East      0.761
#3 K6.ME        0.819
#4 K6.RedSpruce 0.813
#5 K6.West      0.832
#6 K6.WI        0.670

















