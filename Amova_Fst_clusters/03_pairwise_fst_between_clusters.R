rm(list=ls())
setwd("C:/Users/aniaf/Projects/BlackSpruce/928_amova")

library(ggplot2)
library(tidyr)
library(dplyr)
library(dartR)
library(tibble)
library(sf)
library(RColorBrewer)
library(readxl)
library(assigner) # calculating pairwise WC FST




#===================#
#       DATA        #
#===================#

### Genotypes of EPN

# Open genlight
load("../DATA_intermediate/23_filter_EPN_genlight.Rdata")
genl
class(genl)



# Substitute pop with a cluster name
genl$pop
genl@other$cluster
genl$pop <- as.factor(genl@other$cluster)
genl$pop
genl


#===================#
#    Data prep      #
#===================#

# Converting genlight to tidy
gel.tidy <- radiator::tidy_genlight(genl)
gel.tidy$GT_BIN
gel.tidy$GT <- gel.tidy$GT_BIN
head(gel.tidy)
dim(gel.tidy)



#===================#
#   Pairwise Fst    #
#===================#


### Genetic cluster


cluster.fst.wc <- assigner::fst_WC84(gel.tidy,
                           pairwise = TRUE,
                           ci = TRUE,
                           heatmap.fst = FALSE)



# Checking results
head(cluster.fst.wc$pairwise.fst)

cluster.fst.wc$pairwise.fst.full.matrix
cluster.fst.wc$fst.plot
dim(cluster.fst.wc$pairwise.fst.full.matrix)
colnames(cluster.fst.wc$pairwise.fst.full.matrix)


# Saving st results
save(cluster.fst.wc, file = "03_pairwise_fst_between_clusters.RData")





############################################
#          READING PAIRWISE FST            #
############################################

# Read fst output if necessary
#load("03_pairwise_fst_between_clusters.RData")



#===================#
#      Heatmap      #
#===================#

# Getting pairs in order from the full matrix

df.mat <- as.data.frame(cluster.fst.wc$pairwise.fst.full.matrix)
df.mat$POP1 <- row.names(df.mat)
df.long <- df.mat %>% gather(key = "POP2", value = "FST", -POP1)
df.long$FST <- as.numeric(df.long$FST)
head(df.long)
df.order.1 <- data.frame("POP1" = c("West","Central","East","WI","ME"),
                       "lp1" = c(1,2,3,4,5))
df.order.2 <- data.frame("POP2" = c("West","Central","East","WI","ME"),
                       "lp2" = c(1,2,3,4,5))

dm.long <- merge(df.long, df.order.1, by = "POP1", sort=F, all.x=T)
dm.long.2 <- merge(dm.long, df.order.2, by = "POP2", sort=F, all.x=T)
dm.table <- dm.long.2 %>% arrange(lp1, lp2) 
head(dm.table)

# Plotting

dm.table$POP1 <- factor(dm.table$POP1, levels = c("West","Central","East","WI","ME"))
dm.table$POP2 <- factor(dm.table$POP2, levels = c("West","Central","East","WI","ME"))
dm.table

# Selecting pairs
dm.table.sub <- dm.table %>% 
  filter(POP2 != "West") %>%
  filter(POP1 != "ME") %>%
  filter(!(POP1 == "WI" & POP2  %in% c("East","Central"))) %>%
  filter(!(POP1 == "East" & POP2  == "Central"))
dm.table.sub$FST <- ifelse(dm.table.sub$POP1 == dm.table.sub$POP2, NA, dm.table.sub$FST)
dm.table.sub <- na.omit(dm.table.sub)
dm.table.sub



p1 <- ggplot(dm.table.sub) + aes(x = POP1, y = POP2, fill = FST) + 

  geom_tile(color="white") +
  geom_text(aes(label=round(FST,3)), color="black", size=6) +

    coord_equal() +
  scale_fill_gradient(low = "white", high = "grey", name="FST") +
  scale_x_discrete(position = "top") +
  
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(size=20))
p1


png("03_pairwise_fst_between_clusters.png", w=1400,h=1400,res=300)
p1
dev.off()


