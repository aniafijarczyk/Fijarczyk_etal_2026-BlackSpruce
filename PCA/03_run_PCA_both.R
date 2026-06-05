rm(list=ls())
setwd("C:/Users/aniaf/Projects/BlackSpruce/924_PCA")

library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
library(RColorBrewer)
library(SNPRelate)
library(plotly)

brewer.pal(9, "BrBG")



#====================#
#-------DATA---------#
#====================#


#====================#
#      EPN+EPR       #
#====================#


### Metadata
dmeta <- read.csv("../DATA_intermediate/23_filter_clean_indiv_metrics.tsv", sep="\t", header=TRUE)
head(dmeta)
dim(dmeta)
dmeta %>% filter(SPECIES_ID == "EPR")

# Closing current genofile
#snpgdsClose(genofile)

# Opening gds genofile
genofile <- snpgdsOpen("../DATA_intermediate/23_filter_clean.gds")
genofile
head(genofile)
str(genofile)




#====================#
#        PCA         #
#====================#

# Running PCA
pca <- snpgdsPCA(genofile, num.thread=3, autosome.only=FALSE)

# % variance explained
pca$varprop %>% head()
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

# Getting sample ids and population info
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],
                  EV4 = pca$eigenvect[,4],
                  EV5 = pca$eigenvect[,5],
                  EV6 = pca$eigenvect[,6],
                  EV7 = pca$eigenvect[,7],
                  EV8 = pca$eigenvect[,8],
                  stringsAsFactors = FALSE)
head(tab)
tabp <- merge(tab, dmeta, by.x = 'sample.id', by.y = 'id', sort = FALSE)
head(tabp)

unique(tabp$cluster)
tabp$cluster <- factor(tabp$cluster, levels = c("West","Central","East","WI","ME","RedSpruce"))

#====================#
#       Plots        #
#====================#

df.labs <- data.frame("label" = c("East","Central","West","WI","ME","RedSpruce"),
                      "x" = c(-0.06, 0, -0.02, 0.03, -0.09, -0.12),
                      "y" = c(-0.02, 0.003, 0.05, -0.022, 0, 0.06))

df.line <- data.frame(y1 = -0.019, y2 = -0.022, x1 = 0, x2 = 0.02)
df.me <- data.frame(y1 = 0, y2 = 0, x1 = -0.05, x2 = -0.07)

pca.plot_12 <- ggplot(tabp) +
  
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2),colour = "black", data = df.line) +
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2),colour = "black", data = df.me) +
  
  geom_point(data=tabp,  aes(x = EV1, y = EV2, fill = cluster), 
             size=3, pch=21, alpha=0.95) +
  
  scale_fill_manual(values = c("#DFC27D","#C7EAE5","#35978F","black","goldenrod1","firebrick1"), name="Cluster") +
  labs(x = paste0("PC1 = ",round(pc.percent, 2)[1],"%"), y = paste0("PC2 = ",round(pc.percent, 2)[2],"%")) +
  #scale_x_continuous(limits=c(-0.21, 0.1)) +
  geom_text(data = df.labs, aes(x=x,y=y, label=label), size=5) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        #legend.position = c(0.2,0.2),
        legend.position = "none",
        panel.border = element_rect(fill = FALSE),
        axis.text = element_text(size=10),
        axis.title = element_text(size=14)
  )
pca.plot_12

png("03_run_PCA_both_12.png", w=1100, h=900, res=300)
plot_grid(pca.plot_12)
dev.off()



##########################################
### PC1 vs. PC2 by strata
pca.plot__12 <- ggplot(tabp) + aes(x = EV1, y = EV2, fill = STRATA) +
  geom_point(size=3, pch=21) +
  scale_fill_manual(values = c(brewer.pal(9, "BrBG"),"firebrick1")) +
  labs(x = paste0("PC1 = ",round(pc.percent, 2)[1],"%"), y = paste0("PC3 = ",round(pc.percent, 2)[3],"%")) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        panel.border = element_rect(fill = FALSE),
        axis.title = element_text(size=20)
  )
pca.plot__12

### PC1 vs. PC3 by strata
pca.plot_13 <- ggplot(tabp) + aes(x = EV1, y = EV3, fill = STRATA) +
  geom_point(size=3, pch=21) +
  scale_fill_manual(values = c(brewer.pal(9, "BrBG"),"firebrick1")) +
  labs(x = paste0("PC1 = ",round(pc.percent, 2)[1],"%"), y = paste0("PC3 = ",round(pc.percent, 2)[3],"%")) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        panel.border = element_rect(fill = FALSE),
        axis.title = element_text(size=20)
  )
pca.plot_13


# PC3 vs. PC4 by population
pca.plot_34 <- ggplot(tabp) + aes(x = EV3, y = EV4, fill = STRATA) +
  geom_point(size=3, pch=21) +
  scale_fill_manual(values = c(brewer.pal(9, "BrBG"),"firebrick1")) +
  labs(x = paste0("PC3 = ",round(pc.percent, 2)[3],"%"), y = paste0("PC4 = ",round(pc.percent, 2)[4],"%")) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = FALSE),
        legend.position = "none",
        axis.title = element_text(size=20))
pca.plot_34

# PC4 vs. PC5 by population
pca.plot_45 <- ggplot(tabp) + aes(x = EV4, y = EV5, fill = STRATA) +
  geom_point(size=3, pch=21) +
  scale_fill_manual(values = c(brewer.pal(9, "BrBG"),"firebrick1")) +
  labs(x = paste0("PC3 = ",round(pc.percent, 2)[3],"%"), y = paste0("PC4 = ",round(pc.percent, 2)[4],"%")) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = FALSE),
        legend.position = "none",
        axis.title = element_text(size=20))
pca.plot_45




png("03_run_PCA_both_12345.png", h=3000, w=3000, res=300)
plot_grid(pca.plot__12, pca.plot_13, pca.plot_34,pca.plot_45, ncol=2, nrow=2)
dev.off()


