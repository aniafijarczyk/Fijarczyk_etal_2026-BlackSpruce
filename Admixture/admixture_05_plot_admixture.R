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


#######################
#------- DATA---------#
#######################


#######################
#------- PLOT---------#
#######################


################################################################################
### All clusters

head(adm.sorted)
#adm.sorted %>% dplyr::select(STRATA, )

lon_ordered_samples <- adm.sorted %>% filter(Cluster == "K2.Q0") %>% pull(sample)
lon_ordered_samples
length(lon_ordered_samples)
length(unique(lon_ordered_samples))

adm.sorted$sample <- factor(adm.sorted$sample, levels = lon_ordered_samples)

#adm.sorted$K <- factor(adm.sorted$K, levels = c("K1","K2","K3","K4","K5","K6","K7","K8","K9","K10","K11","K12"))
#adm.sorted$Q <- factor(adm.sorted$Q, levels = c("Q0","Q1","Q2","Q3","Q4","Q5","Q6","Q7","Q8","Q9","Q10","Q11"))

adm.sorted$K <- factor(adm.sorted$K, levels = c("K1","K2","K3","K4","K5","K6","K7","K8","K9","K10","K11","K12","K13","K14","K15"))
adm.sorted$Q <- factor(adm.sorted$Q, levels = c("Q0","Q1","Q2","Q3","Q4","Q5","Q6","Q7","Q8","Q9","Q10","Q11","Q12","Q13","Q14"))

set.seed(002)
cols.poly <- sample(as.vector(palette.colors(palette = "Polychrome 36")))
set.seed(003)
cols.alph <- sample(as.vector(palette.colors(palette = "Alphabet")))

p1 <- ggplot(adm.sorted) + aes(x = sample, y = Proportion, fill = Q) +
  geom_bar(position="fill", stat="identity",width=1) +
  facet_wrap(~K, ncol=1, strip.position="right") +
  scale_fill_manual(values = cols.poly) +
  labs(x = "Samples West -> East") +
  theme(panel.background = element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p1




png("admixture_05_plot_admixture_all_K.png", w=2000,h=1200,res=150)
p1
dev.off()



################################################################################
### Clusters 2-9

adm.sorted$group <- ifelse(adm.sorted$SPECIES_ID == 'EPN', adm.sorted$STRATA, adm.sorted$SPECIES_ID)
head(adm.sorted)
adm.sorted <- adm.sorted %>% arrange(group, lon)

lon_ordered_samples <- adm.sorted %>% filter(Cluster == "K2.Q0") %>% pull(sample)
adm.sorted$sample <- factor(adm.sorted$sample, levels = lon_ordered_samples)
adm.selected <- adm.sorted %>% filter(K %in% c("K2","K3","K4","K5","K6","K7","K8","K9"))
adm.selected$K <- factor(adm.selected$K, levels = c("K2","K3","K4","K5","K6","K7","K8","K9"))
adm.selected$Q <- factor(adm.selected$Q, levels = c("Q0","Q1","Q2","Q3","Q4","Q5","Q6","Q7","Q8","Q9","Q10"))
adm.selected$group <- factor(adm.selected$group, levels = c("A","B","C","D","E","F","G","H","I","EPR","UNK"))


set.seed(002)
cols.poly <- sample(as.vector(palette.colors(palette = "Polychrome 36")))
set.seed(016)
cols.alph <- sample(as.vector(palette.colors(palette = "Alphabet")))


p2 <- ggplot(adm.selected) + aes(x = sample, y = Proportion, fill = Q) +
  geom_bar(position="fill", stat="identity",width=1) +
  #facet_wrap(~K, strip.position="right", ncol=1) +
  facet_grid(K~group, scales = "free_x", space='free') +
  #facet_grid(group~K) +
  scale_fill_manual(values = as.vector(palette.colors())) +
  labs(x = "Samples West -> East") +
  theme(panel.background = element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p2



png("admixture_05_plot_admixture_K2_9.png", w=2000,h=1000,res=150)
p2
dev.off()



################################################################################
### Clusters 3-6

adm.sorted <- adm.sorted %>% filter(SPECIES_ID != "UNK")
adm.sorted$group <- ifelse(adm.sorted$SPECIES_ID == 'EPN', adm.sorted$STRATA, adm.sorted$SPECIES_ID)
head(adm.sorted)
adm.sorted <- adm.sorted %>% arrange(group, lon)

lon_ordered_samples <- adm.sorted %>% filter(Cluster == "K2.Q0") %>% pull(sample)
adm.sorted$sample <- factor(adm.sorted$sample, levels = lon_ordered_samples)
adm.selected <- adm.sorted %>% filter(K %in% c("K3","K4","K5","K6"))
adm.selected$K <- factor(adm.selected$K, levels = c("K3","K4","K5","K6"))
adm.selected$Q <- factor(adm.selected$Q, levels = c("Q0","Q1","Q2","Q3","Q4","Q5"))
adm.selected$group <- factor(adm.selected$group, levels = c("A","B","C","D","E","F","G","H","I","EPR","UNK"))


set.seed(002)
cols.poly <- sample(as.vector(palette.colors(palette = "Polychrome 36")))
set.seed(016)
cols.alph <- sample(as.vector(palette.colors(palette = "Alphabet")))

head(adm.selected)


brewer.pal(9, "BrBG")
brewer.pal(11, "Spectral")

colors_K9 <- c("#8C510A","#BF812D", "#DFC27D", "#F6E8C3", "#F5F5F5", "#C7EAE5", "#80CDC1", "#35978F", "#01665E", "firebrick1")

cols_K3to6 <- c("#35978F", "#DFC27D", "firebrick1",
                "#C7EAE5", "firebrick1", "#35978F", "#DFC27D",
                "black", "#C7EAE5","#DFC27D","#35978F", "firebrick1",
                "goldenrod1", "#35978F", "#C7EAE5","#DFC27D","black","firebrick1")

p2 <- ggplot(adm.selected) + aes(x = sample, y = Proportion, fill = Cluster) +
  geom_bar(position="fill", stat="identity",width=1) +
  #facet_wrap(~K, strip.position="right", ncol=1) +
  facet_grid(K~group, scales = "free_x", space='free') +
  #facet_grid(group~K) +
  scale_fill_manual(values = cols_K3to6) +
  labs(x = "Samples West -> East") +
  theme(panel.background = element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")
p2



png("admixture_05_plot_admixture_K3_6.png", w=3000,h=1000,res=150)
p2
dev.off()



################################################################################
### K 6 - clusters

# Removing UNK
adm.sorted <- adm.sorted %>% filter(SPECIES_ID != "UNK")
adm.sorted$group <- ifelse(adm.sorted$SPECIES_ID == 'EPN', adm.sorted$STRATA, adm.sorted$SPECIES_ID)
adm.sorted <- adm.sorted %>% arrange(group, lon)
lon_ordered_samples <- adm.sorted %>% filter(Cluster == "K2.Q0") %>% pull(sample)
adm.sorted$sample <- factor(adm.sorted$sample, levels = lon_ordered_samples)
adm.sorted$K <- factor(adm.sorted$K, levels = c("K1","K2","K3","K4","K5","K6","K7","K8","K9","K10","K11","K12","K13","K14","K15"))
adm.sorted$Q <- factor(adm.sorted$Q, levels = c("Q0","Q1","Q2","Q3","Q4","Q5","Q6","Q7","Q8","Q9","Q10","Q11","Q12","Q13","Q14"))
adm.selected <- adm.sorted %>% filter(K %in% c("K6"))
adm.selected$Q <- factor(adm.selected$Q,
                         levels = c("Q0","Q1","Q2","Q3","Q4","Q5"),
                         labels = c("Q0.ME","Q1.East","Q2.Central","Q3.West","Q4.WI","Q5.Red Spruce"))
adm.selected$group <- factor(adm.selected$group, levels = c("A","B","C","D","E","F","G","H","I","EPR"))

cols_K6 <- c("goldenrod1", "#35978F", "#C7EAE5","#DFC27D","black","firebrick1")

p3 <- ggplot(adm.selected) + aes(x = sample, y = Proportion, fill = Q) +
  geom_bar(position="fill", stat="identity",width=1) +
  facet_grid(K~group, scales = "free_x", space='free') +
  #scale_fill_manual(values = as.vector(palette.colors())) +
  scale_fill_manual(values = cols_K6) +
  labs(x = "Samples West -> East") +
  theme(panel.background = element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")
p3

png("admixture_05_plot_admixture_K6.png", w=2500,h=400,res=150)
p3
dev.off()


################################################################################
### K 5 - clusters

# Removing UNK
adm.sorted <- adm.sorted %>% filter(SPECIES_ID != "UNK")
adm.sorted$group <- ifelse(adm.sorted$SPECIES_ID == 'EPN', adm.sorted$STRATA, adm.sorted$SPECIES_ID)
adm.sorted <- adm.sorted %>% arrange(group, lon)
lon_ordered_samples <- adm.sorted %>% filter(Cluster == "K2.Q0") %>% pull(sample)
adm.sorted$sample <- factor(adm.sorted$sample, levels = lon_ordered_samples)
adm.sorted$K <- factor(adm.sorted$K, levels = c("K1","K2","K3","K4","K5","K6","K7","K8","K9","K10","K11","K12","K13","K14","K15"))
adm.sorted$Q <- factor(adm.sorted$Q, levels = c("Q0","Q1","Q2","Q3","Q4","Q5","Q6","Q7","Q8","Q9","Q10","Q11","Q12","Q13","Q14"))
adm.selected <- adm.sorted %>% filter(K %in% c("K5"))
adm.selected$Q <- factor(adm.selected$Q,
                         levels = c("Q0","Q1","Q2","Q3","Q4"),
                         labels = c("Q0.WI","Q1.Central","Q2.West","Q3.East","Q4.Red Spruce"))
adm.selected$group <- factor(adm.selected$group, levels = c("A","B","C","D","E","F","G","H","I","EPR"))

cols_K5 <- c("black", "#C7EAE5","#DFC27D","#35978F", "firebrick1")

p3 <- ggplot(adm.selected) + aes(x = sample, y = Proportion, fill = Q) +
  geom_bar(position="fill", stat="identity",width=1) +
  facet_grid(K~group, scales = "free_x", space='free') +
  #scale_fill_manual(values = as.vector(palette.colors())) +
  scale_fill_manual(values = cols_K5) +
  labs(x = "Samples West -> East") +
  theme(panel.background = element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")
p3

png("admixture_05_plot_admixture_K5.png", w=2500,h=400,res=150)
p3
dev.off()



################################################################################
### K 4 - clusters

# Removing UNK
adm.sorted <- adm.sorted %>% filter(SPECIES_ID != "UNK")
adm.sorted$group <- ifelse(adm.sorted$SPECIES_ID == 'EPN', adm.sorted$STRATA, adm.sorted$SPECIES_ID)
adm.sorted <- adm.sorted %>% arrange(group, lon)
lon_ordered_samples <- adm.sorted %>% filter(Cluster == "K2.Q0") %>% pull(sample)
adm.sorted$sample <- factor(adm.sorted$sample, levels = lon_ordered_samples)
adm.sorted$K <- factor(adm.sorted$K, levels = c("K1","K2","K3","K4","K5","K6","K7","K8","K9","K10","K11","K12","K13","K14","K15"))
adm.sorted$Q <- factor(adm.sorted$Q, levels = c("Q0","Q1","Q2","Q3","Q4","Q5","Q6","Q7","Q8","Q9","Q10","Q11","Q12","Q13","Q14"))
adm.selected <- adm.sorted %>% filter(K %in% c("K4"))
adm.selected$Q <- factor(adm.selected$Q,
                         levels = c("Q0","Q3","Q2","Q1"),
                         labels = c("Q0.Central","Q3.West","Q2.East","Q1.Red Spruce"))
adm.selected$group <- factor(adm.selected$group, levels = c("A","B","C","D","E","F","G","H","I","EPR"))

cols_K4 <- c("#C7EAE5", "#DFC27D", "#35978F", "firebrick1")

p4 <- ggplot(adm.selected) + aes(x = sample, y = Proportion, fill = Q) +
  geom_bar(position="fill", stat="identity",width=1) +
  facet_grid(K~group, scales = "free_x", space='free') +
  #scale_fill_manual(values = as.vector(palette.colors())) +
  scale_fill_manual(values = cols_K4) +
  labs(x = "Samples West -> East") +
  theme(panel.background = element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")
p4

png("admixture_05_plot_admixture_K4.png", w=2500,h=400,res=150)
p4
dev.off()


adm.k5 <- adm.sorted %>% filter((K %in% c("K5")) & (lon > -75) & (lon < -67))
longitudes <- adm.k5 %>% filter(Cluster %in% c("K5.Q0")) %>% pull(pop)





















################################################################################






############# Group 1
adm.G1 <- adm.sorted %>% filter((K %in% c("K5")) & (pop %in% c(7000,6994,6988,6986,6999,6983,6979,6970)))
p.G1 <- ggplot(adm.G1) + aes(x = sample, y = Proportion, fill = Q) +
  geom_bar(position="fill", stat="identity",width=1) +
  scale_fill_manual(values = as.vector(palette.colors())) +
  theme(panel.background = element_blank(),
        strip.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")
p.G1
png("admixture_05_plot_admixture_group1.png", w=400,h=300,res=150)
p.G1
dev.off()


############# Group 2
adm.G2 <- adm.sorted %>% filter((K %in% c("K5")) & (pop %in% c(6965,6967,6968)))
p.G2 <- ggplot(adm.G2) + aes(x = sample, y = Proportion, fill = Q) +
  geom_bar(position="fill", stat="identity",width=1) +
  scale_fill_manual(values = as.vector(palette.colors())) +
  theme(panel.background = element_blank(),
        strip.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")
p.G2
png("admixture_05_plot_admixture_group2.png", w=400,h=300,res=150)
p.G2
dev.off()


############# Group 3
adm.G3 <- adm.sorted %>% filter((K %in% c("K5")) & (pop %in% c(3268)))
p.G3 <- ggplot(adm.G3) + aes(x = sample, y = Proportion, fill = Q) +
  geom_bar(position="fill", stat="identity",width=1) +
  scale_fill_manual(values = as.vector(palette.colors())) +
  theme(panel.background = element_blank(),
        strip.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")
p.G3
png("admixture_05_plot_admixture_group3.png", w=300,h=300,res=150)
p.G3
dev.off()


############# Group 4
adm.G4 <- adm.sorted %>% filter((K %in% c("K5")) & (pop %in% c(321)))
p.G4 <- ggplot(adm.G4) + aes(x = sample, y = Proportion, fill = Q) +
  geom_bar(position="fill", stat="identity",width=1) +
  scale_fill_manual(values = as.vector(palette.colors())) +
  theme(panel.background = element_blank(),
        strip.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")
p.G4
png("admixture_05_plot_admixture_group4.png", w=300,h=300,res=150)
p.G4
dev.off()


############# Group 5
adm.G5 <- adm.sorted %>% filter((K %in% c("K5")) & (pop %in% c(4360)))
p.G5 <- ggplot(adm.G5) + aes(x = sample, y = Proportion, fill = Q) +
  geom_bar(position="fill", stat="identity",width=1) +
  scale_fill_manual(values = as.vector(palette.colors())) +
  theme(panel.background = element_blank(),
        strip.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")
p.G5
png("admixture_05_plot_admixture_group5.png", w=300,h=300,res=150)
p.G5
dev.off()





adm <- read.csv("admixture_04_combine_output_assigned_clusters.tsv", sep="\t", header=TRUE)
head(adm)
dim(adm)

adm %>% head()
adm %>% filter(BestK6 %in% c("K6.WI")) %>% dim()
adm %>% filter(BestK6 %in% c("K6.ME")) %>% dim()
pops.weird <- adm %>% filter(BestK6 %in% c("K6.WI","K6.ME"))
dim(pops.weird)
pops.weird.selected <- pops.weird %>% dplyr::select("sample","id","SPECIES_ID","SITE_ID","BATCH_ID","POP","STRATA","lat","lon","BestK6","BestK6_Q","call_rate") %>% filter(call_rate > 0.85) %>% arrange(POP)


write.table(pops.weird.selected, "Samples_genetic_structure_outliers.tsv", sep="\t", row.names = FALSE, col.names = TRUE, append=FALSE, quote=FALSE)
