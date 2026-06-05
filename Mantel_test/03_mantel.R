rm(list=ls())
setwd("C:/Users/aniaf/Projects/BlackSpruce/925_Mantel")


library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)
library(sf)

library(RColorBrewer)
library(readxl)
library(units)

library(ade4) # mantel.rtest
library(vegan) # mantel



#######################
#------- DATA---------#
#######################



getDist <- function(filename) {
  
  dgeo <- read.csv(filename, header=T)
  pops <- dgeo$X
  dgeo <- dgeo %>% dplyr::select(-X)
  mat.geo <- as.matrix(dgeo)
  rownames(mat.geo) <- pops
  colnames(mat.geo) <- pops
  dist.geo <- as.dist(mat.geo)
  return(dist.geo)

  }

dist.geo <- getDist("02_inputs_geo_matrix.csv")
dist.geo

dist.fst <- getDist("02_inputs_transfst_matrix.csv")
dist.fst





#########################
#-----Mantel test-------#
#########################


mtest <- mantel.rtest(dist.geo, dist.fst, nrepet = 9999)
mtest

mtest.2 <- mantel(dist.geo, dist.fst, method = "spearman", permutations = 9999, na.rm = TRUE)
mtest.2

mtest.2$statistic
mtest.2$signif
mtest.2$permutation

mantel.test <- data.frame(
  "lab" = c("Mantel's r","p-value"),
  "val" = c(round(mtest.2$statistic,2), mtest.2$signif),
  "set" = rep("all",2)
  )
mantel.test


runMantel <- function(mat.1, mat.2, setname) {
  
  mtest.2 <- mantel(mat.1, mat.2, method = "spearman", permutations = 9999, na.rm = TRUE)
  mantel.test <- data.frame(
    "MantelsR" = c(round(mtest.2$statistic,2)),
    "Pval" = c(mtest.2$signif),
    "Set" = c(setname)
  )
  return(mantel.test)

}


#########################
#       A LOOP          #
#########################


fnames <- list.files(path="./",pattern="02_inputs_*.+csv")
fnames


M <- list()
for (fset in c(1:7)) {
  print(fset)
  dist.geo <- getDist(paste0("02_inputs_geo_matrix_strata_",fset,".csv"))
  dist.fst <- getDist(paste0("02_inputs_transfst_matrix_strata_",fset,".csv"))
  mant.test <- runMantel(dist.geo, dist.fst, fset)
  M[[fset]] <- mant.test
} 
M

dist.geo <- getDist("02_inputs_geo_matrix.csv")
dist.fst <- getDist("02_inputs_transfst_matrix.csv")
M[[8]] <- runMantel(dist.geo, dist.fst, 0)
M

all_results <- bind_rows(M)
all_results$adjusted.p <- p.adjust(all_results$Pval)
all_results

set_names <- data.frame("Set" = c(0:7),
                        "SetName" = c("all","ABC","BCD","CDE","DEF","EFG","FGH","GHI"))
mall_results <- merge(all_results, set_names, by = "Set", sort=F)
mall_results


write.table(mall_results, "03_mantel_results.tsv", sep="\t", col.names = T, row.names = F, quote=F, append=F)















#########################
#---------Plots---------#
#########################



# Function to convert distance matrix into data frame with no repeated pairwise measures and no diagonal values
matrix2long <- function(mat) {
  gd <- as.data.frame(as.matrix(as.dist(mat))) 
  gd$pop1 <- colnames(gd)
  gd2 <- gd %>% gather(key = "pop2", value = "dist", -pop1)
  gd2$pop1 <- as.numeric(gd2$pop1)
  gd2$pop2 <- as.numeric(gd2$pop2)
  gd2$P1 <- apply(gd2[c(1,2)],1, function(x) min(x))
  gd2$P2 <- apply(gd2[c(1,2)],1, function(x) max(x))
  gd3 <- gd2 %>% dplyr::select(-pop1, -pop2) %>% distinct() %>% filter(P1 != P2)
  return(gd3)
}


################################################################################
### All

# Getting long data frames
dist.geo <- getDist("02_inputs_geo_matrix.csv")
dist.fst <- getDist("02_inputs_transfst_matrix.csv")

gd3 <- matrix2long(dist.geo)
fd3 <- matrix2long(dist.fst)

df <- merge(gd3, fd3, by = c("P1","P2"), sort=FALSE)
head(df)

group <- "all"
mantel.result <- mall_results[mall_results$SetName==group,]
mantel.result

x.pos <- min(dist.geo)
y.pos <- max(dist.fst)

mantel.test <- data.frame("x" = c(x.pos, x.pos),
                          "y" = c(y.pos, y.pos + 0.003),
                          "lab" = c("Mantel's r","p-value"),
                          "val" = c(mantel.result$MantelsR, mantel.result$Pval))
mantel.test


p1 <- ggplot(df) + aes(x = dist.x, y = dist.y) +
  geom_point(pch=21) +
  #geom_smooth(method = "lm", colour = "grey20") +
  geom_text(data = mantel.test, aes(x=x, y=y, label=paste0(lab,"=",val)), 
            color = "grey20", hjust=0, size=6) +
  labs(x = "Distance [km]", y = "FST/(1-FST)") +
  theme(panel.background = element_rect(fill=NA, colour="grey20"),
        panel.grid = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size=18))
p1  

png("03_mantel_all.png", w = 800, h = 600, res = 150)
p1
dev.off()






######### Combine long formats

L <- list()
for (fset in c(1:7)) {
  print(fset)
  dist.geo <- getDist(paste0("02_inputs_geo_matrix_strata_",fset,".csv"))
  dist.fst <- getDist(paste0("02_inputs_transfst_matrix_strata_",fset,".csv"))
  gd3 <- matrix2long(dist.geo)
  fd3 <- matrix2long(dist.fst)
  df <- merge(gd3, fd3, by = c("P1","P2"), sort=FALSE)
  df$Set <- fset
  L[[fset]] <- df
} 
L

dist.geo <- getDist("02_inputs_geo_matrix.csv")
dist.fst <- getDist("02_inputs_transfst_matrix.csv")
gd3 <- matrix2long(dist.geo)
fd3 <- matrix2long(dist.fst)
df <- merge(gd3, fd3, by = c("P1","P2"), sort=FALSE)
df$Set <- 0
L[[8]] <- df


all_frames <- bind_rows(L)
set_names <- data.frame("Set" = c(0:7),
                        "SetName" = c("all","ABC","BCD","CDE","DEF","EFG","FGH","GHI"))
mall_frames <- merge(all_frames, set_names, by = "Set", sort=F)
mall_frames






#########################
#      Plot grid        #
#########################


head(mall_frames)
head(mall_results)

max_dist <- mall_frames %>% group_by(Set) %>% summarize(fmax = max(dist.y))
mmall_results <- merge(mall_results, max_dist, by="Set", sort=F)
head(mmall_results)

mall_frames$SetName <- factor(mall_frames$SetName, levels=c("all","ABC","BCD","CDE","DEF","EFG","FGH","GHI"))
mmall_results$SetName <- factor(mmall_results$SetName, levels=c("all","ABC","BCD","CDE","DEF","EFG","FGH","GHI"))

p2 <- ggplot(mall_frames) + 
  geom_point(data = mall_frames, aes(x = dist.x, y = dist.y), pch=21) +
  #geom_smooth(method = "lm", colour = "grey20") +
  geom_text(data = mmall_results, aes(x=0, y=fmax, label=paste0("Mantel's R=",MantelsR,";\nP=",round(adjusted.p,3))), 
            color = "grey20", hjust=0, size=4, vjust=1) +
  facet_wrap(~SetName, scales="free", ncol=2) +
  labs(x = "Distance [km]", y = "FST/(1-FST)") +
  theme(panel.background = element_rect(fill=NA, colour="grey20"),
        panel.grid = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size=18))
p2  


png("03_mantel_grid.png", w = 1200, h = 1400, res = 150)
p2
dev.off()


