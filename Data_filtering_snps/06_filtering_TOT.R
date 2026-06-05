rm(list=ls())
#setwd("C:/Users/aniaf/Projects/BlackSpruce/01_snp_treatment")
dir()


library(tidyverse)
library(dartR)
library(cowplot)
library(readxl)
library(stringr)



############
### DATA ###
############

gel <- gl.load("../DATA_intermediate/04_merging.Rdata")




####################
### BLACK SPRUCE ###
####################

# Black spruce
gl.EPN <- gl.keep.pop(gel, c("EPN"), as.pop="SPECIES_ID")
gl.EPN <- gl.filter.monomorphs(gl.EPN)
gl.EPN <- gl.recalc.metrics(gl.EPN)
gl.EPN
# 2,235 genotypes,  114,462 SNPs , size: 450.5 Mb
# 

gl.EPR <- gl.keep.pop(gel, c("EPR"), as.pop="SPECIES_ID")
gl.EPR <- gl.filter.monomorphs(gl.EPR)
gl.EPR <- gl.recalc.metrics(gl.EPR)

gl.UNK <- gl.keep.pop(gel, c("UNK"), as.pop="SPECIES_ID")
gl.UNK <- gl.filter.monomorphs(gl.UNK)
gl.UNK <- gl.recalc.metrics(gl.UNK)




### Saving raw species files

gl.save(gl.EPN, file="../DATA_intermediate/04_merging_EPN.Rdata")
gl.save(gl.EPR, file="../DATA_intermediate/04_merging_EPR.Rdata")
gl.save(gl.UNK, file="../DATA_intermediate/04_merging_UNK.Rdata")
gl.UNK@other

#################
### FUNCTIONS ###
#################


# A function to obtain per sample call rate
ind_call_rate = function(row, output){
  x = length(row[!is.na(row)])/length(row)
  return(x)
}


# Function for calculating heterozygosity per sample
ind_freq_het = function(row, output){
  hets <- row[row == 2]
  hets_complete <- complete.cases(hets)
  hets_nonan <- hets[hets_complete]
  x = length(hets_nonan)/length(row[!is.na(row)])
  return(x)
}


# Function for calculating heterozygosity per sample
ind_freq_het = function(row, output){
  hets <- row[row == 1]
  hets_complete <- complete.cases(hets)
  hets_nonan <- hets[hets_complete]
  x = length(hets_nonan)/length(row[!is.na(row)])
  return(x)
}




# Cehcking sample call rate
df.test <- as.data.frame(gl.EPN)
df.test$call_rate <- apply(df.test, 1, ind_call_rate)
mean(df.test$call_rate)




############################
### Call rate in species ###
############################

ggplot(gl.EPN@other$loc.metrics) + aes(x = CallRate) +
  geom_histogram(bins=100, fill = "grey70", colour = "white") +
  geom_vline(xintercept = 0.5, colour = "skyblue") +
  #geom_vline(xintercept = 0.75, colour = "skyblue2") +
  #geom_vline(xintercept = 0.85, colour = "skyblue4") +
  theme(panel.background = element_rect(fill =NA, colour = "black"))

ggplot(gl.EPR@other$loc.metrics) + aes(x = CallRate) +
  geom_histogram(bins=100, fill = "grey70", colour = "white") +
  geom_vline(xintercept = 0.5, colour = "skyblue") +
  #geom_vline(xintercept = 0.75, colour = "skyblue2") +
  #geom_vline(xintercept = 0.85, colour = "skyblue4") +
  theme(panel.background = element_rect(fill =NA, colour = "black"))

ggplot(gl.UNK@other$loc.metrics) + aes(x = CallRate) +
  geom_histogram(bins=100, fill = "grey70", colour = "white") +
  geom_vline(xintercept = 0.5, colour = "skyblue") +
  #geom_vline(xintercept = 0.75, colour = "skyblue2") +
  #geom_vline(xintercept = 0.85, colour = "skyblue4") +
  theme(panel.background = element_rect(fill =NA, colour = "black"))



### Finding loci names with poor call rate within pop

lociPoorCallRate <- function(gel, minCR) {
  gel_meta <- gel@other$loc.metrics
  gel_meta$locus <- locNames(gel)
  gel_minCR <- gel_meta %>% filter(CallRate < minCR) %>% pull(locus)
  return(gel_minCR)
}

poor.loci.EPN <- lociPoorCallRate(gl.EPN, 0.5)
poor.loci.EPR <- lociPoorCallRate(gl.EPR, 0.5)
poor.loci.UNK <- lociPoorCallRate(gl.UNK, 0.5)
length(poor.loci.EPN)
length(poor.loci.EPR)
length(poor.loci.UNK)

poor.loci.TOT <- unique(c(poor.loci.EPN, poor.loci.EPR, poor.loci.UNK))
length(poor.loci.TOT)





### Finding loci names with high heterozygosity within pop

head(gel@other$loc.metrics)

lociHighHet <- function(gel, maxHet) {
  gel_meta <- gel@other$loc.metrics
  gel_meta$locus <- locNames(gel)
  gel_maxHet <- gel_meta %>% filter(FreqHets > maxHet) %>% pull(locus)
  return(gel_maxHet)
}

highHet.loci.EPN <- lociHighHet(gl.EPN, 0.5)
highHet.loci.EPR <- lociHighHet(gl.EPR, 0.5)
highHet.loci.UNK <- lociHighHet(gl.UNK, 0.5)
length(highHet.loci.EPN)
length(highHet.loci.EPR)
length(highHet.loci.UNK)


ggplot(gl.EPR@other$loc.metrics) + aes(x = FreqHets) +
  geom_histogram(bins=100, fill = "grey70", colour = "white") +
  labs(y = "N SNPs") +
  geom_vline(xintercept = 0.5, colour = "skyblue") +
  theme(panel.background = element_rect(fill =NA, colour = "black"))

highHet.loci.TOT <- unique(c(highHet.loci.EPN, highHet.loci.EPR, highHet.loci.UNK))
length(highHet.loci.TOT)



### Finding loci names with low MAF (<0.01) within pop

lociLowMAF <- function(gel, minMAF) {
  gel_meta <- gel@other$loc.metrics
  gel_meta$locus <- locNames(gel)
  gel_minHet <- gel_meta %>% filter(maf < minMAF) %>% pull(locus)
  return(gel_minHet)
}

lowMAF.loci.EPN <- lociLowMAF(gl.EPN, 0.01)
lowMAF.loci.EPR <- lociLowMAF(gl.EPR, 0.01)
lowMAF.loci.UNK <- lociLowMAF(gl.UNK, 0.01)
length(lowMAF.loci.EPN)
length(lowMAF.loci.EPR)
length(lowMAF.loci.UNK)

ggplot(gl.UNK@other$loc.metrics) + aes(x = maf) +
  geom_histogram(bins=100, fill = "grey70", colour = "white") +
  labs(y = "N SNPs") +
  geom_vline(xintercept = 0.01, colour = "skyblue") +
  theme(panel.background = element_rect(fill =NA, colour = "black"))

lowMAF.loci.TOT <- unique(c(lowMAF.loci.EPN, lowMAF.loci.EPR, lowMAF.loci.UNK))
length(lowMAF.loci.TOT)


### Filtering bad SNPs

loci.to.drop <- unique(c(poor.loci.TOT, highHet.loci.TOT, lowMAF.loci.TOT))
length(loci.to.drop)


gel.CR <- gl.drop.loc(x = gel, loci.to.drop)
gel.CR <- gl.filter.monomorphs(gel.CR)
gel.CR <- gl.recalc.metrics(gel.CR)
gel.CR

# 2,322 genotypes,  36,712 SNPs



#######################################
# Saving a dataset before MAF filtering

loci.to.drop.2 <- unique(c(poor.loci.TOT, highHet.loci.TOT))
length(loci.to.drop.2)
gel
gel.CR.2 <- gl.drop.loc(x = gel, loci.to.drop.2)
gel.CR.2 <- gl.filter.monomorphs(gel.CR.2)
gel.CR.2 <- gl.recalc.metrics(gel.CR.2)
gel.CR.2

gl.save(gel.CR.2, file="../DATA_intermediate/06_filtering_NoMAF.Rdata")




#################################
### Checking sample call rate ###
#################################


df.tot <- as.data.frame(gel.CR)
df.tot$call_rate <- apply(df.tot, 1, ind_call_rate)
df.tot_CR <- df.tot %>% dplyr::select(call_rate)
df.tot_CR[is.na(df.tot_CR$call_rate),]
head(df.tot_CR)

p1 <- ggplot(df.tot_CR) + aes(x = call_rate) +
  geom_histogram(fill="grey20",color="white") +
  ggtitle("Sample CALL RATE - all species") +
  labs(x = "Call rate", y = "N samples") +
  #scale_x_continuous(limits = c(0,1)) +
  theme(panel.background = element_rect(fill=NA, color="black"),
        panel.grid = element_blank())
p1

df.tot_CR %>% filter(call_rate < 0.2) %>% dim()

png("06_filtering_TOT_sample_call_rate.png",w=800,h=700,res=150)
p1
dev.off()




##############################
### Checking SNP call rate ###
##############################

# 3 different thresholds 0.5, 0.75, and 0.85

p2 <- ggplot(gel.CR@other$loc.metrics) + aes(x = CallRate) +
  geom_histogram(bins=100, fill = "grey70", colour = "white") +
  geom_vline(xintercept = 0.5, colour = "skyblue") +
  geom_vline(xintercept = 0.75, colour = "skyblue2") +
  geom_vline(xintercept = 0.85, colour = "skyblue4") +
  theme(panel.background = element_rect(fill =NA, colour = "black"))
p2

png("06_filtering_TOT_locus_call_rate.png",w=800,h=700,res=150)
p2
dev.off()


######################
### Checking stats ###
######################

head(gel.CR@other$loc.metrics)

loci.gl <- gel.CR@other$loc.metrics %>% dplyr::select(CallRate, OneRatioRef, OneRatioSnp, FreqHomRef, FreqHomSnp, 
                                                                  FreqHets, PICRef, PICSnp, AvgPIC, maf
) %>% gather(key = "Stat", value = "Value")

p3 <- ggplot(loci.gl) + aes(x = Value) +
  geom_histogram(fill="grey80",color="white") +
  facet_wrap(~Stat, scales = "free", ncol = 3) +
  ggtitle("Per locus statistics - all species") +
  labs(x = "Value", y = "N loci") +
  theme(panel.background = element_rect(fill=NA, color="black"),
        panel.grid = element_blank())
p3

png("06_filtering_TOT_locus_stats.png", w=1500, h=2000, res=150)
p3
dev.off()





####################################
### Overall SNP call rate filter ###
####################################

# 50%
gel.CR.CR <- gl.filter.callrate(x = gel.CR,
                                             method = "loc",
                                             threshold = 0.5,
                                             mono.rm = TRUE,
                                             recalc = TRUE,
                                             recursive = TRUE)
gel.CR.CR
# 2,322 genotypes,  36,430 SNPs






######################
### No secondaries ###
######################

#gl.report.secondaries(gel.CR.CR)
#gel.CR.CR.nosec <- gl.filter.secondaries(gel.CR.CR, method = "best", verbose=3)
#gel.CR.CR.nosec <- gl.filter.monomorphs(gel.CR.CR.nosec)
#gel.CR.CR.nosec <- gl.recalc.metrics(gel.CR.CR.nosec)

removeSecondaries <- function(gel) {
  df <- data.frame("locus" = locNames(gel))
  df$CallRate <- gel@other$loc.metrics$CallRate
  df[c('CloneID', 'TagPos', 'Alleles')] <- str_split_fixed(df$locus, '-', 3)
  df_grouped <- df %>% group_by(CloneID) %>% mutate(n = n()) %>% ungroup()
  high.CR.loci <- df_grouped %>% group_by(CloneID) %>% filter(CallRate == max(CallRate)) %>% pull(locus)

  gel.nosec <- gl.keep.loc(x = gel, high.CR.loci)
  gel.nosec <- gl.filter.monomorphs(gel.nosec)
  gel.nosec <- gl.recalc.metrics(gel.nosec)
  return(gel.nosec)
}

gel.CR.CR.nosec <- removeSecondaries(gel.CR.CR)


# Number of loci before
dim(gel.CR.CR)[2]
# 36430
# Number of loci after filtering
dim(gel.CR.CR.nosec)[2]
# 31294






##############
### Saving ###
##############

gl.save(gel.CR.CR.nosec, file="../DATA_intermediate/06_filtering_TOT_CR50_nosec.Rdata")


#gel.CR.CR.nosec <- gl.load("../DATA_intermediate/06_filtering_TOT_CR50_nosec.Rdata")
#gel.CR.CR.nosec
#gl.write.csv(gel.CR.CR.nosec, outfile = "../DATA_intermediate/06_filtering_TOT_CR50_nosec.csv", outpath = tempdir(), verbose = NULL)





################
### Metadata ###
################

# Locus metrics
dmeta <- gel.CR.CR.nosec@other$loc.metrics
head(dmeta)
write.table(dmeta, "06_filtering_TOT_locus_metrics.tsv", sep="\t", row.names = TRUE, col.names = TRUE, quote=FALSE, append=FALSE)

# Sample metrics
dmeta_sample <- gel.CR.CR.nosec@other$ind.metrics
head(dmeta_sample)
df.gel.CR.CR.nosec <- as.data.frame(gel.CR.CR.nosec)
df.gel.CR.CR.nosec$call_rate <- apply(df.gel.CR.CR.nosec, 1, ind_call_rate)
df.gel.CR.CR.nosec$id <- row.names(df.gel.CR.CR.nosec)
df.gel.CR.CR.nosec <- df.gel.CR.CR.nosec %>% dplyr::select(id, call_rate)
head(df.gel.CR.CR.nosec)
dmeta_merged <- merge(dmeta_sample, df.gel.CR.CR.nosec, by = "id", sort=FALSE)
dim(dmeta_merged)
head(dmeta_merged)
write.table(dmeta_merged, "06_filtering_TOT_sample_metrics.tsv", sep="\t", row.names = TRUE, col.names = TRUE, quote=FALSE, append=FALSE)



##################################################################


# Function for calculating heterozygosity per sample
ind_freq_het = function(row, output){
  hets <- row[row == 1]
  hets_complete <- complete.cases(hets)
  hets_nonan <- hets[hets_complete]
  x = length(hets_nonan)/length(row[!is.na(row)])
  return(x)
}


# A function to obtain per sample call rate
ind_call_rate = function(row, output){
  x = length(row[!is.na(row)])/length(row)
  return(x)
}



# Checking heterozygotes
df.test <- as.data.frame(gel)
df.test[1:10, 1:10]
length(gel@ind.names)
row.names(df.test)

df.meta <- gel@other$ind.metrics
df.meta$het <- apply(df.test, 1, ind_freq_het)
df.meta$call_rate <- apply(df.test, 1, ind_call_rate)
head(df.meta)

p1 <- ggplot(df.meta) + 
  geom_point(aes(x = call_rate, y = het, color = SPECIES_ID))
p1

png("06_filtering_TOT_heterozygosity.png", w = 800, h=600, res=150)
p1
dev.off()


