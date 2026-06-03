rm(list=ls())
setwd("C:/Users/aniaf/Projects/BlackSpruce/929_heterozygosity")

library(ggplot2)
library(tidyr)
library(dplyr)
library(dartR)
library(stringr)


#===================#
#       DATA        #
#===================#


### Genotypes: no MAF filter, snp call rate about 50%, and no high heterozygosity snps (>0.5), no filter for secondaries

# Open dartR dataset
gel <- gl.load("../DATA_intermediate/06_filtering_NoMAF.Rdata")
class(gel)
#samples <- gel@other$ind.metrics$id
#length(samples)
gel
dim(gel)


# Select SNPs with a call rate > 50%
gel.flt <- gl.filter.callrate(x = gel,
                                method = "loc",
                                threshold = 0.5,
                                mono.rm = TRUE,
                                recalc = TRUE,
                                recursive = TRUE)
gel.flt

# Removing secondary SNPs
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

gel.flt.nosec <- removeSecondaries(gel.flt)
gel.flt.nosec



#===================#
#      Checks       #
#===================#

######## checking if "1" corresponds to heterozygotes

gel.flt.nosec@gen
gel.flt.nosec@other$loc.metrics %>% head()
gel.flt.nosec@other$loc.metrics$FreqHets %>% length()
gel.flt.nosec

tab <- as.data.frame(as.matrix(gel.flt.nosec))
class(tab)
tab[c(1:10), c(1:10)]


cols_1 <- colSums(tab == 1, na.rm = TRUE)
df_cols_1 <- as.data.frame(cols_1)
colnames(df_cols_1) <- c("hets_1")

cols_2 <- colSums(tab == 2, na.rm = TRUE)
df_cols_2 <- as.data.frame(cols_2)
colnames(df_cols_2) <- c("hets_2")

colna_nona <- colSums(!is.na(tab))
df_cols_nona <- as.data.frame(colna_nona)
colnames(df_cols_nona) <- c("nonans")

df_cols <- cbind(df_cols_1, df_cols_2, df_cols_nona)
dim(df_cols)
head(df_cols)
df_cols$freqHets1 <- df_cols$hets_1/df_cols$nonans
df_cols$freqHets2 <- df_cols$hets_2/df_cols$nonans
df_cols$FreqHets <- gel.flt.nosec@other$loc.metrics$FreqHets

head(df_cols)
plot(df_cols$FreqHets, df_cols$freqHets1)

###############################################


#====================#
#   Filter samples   #
#====================#


# Selecting individuals - clean set of EPN
dmeta <- read.csv("../DATA_intermediate/23_filter_EPN_indiv_metrics.tsv", sep="\t")
dim(dmeta)
head(dmeta)
unique(dmeta$SPECIES_ID)

# Selecting pops with > 3 ind (7 or more in reality)
dpops <- dmeta %>% group_by(POP) %>% dplyr::summarise(n=n())
dpops <- dpops %>% arrange(n)
dpops
dpops_sub <- dpops %>% filter(n>=5)
hist(dpops_sub$n, breaks=20)
mean(dpops_sub$n)
summary(dpops_sub$n)
samples_flt <- dmeta %>% filter(POP %in% dpops_sub$POP) %>% pull(id)
length(samples_flt)
fmeta <- dmeta %>% filter(POP %in% dpops_sub$POP)
dim(fmeta)
head(fmeta)


# Filtering genotypes
gel.clean <- gl.keep.ind(gel.flt.nosec, ind.list = samples_flt, recalc = TRUE, mono.rm = TRUE)
gel.clean



#=======================#
#   Creating genlight   #
#=======================#


# Genotypes and metadata
gts <- data.frame(gel.clean)
gts[c(1:5), c(1:5)]
dmeta.clean <- gel.clean@other$loc.metrics
imeta.clean <- gel.clean@other$ind.metrics
imeta.clean <- merge(imeta.clean, dmeta[c("id","cluster","group")], by = "id", sort=F, all.x=T)
imeta.clean %>% head()
dmeta.clean %>% head()
dmeta.clean$locus <- row.names(dmeta.clean)
dmeta.clean[c('CloneID', 'TagPos', 'Alleles')] <- str_split_fixed(dmeta.clean$locus, '-', 3)
head(dmeta.clean)

strata_list <- list('region' = imeta.clean$STRATA,
                    'group' = imeta.clean$group,
                    'cluster' = imeta.clean$cluster)
strata_list

genl <- new("genlight", gts, ploidy=2,
            ind.names=imeta.clean$id,
            loc.names=dmeta.clean$locus,
            loc.all=dmeta.clean$Alleles,
            chromosome=dmeta.clean$CloneID,
            position=dmeta.clean$TagPos,
            pop=imeta.clean$POP,
            other=strata_list)

genl
genl@other


save(genl, file="../DATA_intermediate/06_filtering_NoMAF_genlight.Rdata")

genl@pop
genl@ind.names




#=======================#
#     Creating csv      #
#=======================#

# SNP matrix
genl@gen[[1]]
mat <- as.data.frame(as.matrix(genl))
mat[1:10,1:10]

# POP and id
df.pops <- data.frame("id" = genl@ind.names, "POP" = genl@pop)
row.names(df.pops) <- df.pops$id
head(df.pops)

# Combine
dout <- cbind(df.pops, mat)
dout[1:10,1:10]
dim(dout)

# Save to csv
write.csv(dout, file = "01_prep_inputs.csv", row.names = F)




