rm(list=ls())
setwd("C:/Users/aniaf/Projects/BlackSpruce/928_amova")

library(ggplot2)
library(tidyr)
library(dplyr)
library(dartR)
library(poppr)
library(stringr)



#===================#
#       DATA        #
#===================#

# Open genlight
load("../DATA_intermediate/23_filter_EPN_genlight.Rdata")
class(genl)
genl


#===================#
#     Data prep     #
#===================#

### Subset populations with more than 5 trees
dmeta <- read.csv("../DATA_intermediate/23_filter_EPN_indiv_metrics.tsv", sep="\t")
head(dmeta)
dim(dmeta)
unique(dmeta$SPECIES_ID)

dpops <- dmeta %>% group_by(POP) %>% dplyr::summarise(n=n())
dpops <- dpops %>% arrange(n)
dpops_sub <- dpops %>% filter(n>=5)
samples_flt <- dmeta %>% filter(POP %in% dpops_sub$POP) %>% pull(id)
fmeta <- dmeta %>% filter(POP %in% dpops_sub$POP)
dim(fmeta)
head(fmeta)

# Subsetting samples from genotype matrix
ngel <- genl[indNames(genl) %in% samples_flt]
ngel
nPop(ngel)
dim(ngel)



#===============================#
#  Setting hierarchical groups  #
#===============================#

# ~cluster/pop

# Adding pop names to strata structure
new.strata <- data.frame(other(ngel))
head(new.strata)
new.strata$pop <- pop(ngel)
new.strata <- new.strata %>% dplyr::select(cluster, pop)
head(new.strata)

# Setting new strata
strata(ngel) <- new.strata
nameStrata(ngel) <- ~cluster/pop
head(strata(ngel, ~cluster/pop))




#===================#
#       AMOVA       #
#===================#

table(strata(ngel, ~pop))  # Populations
table(strata(ngel, ~cluster))

amova.group <- poppr.amova(ngel, ~cluster/pop)
amova.group

save(amova.group, file="02_amova_cluster.Rdata")


# Testing results
# Import if necessary
load("02_amova_cluster.Rdata")

amova.test <- randtest(amova.group, nrepet = 999)
amova.test
save(amova.test, file="02_amova_cluster_test.Rdata")

# Saving output as a table
dres <- as.data.frame(amova.group$results)
head(dres)
dvar <- as.data.frame(amova.group$componentsofcovariance)
head(dvar)
dres$Sigma <- dvar$Sigma
dres$Var_perc <- dvar$`%`
head(dres)

write.table(dres, "02_amova_cluster_result.tsv", sep="\t",col.names=T, row.names=T, quote=F, append=F)


