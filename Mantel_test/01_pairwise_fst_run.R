rm(list=ls())
setwd("C:/Users/aniaf/Projects/BlackSpruce/925_Mantel")


library(tidyr)
library(dplyr)
library(SNPRelate)
library(assigner) # calculating pairwise WC FST



#######################
#------- DATA---------#
#######################

load("01_pairwise_fst_genlight.Rdata")

### Saving pop order
write.table(data.frame("pop"=levels(gl_random_filtered$pop)),
            "01_pairwise_fst_pop_order.tsv", sep="\t", col.names=T,row.names=F,append=F, quote=F)


### Converting genlight to tidy
gel.tidy <- radiator::tidy_genlight(gl_random_filtered)
gel.tidy$GT_BIN
gel.tidy$GT <- gel.tidy$GT_BIN
gel.tidy
dim(gel.tidy)


###############################################################################################
## Calculating pairwise Fst

### All data
pairwise.fst.wc <- fst_WC84(gel.tidy, snprelate = FALSE, pairwise = TRUE,
                            ci = TRUE, parallel.core = 7,
                            heatmap.fst = FALSE, 
                            strata = NULL)

save(pairwise.fst.wc, file = "01_pairwise_fst_run.RData")

