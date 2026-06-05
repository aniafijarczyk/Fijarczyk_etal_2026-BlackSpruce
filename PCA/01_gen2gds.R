rm(list=ls())
setwd("C:/Users/aniaf/Projects/BlackSpruce/924_PCA")

library(tidyr)
library(dplyr)
library(dartR)
library(stringr)



#====================#
#-------DATA---------#
#====================#


#===================#
#     Function      #
#===================#


convert_gen2gds <- function(input_file, output_file) {
  
  # Importing data
  gel <- gl.load(input_file)
  
  # Metadata
  dmeta <- as.data.frame(gel@other$ind.metrics)
  
  # Setting chromosome names 
  df.loci <- gel$other$loc.metrics
  df.loci$locus <- locNames(gel)
  df.loci[c('CloneID', 'TagPos', 'Alleles')] <- str_split_fixed(df.loci$locus, '-', 3)
  
  # Chromosomes
  chroms <- df.loci$CloneID
  gel$chromosome <- as.factor(chroms)
  
  # Positions
  positions <- df.loci$TagPos
  gel$position <- as.factor(positions)
  
  # Populations
  pops <- dmeta$POP
  npops <- replace(pops, is.na(pops), 'popX')
  gel$pop <- as.factor(npops)
  
  # Loci metrics
  gel$other$loc.metrics <- df.loci
  
  gds <-  gl2gds(gel, outfile = output_file,
                 outpath = "./",
                 snp_pos = "TagPos",
                 snp_chr = "CloneID",
                 chr_format = "character",
                 verbose = NULL
  )
  
}


#========================#
#     Saving as gds      #
#========================#


### Dataset including all EPN and EPR
convert_gen2gds("../DATA_intermediate/23_filter_clean.Rdata", "23_filter_clean.gds")



### Dataset including only EPN
convert_gen2gds("../DATA_intermediate/23_filter_EPN.Rdata", "23_filter_EPN.gds")






