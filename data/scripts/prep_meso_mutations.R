library(tidyverse)
library(readxl)
library(GenomicRanges)

load('../mutation_data/hg19.1Mb.ranges.Polak.Nature2015.RData')
rownames = names(interval.ranges)
biphasic_mesomics = read.csv("../mutation_data/Mesothelioma_MMB_MutCountperBin.txt",
                             sep="\t")
sarcomatoid_mesomics = read.csv("../mutation_data/Mesothelioma_MMS_MutCountperBin.txt",
                             sep="\t")
epithelioid_mesomics = read.csv("../mutation_data/Mesothelioma_MME_MutCountperBin.txt",
                                sep="\t")

aggregate_mutations <- function(mut_df, name) {
  agg = rowSums(mut_df[, 4:ncol(mut_df)])
  agg = data.frame(cancer_type = agg)
  colnames(agg) = name
  rownames(agg) = rownames
  return(agg)
}

agg_epithelioid = aggregate_mutations(epithelioid_mesomics, "epithelioid_mesomics")
agg_sarcomatoid = aggregate_mutations(sarcomatoid_mesomics, "sarcomatoid_mesomics")
agg_biphasic = aggregate_mutations(biphasic_mesomics, "biphasic_mesomics")

waddell = read.csv("../mutation_data/mesothelioma_WGS_Waddell.csv")

# mut_df = mut_df[, 1:3]

# mut_df["agg"] = agg
# 
# mut_overlaps = as.data.frame(countOverlaps(interval.ranges, mut_granges))



