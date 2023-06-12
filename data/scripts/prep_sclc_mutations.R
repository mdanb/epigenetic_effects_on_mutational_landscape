library(tidyverse)
library(readxl)
library(GenomicRanges)

load('/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/mutation_data/hg19.1Mb.ranges.Polak.Nature2015.RData')
rownames = names(interval.ranges)
mut_df = read.csv("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/mutation_data/SCLC_MutCountperBin.txt",
                  sep="\t")
agg = rowSums(mut_df[, 4:ncol(mut_df)])
agg = data.frame(SCLC = agg)
rownames(agg) = rownames
# mut_df = mut_df[, 1:3]

# mut_df["agg"] = agg
# 
# mut_overlaps = as.data.frame(countOverlaps(interval.ranges, mut_granges))
write.csv(agg, paste(
                "/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/processed_data",
                "sclc_count_overlaps.csv", sep="/"))
