library(tidyverse)
library(readxl)
library(GenomicRanges)

load('/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/mutation_data/hg19.1Mb.ranges.Polak.Nature2015.RData')

mut_df = read_excel("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/mutation_data/sclc.xlsx",
                sheet=3, skip=1, na=c("", "NA"))

mut_df = mut_df[, 1:23]

# remove weird mutation
mut_df = mut_df[-which(mut_df["Start"] > mut_df["End"]), ]
mut_granges = GRanges(mut_df[["Chromosome"]], 
                      ranges=IRanges(start=mut_df[["Start"]], 
                      end=mut_df[["End"]]))

mut_overlaps = as.data.frame(countOverlaps(interval.ranges, mut_granges))
colnames(mut_overlaps) = "SCLC"
write.csv(as.data.frame(mut_overlaps), paste(
                "/ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/data/processed_data",
                "sclc_count_overlaps.csv", sep="/"))
