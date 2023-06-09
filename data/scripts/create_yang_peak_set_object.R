library(GenomicRanges)
source("count_overlaps_utils.R")
library(rtracklayer)                                                            

df = read.csv("../pan_cancer_peak_set.txt", sep="\t")
interval.ranges = GRanges(seqnames=df[["seqnames"]],
                           ranges = IRanges(
                                        start=df[["start"]],
                                        end=df[["end"]]
                                        )
                          )
save(interval.ranges, file="../peak_set_yang.Rdata")
