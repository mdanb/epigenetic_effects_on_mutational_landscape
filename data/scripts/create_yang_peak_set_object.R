library(GenomicRanges)

df = read.csv("../pan_cancer_peak_set.txt", sep="\t")
interval.ranges = GRanges(seqnames=df[["seqnames"]],
                           ranges = IRanges(
                                        start=df[["start"]],
                                        end=df[["end"]]
                                        )
                          )
save(interval.ranges, file="../peak_set_yang.Rdata")