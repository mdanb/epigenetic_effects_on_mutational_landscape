library(GenomicRanges)
library(dplyr)

df = readRDS("../processed_data/MSI_SNP_consensus_hg19.RDS")
ir = GRanges(ranges=IRanges(start=df[["start"]], end=df[["end"]]), 
             seqnames=df[["seqnames"]])
load("../mutation_data/hg19.1Mb.ranges.Polak.Nature2015.RData")
msi_high = data.frame(msi_high = countOverlaps(interval.ranges, ir))
msi_high["X"] = rownames(msi_high)
colorect_adeno = read.csv("../processed_data/per_patient_mutations/ColoRect-AdenoCA_per_donor.csv")
extra_donor = colorect_adeno[c("X", "DO44090")]
msi_high = left_join(msi_high, extra_donor)
msi_high["msi_high"] = msi_high["msi_high"] + msi_high["DO44090"]
msi_high = msi_high[c("X", "msi_high")]
mss = cbind(colorect_adeno["X"], 
            data.frame(mss=rowSums(colorect_adeno[, !grepl("X|DO44090",
                                                           colnames(colorect_adeno))])))
msi_high = left_join(msi_high, mss)

rownames(msi_high) = msi_high[["X"]]
msi_high = msi_high[c("msi_high", "mss")]

write.csv(msi_high, "../processed_data/msi_high.csv")