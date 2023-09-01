library(GenomicRanges)


load('../mutation_data/hg19.1Mb.ranges.Polak.Nature2015.RData')
waddell = readRDS("../mutation_data/mutation_x_patient_grange.rds")
count_overlaps = data.frame(lapply(waddell, 
                            function(x) countOverlaps(interval.ranges, x)))
write.csv(count_overlaps, "../mutation_data/waddell_per_patient.csv")