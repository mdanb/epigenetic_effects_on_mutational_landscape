library(optparse)
library(rtracklayer)

load('../mutation_data/hg19.1Mb.ranges.Polak.Nature2015.RData')

option_list <- list( 
  make_option("--cancer_type", type="character"),
  make_option("--subsampled", action="store_true", default=FALSE)
)

args = parse_args(OptionParser(option_list=option_list))

cancer_type = args$cancer_type


if (subsampled) {
  for (dir in list.dirs(paste("../mutation_data/bed_files", paste(cancer_type, 
                                                                  "subsampled", 
                                                                  sep="_"), sep="/"),
                        recursive=F)) {
    mutations_df = data.frame(row.names = seq(1,2128))
    for (i in seq(0, 99)) {
      sample = paste0("S", i, ".bed")
      counts = import(paste(dir, "IntersectedCount_paz_Cancergroup",
                            paste("IntersectedCount", cancer_type, sample,
                                  sep="_"), sep="/"))
      # seqlevelsStyle(counts) <- "UCSC" # This will prefix "chr" to sequence names if not already present
      # overlaps = findOverlaps(interval.ranges, counts)
      # names(counts)[subjectHits(overlaps)] = names(interval.ranges)[queryHits(overlaps)]
      mutations_df = cbind(mutations_df, counts$score)
      seed = floor(i / 10) + 1
      fold = i %% 10 + 1
      colnames(mutations_df)[ncol(mutations_df)] = paste(cancer_type, "n", basename(dir),
                                                         "seed", seed, "fold", fold,
                                     sep="_")
      
    }
    seqlevelsStyle(counts) <- "UCSC" # This will prefix "chr" to sequence names if not already present
    overlaps = findOverlaps(interval.ranges, counts)
    names(counts)[subjectHits(overlaps)] = names(interval.ranges)[queryHits(overlaps)]
    rownames(mutations_df) = names(counts)
    fn = paste0(paste(cancer_type, "n", basename(dir), sep="_"), ".csv")
    write.csv(mutations_df, paste("../processed_data", fn, sep="/"))
  }
} else {
  mutations_df = data.frame(row.names = seq(1,2128))
  counts = import(paste(dir, "IntersectedCount_paz_Cancergroup",
                        paste("IntersectedCount", cancer_type, sep="_"), 
                        sep="/"))
  mutations_df = cbind(mutations_df, counts$score)
  colnames(mutations_df)[ncol(mutations_df)] = cancer_type
  seqlevelsStyle(counts) <- "UCSC" # This will prefix "chr" to sequence names if not already present
  overlaps = findOverlaps(interval.ranges, counts)
  names(counts)[subjectHits(overlaps)] = names(interval.ranges)[queryHits(overlaps)]
  rownames(mutations_df) = names(counts)
  fn = paste0(cancer_type, ".csv")
  write.csv(mutations_df, paste("../processed_data", fn, sep="/"))
}