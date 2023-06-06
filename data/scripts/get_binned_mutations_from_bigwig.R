library(rtracklayer)
library(tibble)
library(dplyr)
library(optparse)
library(parallel)

load('../mutation_data/hg19.1Mb.ranges.Polak.Nature2015.RData')
chr_ranges = read.csv("../processed_data/chr_ranges.csv")
dir.create("../processed_data/binned_atac")

option_list <- list( 
  make_option("--cohort", type="character"),
  make_option("--cores", type="integer")
)

create_binned_atac <- function(file, interval_ranges, chr_ranges) {
  filename = unlist(strsplit(file, split="/"))[4]
  filename = gsub("\\.insertions\\.bw", "", filename)
  filename = paste(filename, "rds", sep=".")
  file_path = paste(dirpath, filename, sep="/")
  
  if (!file.exists(file_path)) {
    print(paste0("Importing ", file, "..."))
    bigwig = import(file)
    print("Done!")
    
    print("Computing overlaps...")
    overlaps = as_tibble(findOverlaps(interval.ranges, bigwig))
    print("Done!")
    
    score = bigwig$score
    overlaps["hit_score"] = score[overlaps$subjectHits]
    overlaps = overlaps %>% 
      group_by(queryHits) %>%
      summarise(sum(hit_score))
    colnames(overlaps) = c("bin", "counts")
    binned_overlaps = data.frame(counts = rep(0, 2998))
    rownames(binned_overlaps) = chr_ranges[["x"]]
    binned_overlaps[overlaps[["bin"]], "counts"] = overlaps[["counts"]]
    dirpath = "../processed_data/binned_atac"
    saveRDS(binned_overlaps, file_path)
  }
}

args = parse_args(OptionParser(option_list=option_list))
cores = args$cores
cohort = args$cohort
files = list.files(paste0("../bigwig/", cohort), full.names = T)

  
mclapply(files, create_binned_atac,
         interval_ranges=interval.ranges,
         chr_ranges=chr_ranges,
         mc.cores=cores)


