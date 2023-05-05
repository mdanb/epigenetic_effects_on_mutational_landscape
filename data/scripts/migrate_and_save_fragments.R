library(parallel)
library(rtracklayer)
library(optparse)
source("count_overlaps_utils.R")

option_list <- list(
  make_option("--dataset", type="character"),
  make_option("--cores", type="integer"),
  make_option("--tissue", type="character", default="all")
)

args = parse_args(OptionParser(option_list=option_list))

dataset = args$dataset
cores = args$cores
tissue = args$tissue
ch = import.chain("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/hg38ToHg19.over.chain")

helper <- function(files, migrated_filepaths, ch, cores) {
  # files = lapply(strsplit(files, split = "/"), "[", 9)
  print("IMPORTING")
  print(files)
  fragments = mclapply(files, import, format="bed", mc.cores=cores)
  
  print("MIGRATING")
  print(files)
  migrated_fragments = mclapply(fragments, migrate_bed_file_to_hg37, ch, 
                                mc.cores=cores)
  print("EXPORTING")
  print(files)
  mclapply(seq_along(migrated_fragments), function(i) {export(migrated_fragments[[i]],
                                                           migrated_filepaths[i],
                                                           format="bed",
                                                           index=T)})
}

get_files_not_done <- function(files, dir_path) {
  # migrated_filepaths_bgz = paste(dir_path, gsub(".gz", ".bgz", 
  #                                               lapply(strsplit(files, split = "/"), "[", 9)), 
  #                                sep = "/")
  migrated_filepaths_tbi = paste(files, "tbi", sep=".")
  migrated_filepaths = paste(dir_path, 
                             lapply(strsplit(files, split = "/"), "[", 9),
                             sep = "/")
  files_not_done = c()
  migrated_filepaths_not_done = c()
  for (i in 1:length(migrated_filepaths)) {
    bgz_file = files[i]
    tbi_file = migrated_filepaths_tbi[i]
    if (!file.exists(bgz_file) || !file.exists(tbi_file)) {
      migrated_filepaths_not_done = append(migrated_filepaths_not_done, 
                                           migrated_filepaths[i])
      files_not_done = append(files_not_done, files[i])
    }
  }
  return(list(files_not_done, migrated_filepaths_not_done))
}

if (dataset == "Bingren") {
  dir_path = "/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/bed_files/bingren_scATAC/migrated_to_hg19"
  dir.create(dir_path)
  files = list.files("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/bed_files/bingren_scATAC",
                     pattern="bed.bgz",
                     full.names=TRUE)
} else if (dataset == "Tsankov") {
  dir_path = "/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/bed_files/Tsankov_scATAC/migrated_to_hg19"
  dir.create(dir_path)
  files = list.files("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/bed_files/Tsankov_scATAC",
                     pattern="tsv.gz",
                     full.names=TRUE)
} else if (dataset == "Greenleaf_brain") {
  dir_path = "/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/bed_files/greenleaf_brain_scATAC/migrated_to_hg19"
  dir.create(dir_path)
  files = list.files("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/bed_files/greenleaf_brain_scATAC",
                     pattern="tsv.gz",
                     full.names=TRUE)
} else if (dataset == "Greenleaf_pbmc_bm") {
  dir_path = "/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/bed_files/greenleaf_pbmc_bm_scATAC/migrated_to_hg19"
  dir.create(dir_path)
  files = list.files("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/bed_files/greenleaf_pbmc_bm_scATAC",
                     pattern="tsv.gz",
                     full.names=TRUE)
}

if (tissue == "all") {
   tissue = "*"
}
files = files[grepl(tissue, files)]
filepaths = get_files_not_done(files, dir_path)
files = filepaths[[1]]
migrated_filepaths = filepaths[[2]]
split_files = split(files, ceiling(seq_along(files)/cores))
split_migrated_filepaths = split(migrated_filepaths, 
                                 ceiling(seq_along(files)/cores))
if (!is.null(files)) {
  for (i in 1:length(split_files)) {
    helper(split_files[[i]], split_migrated_filepaths[[i]], ch, cores)
  }
}
