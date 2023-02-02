library(parallel)
library(rtracklayer)
library(optparse)
source("count_overlaps_utils.R")

option_list <- list(
  make_option("--dataset", type="character"),
  make_option("--cores", type="integer")
)

args = parse_args(OptionParser(option_list=option_list))

dataset = args$dataset
cores = args$cores
ch = import.chain("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/hg38ToHg19.over.chain")

helper <- function(files, dir_path, ch, cores) {
  fragments = mclapply(files, import, format="bed", mc.cores=cores)
  migrated_fragments = mclapply(fragments, migrate_bed_file_to_hg37, ch, 
                                mc.cores=cores)
  files = lapply(strsplit(files, split = "/"), "[", 9)
  migrated_filepaths = paste(dir_path, gsub(".gz", "", files),
                             sep = "/")
  
  mclapply(seq_along(migrated_fragments), function(i) {export(migrated_fragments[[i]],
                                                           migrated_filepaths[i],
                                                           format="bed",
                                                           index=T)})
}

if (dataset == "bing_ren") {
  dir_path = "/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/bed_files/bingren_scATAC/migrated_to_hg19"
  dir.create(dir_path)
  files = list.files("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/bed_files/bingren_scATAC",
                     pattern="bed.gz",
                     full.names=TRUE)
  for (l in split(files, ceiling(seq_along(files)/8))) {
    helper(l, dir_path, ch, cores)
  }
  
} else if (dataset == "tsankov") {
  dir_path = "/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/bed_files/Tsankov_scATAC/migrated_to_hg19"
  dir.create(dir_path)
  files = list.files("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/bed_files/Tsankov_scATAC",
                     pattern="tsv.gz",
                     full.names=TRUE)
  for (l in split(files, ceiling(seq_along(files)/8))) {
    helper(l, dir_path, ch, cores)
  }
  # fragments = mclapply(files, import, format="bed", mc.cores=cores)
  # migrated_fragments = mclapply(fragments, migrate_bed_file_to_hg37, ch, 
  #                               mc.cores=cores)
} else if (dataset == "greenleaf_brain") {
  dir_path = "/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/bed_files/greenleaf_brain_scATAC/migrated_to_hg19"
  dir.create(dir_path)
  files = list.files("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/bed_files/greenleaf_brain_scATAC",
                     pattern="tsv.gz",
                     full.names=TRUE)
  for (l in split(files, ceiling(seq_along(files)/8))) {
    helper(l, dir_path, ch, cores)
  }
  # fragments = mclapply(files, import, format="bed", mc.cores=cores)
  # migrated_fragments = mclapply(fragments, migrate_bed_file_to_hg37, ch, 
  #                               mc.cores=cores)
} else if (dataset == "greenleaf_pbmc_bm") {
  dir_path = "/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/bed_files/greenleaf_pbmc_bm_scATAC/migrated_to_hg19"
  dir.create(dir_path)
  files = list.files("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/bed_files/greenleaf_pbmc_bm_scATAC",
                     pattern="tsv.gz",
                     full.names=TRUE)
  for (l in split(files, ceiling(seq_along(files)/2))) {
    helper(l, dir_path, ch, cores)
  }
  # fragments = mclapply(files, import, format="bed", mc.cores=cores)
  # migrated_fragments = mclapply(fragments, migrate_bed_file_to_hg37, ch, 
  #                               mc.cores=cores)
}

