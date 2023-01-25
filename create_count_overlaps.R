library(tools)                                                                  
library(stringr)                                                                
library(data.table)                                                             
library(dplyr)                                                                  
library(tibble)                                                                 
library(parallel)
library(optparse)
library(rtracklayer)
library(readxl)
source("create_count_overlaps_utils.R")

option_list <- list( 
  make_option("--dataset", type="character"),
  make_option("--cell_number_filter", type="integer"),
  make_option("--cores", type="integer")
)

args = parse_args(OptionParser(option_list=option_list))
#args = parse_args(OptionParser(option_list=option_list), args =
#                  c("--dataset=shendure",
#                    "--cell_number_filter=1",
#                    "--cores=1"))

cell_number_filter = args$cell_number_filter
cores = args$cores
dataset = args$dataset

get_and_save_num_cells_per_sample <- function(sample, sample_file_name) {
  sample_file_name = paste0("cell_counts_", file_path_sans_ext(sample_file_name, 
                                                               TRUE), ".rds")
  counts_per_cell_type = get_num_cells_per_sample(sample)
  path = "processed_data/cell_counts_per_sample"
  file_path = paste(path, sample_file_name, sep="/")
  saveRDS(counts_per_cell_type, file_path)
  return(counts_per_cell_type)
}

compute_count_overlaps <- function(sample, interval_ranges) {
  grl_in = sample %>%                                             
    group_split(cell_type) %>%                             
    lapply(makeGRangesFromDataFrame, keep.extra.columns=TRUE)
  grl = GRangesList(grl_in)
  count_overlaps = lapply(grl, function(x) countOverlaps (interval_ranges, x))
  cell_types = unlist(lapply(grl, function(x) unique(x$cell_type)))
  names(count_overlaps) = cell_types
  return(count_overlaps)
}

create_count_overlaps_files <- function(file, cell_number_filter, metadata,
                                  interval_ranges, chain, dataset) {
  filename = get_sample_filename(file, dataset)
  filepath = paste("processed_data/count_overlap_data", filename, sep="/")
    if (!file.exists(filepath)) {
      print(paste("Processing", file, sep= " "))
      sample = import_sample(file, dataset)
      sample_name = get_sample_name(file, dataset)
      filtered_metadata = filter_metadata_by_sample_name(sample_name, metadata)
      if (nrow(filtered_metadata) == 0) {
        print(paste(file, "has no high quality fragments", sep=" "))
        return()
      }
      
      sample_barcodes_in_metadata = get_sample_barcodes_in_metadata(filtered_metadata,
                                                                    dataset)
      
      if (dataset == "tsankov" || dataset == "greenleaf_brain") {
        sample$name = substr(sample$name, 1, 16)
        sample = migrate_bed_file_to_hg37(sample, chain)
      }
      else if (dataset == "bingren") {
        sample = migrate_bed_file_to_hg37(sample, chain)
      }
      
      # Because the function applied is used in another place over multiple 
      # samples, not just one, as is done here.
      sample <- unlist(lapply(c(1), 
                              filter_samples_to_contain_only_cells_in_metadata,
                              list(sample),
                              list(sample_barcodes_in_metadata)))
      
      sample <- get_sample_cell_types(sample[[1]], sample_barcodes_in_metadata,
                                      filtered_metadata)
      counts_per_cell_type <- get_and_save_num_cells_per_sample(sample, file)
      sample <- filter_sample_by_cell_number(sample,
                                             counts_per_cell_type, 
                                             cell_number_filter)
      
      count_overlaps <- compute_count_overlaps(sample, interval_ranges)
      saveRDS(count_overlaps, filepath)
    }
}

filter_sample_by_cell_number <- function(sample, 
                                         counts_per_cell_type, 
                                         cell_number_filter) {
  cell_type_keep = counts_per_cell_type[counts_per_cell_type$n_cells >= 
                                          cell_number_filter, ]$cell_type
  sample = sample %>%                                             
           filter(cell_type %in% cell_type_keep)
  return(sample)
}

get_num_cells_per_sample <- function(sample) {
  counts_per_cell_type = sample %>%                               
                         group_by(cell_type) %>%                 
                         summarize(n_cells = n_distinct(name))
  return(counts_per_cell_type)
}

get_sample_cell_types <- function(fragments, 
                                  sample_barcodes_in_metadata,
                                  filtered_metadata) {
  sample_idx_in_metadata = match(fragments$name, 
                                 sample_barcodes_in_metadata)
  fragments$cell_type = filtered_metadata[unlist(sample_idx_in_metadata), 
                                          "cell_type"]
  return(as.tibble(fragments))
}


get_sample_cell_types_bingren <- function(fragments, 
                                  sample_barcodes_in_metadata,
                                  filtered_metadata) {
  sample_idx_in_metadata = match(fragments$name, 
                                 sample_barcodes_in_metadata)
  fragments$cell_type = filtered_metadata[unlist(sample_idx_in_metadata), 
                                           "cell.type"]
  return(as.tibble(fragments))
}

get_sample_cell_types_shendure <- function(sample, 
                                           sample_barcodes_in_metadata,
                                           filtered_metadata) {
  sample_idx_in_metadata = match(sample$name, sample_barcodes_in_metadata)
  sample$cell_type = filtered_metadata[unlist(sample_idx_in_metadata), 
                                       "cell_type"]
  return(as.tibble(sample))
}

get_sample_cell_types_tsankov <- function(sample,
                                          sample_barcodes_in_metadata,
                                          metadata) {
  sample_idx_in_metadata = match(sample$name, sample_barcodes_in_metadata)
  sample$cell_type = metadata[unlist(sample_idx_in_metadata),
                              "celltypes"]
  return(as.tibble(sample))
}


load('/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/mutation_data/hg19.1Mb.ranges.Polak.Nature2015.RData')

dir.create("processed_data/count_overlap_data", recursive=TRUE)                                       
dir.create("processed_data/cell_counts_per_sample")                                       
ch = import.chain("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/hg38ToHg19.over.chain")

if (dataset == "bingren") {
  metadata_bingren = read.table("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/metadata/GSE184462_metadata.tsv", 
                        sep="\t",
                        header=T)
  colnames(metadata_bingren)[grepl("cell.type", colnames(metadata_bingren))] = "cell_type"
  files_bingren = list.files("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/bed_files/",
                              pattern=".*fragments\\.bed\\.gz")
  
  mclapply(files_bingren, create_count_overlaps_files,
           cell_number_filter=cell_number_filter,
           metadata=metadata_bingren,
           interval_ranges=interval.ranges,
           chain=ch,
           dataset=dataset,
           mc.cores=cores)
} else if (dataset == "shendure") {
  metadata_Shendure = read.table("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/metadata/GSE149683_File_S2.Metadata_of_high_quality_cells.txt",
                                  sep="\t",
                                  header=TRUE)
  files_Shendure = list.files("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/bed_files/JShendure_scATAC/",
                              pattern = ".*fragments\\.txt\\.gz")
  
  mclapply(files_Shendure, create_count_overlaps_files,
         cell_number_filter=cell_number_filter,
         metadata=metadata_Shendure,
         interval_ranges=interval.ranges,
         chain=ch,
         dataset=dataset,
         mc.cores=cores)
} else if (dataset == "tsankov") {
  metadata_tsankov_proximal = 
    read.csv("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/metadata/tsankov_lung_proximal_barcode_annotation.csv")
  metadata_tsankov_distal = 
    read.csv("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/metadata/tsankov_lung_distal_barcode_annotation.csv")
  colnames(metadata_tsankov_proximal)[grepl("celltypes",
                                            colnames(metadata_tsankov_proximal))] = "cell_type"
  colnames(metadata_tsankov_proximal)[grepl("Sample",
                                            colnames(metadata_tsankov_proximal))] = "sample"
  colnames(metadata_tsankov_distal)[grepl("celltypes",
                                            colnames(metadata_tsankov_distal))] = "cell_type"
  colnames(metadata_tsankov_distal)[grepl("Sample",
                                            colnames(metadata_tsankov_distal))] = "sample"
  

  files_Tsankov_distal = list.files("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/bed_files/Tsankov_scATAC/", 
                                    pattern="RPL")
  files_Tsankov_proximal = list.files("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/bed_files/Tsankov_scATAC/", 
                                      pattern="IC")
  mclapply(files_Tsankov_proximal, 
           create_count_overlaps_files,
           cell_number_filter=cell_number_filter,
           metadata=metadata_tsankov_proximal,
           interval_ranges=interval.ranges,
           chain=ch,
           dataset=dataset,
           mc.cores=cores)
  
  mclapply(files_Tsankov_distal,
           create_count_overlaps_files,
            cell_number_filter=cell_number_filter,
            metadata=metadata_tsankov_distal,
            interval_ranges=interval.ranges,
            chain=ch,
            dataset=dataset,
            mc.cores=cores)
} else if (dataset == "greenleaf_brain") {
    if (!(file.exists("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/metadata/GSE162170_atac_cell_metadata_with_cell_names.txt"))) {
      metadata_greenleaf_brain =
        read.csv("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/metadata/GSE162170_atac_cell_metadata.txt.gz",
                 sep="\t")
      cluster_to_cell_names = read.csv("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/scATAC_Colors_greenleaf_brain.txt",
                                       sep = "\t") 
      cell_name_idx = match(metadata_greenleaf_brain[["Iterative.LSI.Clusters"]], 
                         cluster_to_cell_names[["cluster"]])
      cell_types = cluster_to_cell_names[cell_name_idx, "class"]
      metadata_greenleaf_brain["cell_type"] = cell_types
      write.csv(metadata_greenleaf_brain,
                "/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/metadata/GSE162170_atac_cell_metadata_with_cell_names.txt")
    }
    else {
      metadata_greenleaf_brain =
        read.csv("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/metadata/GSE162170_atac_cell_metadata_with_cell_names.txt")
    }

  colnames(metadata_greenleaf_brain)[grepl("Sample.ID", colnames(metadata_greenleaf_brain))] = "sample"
  files_greenleaf_brain = list.files("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/bed_files/greenleaf_brain_scATAC/", 
                                     pattern=".*fragments\\.tsv\\.gz")
  # mclapply(files_greenleaf_brain, create_count_overlaps_files,
  #          cell_number_filter=cell_number_filter,
  #          metadata=metadata_greenleaf_brain,
  #          interval_ranges=interval.ranges,
  #          chain=ch,
  #          dataset=dataset,
  #          mc.cores=cores)
  lapply(files_greenleaf_brain, create_count_overlaps_files,
           cell_number_filter=cell_number_filter,
           metadata=metadata_greenleaf_brain,
           interval_ranges=interval.ranges,
           chain=ch,
           dataset=dataset)
} else if (dataset == "greenleaf_pbmc_bm") {
    if (!(file.exists("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/metadata/greenleaf_pbmc_bm.txt"))) {
      metadata_greenleaf_pbmc_bm = readRDS("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/metadata/scATAC-Healthy-Hematopoiesis-191120.rds")
      metadata_greenleaf_pbmc_bm = data.frame(barcode = metadata_greenleaf_pbmc_bm$Barcode, 
                                              cell_type = unlist(lapply(
                                              strsplit(metadata_greenleaf_pbmc_bm$BioClassification, 
                                                       split="_"), "[", 2)))
      
      write.csv(metadata_greenleaf_pbmc_bm, 
                "/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/metadata/greenleaf_pbmc_bm.txt")
    }
    else {
      metadata_greenleaf_pbmc_bm = 
        read.csv("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/metadata/greenleaf_pbmc_bm.txt")
    }
    files_greenleaf_pbmc_bm = list.files("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/bed_files/greenleaf_pbmc_bm_scATAC", 
                                         pattern=".*fragments\\.tsv\\.gz")
    mclapply(files_greenleaf_pbmc_bm, 
             create_count_overlaps_files,
             cell_number_filter=cell_number_filter,
             metadata=metadata_greenleaf_pbmc_bm,
             interval_ranges=interval.ranges,
             chain=ch,
             dataset=dataset,
             mc.cores=cores)
}
