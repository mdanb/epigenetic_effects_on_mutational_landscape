library(tools)                                                                  
library(rtracklayer)                                                            
library(stringr)                                                                
library(data.table)                                                             
library(dplyr)                                                                  
library(tibble)                                                                 
library(exomeCopy)
library(parallel)
library(liftOver)

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 1) {
  stop("One argument must be supplied", call.=FALSE)
} else {
  CELL_NUMBER_FILTER = as.integer(args[1])
} 

get_sample_name <- function(file) {
  if (grepl("frontal_cortex", file)) {
    sample_name = "frontal_cortex"
  } 
  else {
    sample_name = str_extract(file, pattern="_.*_SM")              
    sample_name = substr(sample_name, 2, nchar(sample_name) - 3)    
  }
  return(sample_name)
}

get_sample_name_shendure <- function(file) {
  sample_name = unlist(strsplit(file, split="_"))[3]
  sample_name = unlist(strsplit(sample_name, split="[.]"))[1]
  return(sample_name)
}

filter_sample_by_cell_number <- function(sample, cell_number_filter) {
  counts_per_cell_type = sample %>%                               
    group_by(cell_type) %>%                 
    summarize(n_cells = n_distinct(name))    
  
  cell_type_keep = counts_per_cell_type[counts_per_cell_type$n_cells >= 
                                          CELL_NUMBER_FILTER, ]$cell_type
  sample = sample %>%                                             
    filter(cell_type %in% cell_type_keep)
  return(sample)
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

get_sample_cell_types <- function(sample, sample_name, metadata) {
  filtered_metadata = metadata[metadata$sample %like% sample_name, ]
  sample_cellIDs = filtered_metadata[["cellID"]]                  
  sample_barcodes_in_metadata = unlist(lapply(sample_cellIDs, 
                                              function(x) unlist(strsplit(x, 
                                                                    '\\+'))[2]))
  sample = sample[sample$name %in% sample_barcodes_in_metadata]   
  sample_idx_in_metadata = match(sample$name, sample_barcodes_in_metadata)
  sample$cell_type = filtered_metadata[unlist(sample_idx_in_metadata), 
                                       "cell.type"]
  return(as.tibble(sample))
}

get_sample_cell_types_shendure <- function(sample, sample_name, metadata) {
  filtered_metadata = metadata[metadata$sample %like% sample_name, ]
  sample_barcodes_in_metadata = filtered_metadata[["cell"]]                  
  sample = sample[sample$name %in% sample_barcodes_in_metadata]   
  sample_idx_in_metadata = match(sample$name, sample_barcodes_in_metadata)
  sample$cell_type = filtered_metadata[unlist(sample_idx_in_metadata), 
                                       "cell_type"]
  return(as.tibble(sample))
}

get_sample_cell_types_tsankov <- function(sample, metadata) {
  sample_barcodes_in_metadata = unlist(lapply(metadata[["X"]],
                                              function(x) unlist(strsplit(x,
                                                                    '#'))[2]))
  sample = sample[sample$name %in% sample_barcodes_in_metadata]
  sample_idx_in_metadata = match(sample$name, sample_barcodes_in_metadata)
  sample$cell_type = metadata[unlist(sample_idx_in_metadata),
                                       "celltypes"]
  return(as.tibble(sample))
}

migrate_bed_file_to_hg37 <- function(bed_sample, chain) {
  seqlevelsStyle(bed_sample) = "UCSC"  # necessary
  bed_sample = unlist (liftOver(bed_sample, ch))
  return(bed_sample)
}

create_count_overlaps_file <- function(file, cell_number_filter, metadata,
                                       interval_ranges, chain) {
  filename = paste("count_filter", CELL_NUMBER_FILTER,  
                   "count_overlaps", paste(file_path_sans_ext(file, TRUE),
                                           "rds", sep="."), sep="_")
  subdivided_filename = paste("count_filter", CELL_NUMBER_FILTER,  
                              "subdivided_count_overlaps", 
                              paste(file_path_sans_ext(file, TRUE),
                                    "rds", sep="."), sep="_")
  filepath = paste("processed_data/count_overlap_data", filename, sep="/")
  subdivided_filepath = paste("processed_data/count_overlap_data", 
                              subdivided_filename, 
                              sep="/")
  if (!file.exists(filepath) || !file.exists(subdivided_filepath)) {
    print(paste("Processing", file, sep= " "))
    sample = import(paste("raw_dir", "bed_files", file, sep="/"), format="bed")
    sample = migrate_bed_file_to_hg37(sample, chain)
    sample_name = get_sample_name(file)
    sample <- get_sample_cell_types(sample, sample_name, metadata)
    sample <- filter_sample_by_cell_number(sample, CELL_NUMBER_FILTER)
    if (!file.exists(filename)) {
      count_overlaps <- compute_count_overlaps(sample, interval_ranges)
      saveRDS(count_overlaps, filepath)
    }
    if (!file.exists(subdivided_filename)) {
      subdivided_interval_ranges = subdivideGRanges(interval.ranges, 
                                                    subsize = 200000)
        # names(subdivided_interval_ranges) = rep(colnames(count_overlaps), each=5)
      
      count_overlaps_subdivided <- compute_count_overlaps(sample, 
                                                     subdivided_interval_ranges)
      # colnames(count_overlaps) = rep(colnames(count_overlaps), each=5)
      
      saveRDS(count_overlaps_subdivided, subdivided_filepath)
    }
  }
}

create_count_overlaps_file_shendure <- function(file, cell_number_filter, 
                                                metadata, interval_ranges) {
  filename = paste("Shendure_count_filter", CELL_NUMBER_FILTER,  
                   "count_overlaps", paste(file_path_sans_ext(file, TRUE),
                                           "rds", sep="."), sep="_")
  filepath = paste("processed_data/count_overlap_data", filename, sep="/")
  if (!file.exists(filepath)) {
    print(paste("Processing", file, sep= " "))
    sample = import(paste("raw_dir", "bed_files", "JShendure_scATAC", file, 
                          sep="/"), format="bed")
    sample_name = get_sample_name_shendure(file)
    sample <- get_sample_cell_types_shendure(sample, sample_name, metadata)
    sample <- filter_sample_by_cell_number(sample, CELL_NUMBER_FILTER)
    if (!file.exists(filename)) {
      count_overlaps <- compute_count_overlaps(sample, interval_ranges)
      saveRDS(count_overlaps, filepath)
    }
  }
}

create_count_overlaps_file_tsankov <- function(file, cell_number_filter, metadata,
                                       interval_ranges, chain) {
  filename = paste(file_path_sans_ext(file, TRUE), "rds", sep=".")
  filename = unlist(strsplit(filename, split = "_"))
  filename = paste(filename[2:length(filename)], collapse="_")
  filename = paste("Tsankov_count_filter", CELL_NUMBER_FILTER,  
                   "count_overlaps", filename, sep="_")
  filepath = paste("processed_data/count_overlap_data", filename, sep="/")
  if (!file.exists(filepath)) {
    print(paste("Processing", file, sep= " "))
    sample = import(paste0(paste("raw_dir", "bed_files", "Tsankov_scATAC", 
                          substr(file, 1, nchar(file) - 3), sep="/"), "tsv"),
                          format="bed")
    sample = migrate_bed_file_to_hg37(sample, chain)
    sample <- get_sample_cell_types_tsankov(sample, metadata)
    sample = filter_sample_by_cell_number(sample, CELL_NUMBER_FILTER)
    if (!file.exists(filename)) {
      count_overlaps = compute_count_overlaps(sample, interval_ranges)
      saveRDS(count_overlaps, filepath)
    }
  }
}

metadata = read.table("raw_dir/metadata/GSE184462_metadata.tsv", sep="\t", 
                      header=TRUE)
metadata_Shendure = read.table("raw_dir/metadata/GSE149683_File_S2.Metadata_of_high_quality_cells.txt", 
                               sep="\t", 
                               header=TRUE)
metadata_tsankov_proximal = read.table("raw_dir/metadata/tsankov_lung_proximal_barcode_annotation.csv", 
                                       sep=",", 
                                       header=TRUE)
metadata_tsankov_distal = read.table("raw_dir/metadata/tsankov_lung_distal_barcode_annotation.csv", 
                                       sep=",", 
                                       header=TRUE)
colnames(metadata_tsankov_proximal)[2] <- "celltypes"

load('raw_dir/mutation_data/hg19.1Mb.ranges.Polak.Nature2015.RData')
interval.ranges

dir.create("processed_data/count_overlap_data", recursive=TRUE)                                       

files = setdiff(list.files("raw_dir/bed_files/"), 
                list.dirs("raw_dir/bed_files", recursive = FALSE, 
                          full.names = FALSE))
files_Shendure = setdiff(list.files("raw_dir/bed_files/JShendure_scATAC/"), 
                         list.dirs("raw_dir/bed_files/JShendure_scATAC/", 
                                   recursive = FALSE, 
                                   full.names = FALSE))

files_Tsankov_distal = list.files("raw_dir/bed_files/Tsankov_scATAC/", pattern=".*distal.*")
files_Tsankov_proximal = list.files("raw_dir/bed_files/Tsankov_scATAC/", pattern=".*proximal*")

hg38_path = system.file (package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch = import.chain(hg38_path)

#lapply(files, create_count_overlaps_file,
#                cell_number_filter=CELL_NUMBER_FILTER,
#                metadata=metadata,
#                interval_ranges=interval.ranges,
#                chain=ch)
#                mc.cores = 1)

lapply(files_Shendure, create_count_overlaps_file_shendure,
         cell_number_filter=CELL_NUMBER_FILTER,
         metadata=metadata_Shendure,
         interval_ranges=interval.ranges)
#         mc.cores = 1)

#lapply(files_Tsankov_proximal, create_count_overlaps_file_tsankov,
#         cell_number_filter=CELL_NUMBER_FILTER,
#         metadata=metadata_tsankov_proximal,
#         interval_ranges=interval.ranges)
#         mc.cores = 1)

#lapply(files_Tsankov_distal, create_count_overlaps_file_tsankov, 
#         cell_number_filter=CELL_NUMBER_FILTER,
#         metadata=metadata_tsankov_distal,
#         interval_ranges=interval.ranges)
#         #mc.cores = 1)
