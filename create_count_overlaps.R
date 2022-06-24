library(tools)                                                                  
library(rtracklayer)                                                            
library(stringr)                                                                
library(data.table)                                                             
library(dplyr)                                                                  
library(tibble)                                                                 
library(exomeCopy)
library(parallel)

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

create_count_overlaps_file <- function(file, cell_number_filter, metadata,
                                       interval_ranges) {
  filename = paste("count_filter", CELL_NUMBER_FILTER,  
                   "count_overlaps", paste(file_path_sans_ext(file, TRUE),
                                           "rds", sep="."), sep="_")
  subdivided_filename = paste("count_filter", CELL_NUMBER_FILTER,  
                              "subdivided_count_overlaps", 
                              paste(file_path_sans_ext(file, TRUE),
                                    "rds", sep="."), sep="_")
  filepath = paste("count_overlap_data", filename, sep="/")
  subdivided_filepath = paste("count_overlap_data", 
                              subdivided_filename, 
                              sep="/")
  if (!file.exists(filepath) || !file.exists(subdivided_filepath)) {
    print(paste("Processing", file, sep= " "))
    sample = import(paste("raw_dir", "bed_files", file, sep="/"), format="bed")
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
      count_overlaps_subdivided <- compute_count_overlaps(sample, 
                                                     subdivided_interval_ranges)
      saveRDS(count_overlaps_subdivided, subdivided_filepath)
    }
  }
}

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 1) {
  stop("One argument must be supplied", call.=FALSE)
} else {
  CELL_NUMBER_FILTER = as.integer(args[1])
} 

metadata = read.table("raw_dir/metadata/GSE184462_metadata.tsv", sep="\t", 
                      header=TRUE)  
load('raw_dir/mutation_data/hg19.1Mb.ranges.Polak.Nature2015.RData')
interval.ranges

dir.create("count_overlap_data")                                                

files = setdiff(list.files("raw_dir/bed_files/"), 
                list.dirs("raw_dir/bed_files", recursive = FALSE, 
                          full.names = FALSE))

mclapply(files, create_count_overlaps_file, 
                cell_number_filter=CELL_NUMBER_FILTER,
                metadata=metadata,
                interval_ranges=interval.ranges,
                mc.cores = 8)
