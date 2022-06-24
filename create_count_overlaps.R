library(tools)                                                                  
library(rtracklayer)                                                            
library(stringr)                                                                
library(data.table)                                                             
library(dplyr)                                                                  
library(tibble)                                                                 
library(exomeCopy)

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

get_count_overlaps <- function(sample, interval_ranges) {
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

#load('/ahg/regevdata/projects/ICA_Lung/Bruno/celloforigin_prj/hg19.1Mb.ranges.Polak.Nature2015.RData')
load('raw_dir/mutation_data/hg19.1Mb.ranges.Polak.Nature2015.RData')

interval.ranges
subdivided_interval_ranges = subdivideGRanges(interval.ranges, subsize = 200000)
dir.create("count_overlap_data")                                                

metadata = read.table("raw_dir/metadata/GSE184462_metadata.tsv", sep="\t", 
                      header=TRUE)  

for (file in setdiff(list.files("raw_dir/bed_files/"), 
                     list.dirs("raw_dir/bed_files", recursive = FALSE, 
                               full.names = FALSE))) {                                    
  filename = paste("count_filter", CELL_NUMBER_FILTER,  
                   "count_overlaps", paste(file_path_sans_ext(file, TRUE),
                                           "rds", sep="."), sep="_")
  subdivided_filename = paste("count_filter", CELL_NUMBER_FILTER,  
                              "subdivided_count_overlaps", paste(file_path_sans_ext(file, TRUE),
                                                                 "rds", sep="."), sep="_")
  if (!file.exists(filename) || !file.exists(subdivided_filename)) {
    print(paste("Processing", file, sep= " "))
    sample = import(paste("raw_dir", "bed_files", file, sep="/"), format="bed")
    sample_name = get_sample_name(file)
    sample <- get_sample_cell_types(sample, sample_name, metadata)
    sample <- filter_sample_by_cell_number(sample, CELL_NUMBER_FILTER)
    count_overlaps <- get_count_overlaps(sample, interval.ranges)
    count_overlaps_subdivided <- get_count_overlaps(sample, 
                                                    subdivided_interval_ranges)
    saveRDS(count_overlaps, paste("count_overlap_data", filename, sep="/"))
    saveRDS(count_overlaps_subdivided, paste("count_overlap_data", 
                                             subdivided_filename, 
                                             sep="/"))
  }
}
