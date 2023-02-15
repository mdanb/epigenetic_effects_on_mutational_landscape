library(tidyverse)
library(optparse)
source("count_overlaps_utils.R")

load('/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/mutation_data/hg19.1Mb.ranges.Polak.Nature2015.RData')

option_list <- list( 
  make_option("--datasets", type = "character"),
  make_option("--cell_number_filter", type="integer"),
  make_option("--annotation", type="character")
)

args = parse_args(OptionParser(option_list=option_list))
# args = parse_args(OptionParser(option_list=option_list), args =
#                     c("--datasets=Greenleaf_pbmc_bm",
#                       "--cell_number_filter=1"))

annotation = args$annotation
cell_number_filter = args$cell_number_filter
datasets = unlist(strsplit(args$datasets, split = ","))

get_cell_counts_df <- function(count_overlaps_filename, annotation) {
  cell_counts_filename = unlist(strsplit(unlist(
                                         strsplit(count_overlaps_filename,
                                         split="/"))[6], 
                                         "_"))
  start_index = grep("overlaps", cell_counts_filename) + 1
  cell_counts_filename = 
    cell_counts_filename[start_index:length(cell_counts_filename)]
  cell_counts_filename = paste(cell_counts_filename, collapse="_")
  if (grepl("IC", cell_counts_filename) || grepl("RPL", 
                                                       cell_counts_filename)) {
    cell_counts_filename = paste(unlist(strsplit(cell_counts_filename, "[.]"))[1], 
                                 "fragments.rds", sep = "_")
    # cell_counts_filename = paste("fragments", cell_counts_filename, sep="_")
  }
  cell_counts_filename = paste("cell_counts", cell_counts_filename, sep="_")
  cell_counts_path = paste("..", "processed_data", 
                           "cell_counts_per_sample", annotation,
                           cell_counts_filename, sep = "/")
  df = readRDS(cell_counts_path)
  return(df)  
}

add_to_combined_dataframes <- function(count_overlaps, combined_count_overlaps,
                                       tissue_name, cell_types, cell_counts,
                                       combined_count_overlaps_metadata, 
                                       f) {
  if (nrow(count_overlaps) == 0) {
    print(paste(f, "has no count overlaps"), sep=" ")
    return(list(combined_count_overlaps, combined_count_overlaps_metadata))
  }
  for (i in 1:nrow(count_overlaps)) {
    count_overlaps_exists = rownames(count_overlaps)[i] %in%
                            rownames(combined_count_overlaps)
    # tryCatch(
    if (count_overlaps_exists) {
      idx = match(rownames(count_overlaps)[i],
                  rownames(combined_count_overlaps))
      combined_count_overlaps[idx, ] = combined_count_overlaps[idx, ] +
                                       count_overlaps[i, ]
      idx_tissue = match(tissue_name[1],
                         combined_count_overlaps_metadata[, "tissue_name"])
      idx_cell_type = match(cell_types[i],
                            combined_count_overlaps_metadata[, "cell_type"])
      num_cells = cell_counts[cell_counts["cell_type"] == cell_types[i], 
                              "n_cells"] %>% pull(n_cells)
      n_cells_idx = intersect(idx_tissue, idx_cell_type)
      combined_count_overlaps_metadata[n_cells_idx, "num_cells"] = 
        combined_count_overlaps_metadata[n_cells_idx, "num_cells"] + num_cells
    }
    else {
      combined_count_overlaps = rbind(combined_count_overlaps,
                                      count_overlaps[i, ])
      combined_count_overlaps_metadata[dim(combined_count_overlaps_metadata)[1] + 1, 
                                       "tissue_name"] = tissue_name[1]
      combined_count_overlaps_metadata[dim(combined_count_overlaps_metadata)[1], 
                                       "cell_type"] = cell_types[i]
      n_cells = cell_counts[cell_counts["cell_type"] == cell_types[i], 
                            "n_cells"]
      combined_count_overlaps_metadata[dim(combined_count_overlaps_metadata)[1], 
                                       "num_cells"] = n_cells
    }
      # },
      # error = function(err) {
      #   print(paste(f, "has no count overlaps"), sep=" ")
      # })
  }
  return(list(combined_count_overlaps, combined_count_overlaps_metadata))
}

save_combined_overlaps <- function(filepaths,
                                   combined_filepath, 
                                   dataset, annotation) {
  # unsquashed_overlaps_filepath=NULL,
  # unsquashed_tissue_names_filepath=NULL) {
  combined_count_overlaps = data.frame()
  combined_count_overlaps_metadata = data.frame()
  for (f in filepaths) {
    print(paste("Processing count overlaps for file: ", f, sep=""))
    count_overlaps = readRDS(f)
    tissue_name <- get_tissue_name(f, dataset)
    cell_types = names(count_overlaps)
    count_overlaps = as.data.frame(do.call(rbind, count_overlaps),
                                   row.names = paste(tissue_name,
                                                     cell_types))
    cell_counts = get_cell_counts_df(f, annotation)  
    dfs = add_to_combined_dataframes(count_overlaps, combined_count_overlaps,
                                     tissue_name, cell_types, cell_counts,
                                     combined_count_overlaps_metadata,
                                     f)    
    combined_count_overlaps = dfs[[1]]
    combined_count_overlaps_metadata = dfs[[2]]
    
    # if (grepl("subdivided", f)) {
    #   colnames(count_overlaps) = rep(colnames(count_overlaps), each=5)
    # }
    
    # if (!is.null(unsquashed_overlaps_filepath)) {
    #   unsquashed_tissue_names = append(unsquashed_tissue_names, rep(tissue_name,
    #                                                     dim(count_overlaps)[1]))
    #   unsquashed_count_overlaps = rbind(unsquashed_count_overlaps, 
    #                                     count_overlaps)
    # }
  }
  saveRDS(combined_count_overlaps, combined_filepath)
  temp = unlist(strsplit(combined_filepath, "/"))[7]
  metadata_filename = paste(unlist(strsplit(temp, "[.]"))[1], "metadata.rds", 
                            sep="_")
  metadata_filepath = paste("..", "processed_data", "count_overlap_data", 
                            "combined_count_overlaps", annotation, 
                            metadata_filename, sep="/")
  saveRDS(combined_count_overlaps_metadata, metadata_filepath)
}

combined_data_path = 
  "../processed_data/count_overlap_data/combined_count_overlaps"

# unsquashed_overlaps_filepath = paste(combined_data_path,
#                                      "count_filter_",
#                                      cell_number_filter,
#                                      "_unsquashed_count_overlaps.rds", sep="")
# 
# unsquashed_tissue_names_filepath = paste(combined_data_path,
#                                          "count_filter_",
#                                          cell_number_filter,
#                                          "_unsquashed_tissue_names.rds", sep="")
# 
# combined_subdivided_filepath = paste(combined_data_path,
#                                      "count_filter_",
#                                      cell_number_filter,
#                                      "_subdivided_combined_count_overlaps.rds",
#                                      sep="")

# dir.create("../../data/processed_data/cell_counts_per_sample/combined_cell_counts")
for (dataset in datasets) {
  combined_filepath = paste(combined_data_path, annotation, sep="/")
  dir.create(combined_filepath)
  combined_filename = paste(dataset,
                            "count_filter",
                            cell_number_filter, 
                            "combined_count_overlaps.rds", sep="_")
  combined_filepath = paste(combined_filepath, combined_filename, sep="/")
  if (!file.exists(combined_filepath)) { # || 
      #!file.exists(unsquashed_overlaps_filepath)) {
    pattern = paste(dataset, "count_filter", cell_number_filter, "count_overlaps", 
                    sep="_")
    co_fp = paste("../processed_data/count_overlap_data", annotation,
                  sep = "/")
    filepaths = list.files(co_fp, 
                           pattern = pattern, 
                           full.names = TRUE)
    save_combined_overlaps(filepaths, combined_filepath, dataset)
  }
}
