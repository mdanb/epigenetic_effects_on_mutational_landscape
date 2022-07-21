library(stringr)
source("utils.R")

load('raw_dir/mutation_data/hg19.1Mb.ranges.Polak.Nature2015.RData')
CELL_NUMBER_FILTER = 1

save_combined_overlaps <- function(filepaths,
                                   combined_filepath,
                                   unsquashed_overlaps_filepath=NULL,
                                   unsquashed_tissue_names_filepath=NULL) {
  combined_count_overlaps = data.frame()
  combined_count_overlaps_metadata = data.frame()
  unsquashed_count_overlaps = data.frame()
  unsquashed_tissue_names = c()
  for (f in filepaths) {
    print(paste("Processing count overlaps for file: ", f, sep=""))
    if (grepl("frontal_cortex", f)) {
      next
    }
    tissue_name <- get_tissue_name(f)
    count_overlaps = readRDS(f)
    cell_types = names(count_overlaps)
    count_overlaps = as.data.frame(do.call(rbind, count_overlaps),
                                   row.names = paste(tissue_name,
                                                     cell_types))
    # if (grepl("subdivided", f)) {
    #   colnames(count_overlaps) = rep(colnames(count_overlaps), each=5)
    # }
    
    if (!is.null(unsquashed_overlaps_filepath)) {
      unsquashed_tissue_names = append(unsquashed_tissue_names, rep(tissue_name,
                                                                    dim(count_overlaps)[1]))
      unsquashed_count_overlaps = rbind(unsquashed_count_overlaps, 
                                        count_overlaps)
    }
    
    for (i in 1:nrow(count_overlaps)) {
      count_overlaps_exists = rownames(count_overlaps)[i] %in%
        rownames(combined_count_overlaps)
      tryCatch(
        if (count_overlaps_exists) {
          idx = match(rownames(count_overlaps)[i],
                      rownames(combined_count_overlaps))
          combined_count_overlaps[idx, ] = combined_count_overlaps[idx, ] +
            count_overlaps[i, ]
        }
        else {
          combined_count_overlaps = rbind(combined_count_overlaps,
                                          count_overlaps[i, ])
          combined_count_overlaps_metadata[dim(combined_count_overlaps_metadata)[1] + 1, 
                                           "tissue_name"] = tissue_name[1]
          combined_count_overlaps_metadata[dim(combined_count_overlaps_metadata)[1], 
                                           "cell_type"] = cell_types[i]
        },
        error = function(err) {
          print(paste(f, "has no count overlaps"), sep=" ")
        })
    }
  }
  saveRDS(combined_count_overlaps, combined_filepath)
  saveRDS(combined_count_overlaps_metadata, 
          paste(unlist(strsplit(combined_filepath, "[.]"))[1], "metadata.rds", 
                sep="_"))
  if (!is.null(unsquashed_overlaps_filepath)) {
    saveRDS(unsquashed_count_overlaps, unsquashed_overlaps_filepath)
    saveRDS(unsquashed_tissue_names, unsquashed_tissue_names_filepath)
  }
}

save_combined_overlaps_shendure_tsankov <- function(filepaths,
                                                    combined_filepath) {
  combined_count_overlaps = data.frame()
  for (f in filepaths) {
    print(paste("Processing count overlaps for file: ", f, sep=""))
    if (grepl("Tsankov", f)) {
      tissue_name <- get_tissue_name_tsankov(f)
    }
    else {
      tissue_name <- get_tissue_name_shendure(f)
    }
    count_overlaps = readRDS(f)
    cell_types = names(count_overlaps)
    count_overlaps = as.data.frame(do.call(rbind, count_overlaps),
                                   row.names = paste(tissue_name,
                                                     cell_types))
    for (i in 1:nrow(count_overlaps)) {
      count_overlaps_exists = rownames(count_overlaps)[i] %in%
        rownames(combined_count_overlaps)
      tryCatch(
        if (count_overlaps_exists) {
          idx = match(rownames(count_overlaps)[i],
                      rownames(combined_count_overlaps))
          combined_count_overlaps[idx, ] = combined_count_overlaps[idx, ] +
            count_overlaps[i, ]
        }
        else {
          combined_count_overlaps = rbind(combined_count_overlaps,
                                          count_overlaps[i, ])
        },
        error = function(err) {
          print(paste(f, "has no count overlaps"), sep=" ")
        })
    }
  }
  saveRDS(combined_count_overlaps, combined_filepath)
}

combined_data_path = "processed_data/count_overlap_data/combined_count_overlaps/"
dir.create(combined_data_path)

combined_filepath = paste(combined_data_path, 
                          "count_filter_",
                          CELL_NUMBER_FILTER, 
                          "_combined_count_overlaps.rds", sep="")

combined_filepath_shendure = paste(combined_data_path, 
                                   "shendure_count_filter_",
                                   CELL_NUMBER_FILTER, 
                                   "_combined_count_overlaps.rds", sep="")

combined_filepath_tsankov = paste(combined_data_path, 
                                  "tsankov_count_filter_",
                                  CELL_NUMBER_FILTER, 
                                  "_combined_count_overlaps.rds", sep="")

unsquashed_overlaps_filepath = paste(combined_data_path, 
                                     "count_filter_",
                                     CELL_NUMBER_FILTER, 
                                     "_unsquashed_count_overlaps.rds", sep="")

unsquashed_tissue_names_filepath = paste(combined_data_path, 
                                         "count_filter_",
                                         CELL_NUMBER_FILTER, 
                                         "_unsquashed_tissue_names.rds", sep="")

combined_subdivided_filepath = paste(combined_data_path, 
                                     "count_filter_",
                                     CELL_NUMBER_FILTER, 
                                     "_subdivided_combined_count_overlaps.rds", sep="")

if (!file.exists(combined_filepath) || 
    !file.exists(unsquashed_overlaps_filepath)) {
  pattern = paste("^count_filter_", CELL_NUMBER_FILTER, "_count_overlaps", 
                  sep="")
  filepaths = list.files("processed_data/count_overlap_data", pattern = pattern, 
                         full.names = TRUE)
  save_combined_overlaps(filepaths, combined_filepath, 
                         unsquashed_overlaps_filepath,
                         unsquashed_tissue_names_filepath)
}

if (!file.exists(combined_filepath_shendure)) {
  pattern = paste("Shendure_count_filter_", CELL_NUMBER_FILTER, 
                  "_count_overlaps", sep="")
  filepaths = list.files("processed_data/count_overlap_data", pattern = pattern, 
                         full.names = TRUE)
  save_combined_overlaps_shendure_tsankov(filepaths, combined_filepath_shendure)
}

if (!file.exists(combined_filepath_tsankov)) {
  pattern = paste("Tsankov_count_filter_", CELL_NUMBER_FILTER, 
                  "_count_overlaps", sep="")
  filepaths = list.files("processed_data/count_overlap_data", pattern = pattern, 
                         full.names = TRUE)
  save_combined_overlaps_shendure_tsankov(filepaths, combined_filepath_tsankov)
}

# if (!file.exists(combined_subdivided_filepath)) {
#   pattern_subdivided = paste("count_filter_", CELL_NUMBER_FILTER, "_subdivided", 
#                              sep="")
#   filepaths_subdivided = list.files("processed_data/count_overlap_data", 
#                                     pattern = pattern_subdivided, 
#                                     full.names = TRUE)
#   save_combined_overlaps(filepaths_subdivided, combined_subdivided_filepath)
# }
