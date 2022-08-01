library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(tibble)
library(tidyr)
library(stringr)
library(exomeCopy)
library(parallel)
source("utils.R")

load('raw_dir/mutation_data/hg19.1Mb.ranges.Polak.Nature2015.RData')
dir.create("figures") 
CELL_NUMBER_FILTER = 1

load_mutation_data <- function() {
  mut = readRDS('raw_dir/mutation_data/mutations.aggregated.PCAWG.RK.20181215.Rds')
  # get all mutations, without distinctions e.g clonal vs subclonal 
  mut = mut[, 1:37]
  return(mut)
}

get_per_cancer_mut_data <- function(all_mutations, interval_ranges) {
  cancer_types_mut_data <- list()
  for (cancer_type in colnames(all_mutations)) {
    idx_cancer_type = match(rownames(as.data.frame(interval_ranges)), 
                            rownames(all_mutations[cancer_type]))
    cancer_types_mut_data = append(cancer_types_mut_data,
                                   list(mut[cancer_type][idx_cancer_type, ]))
  }
  return(cancer_types_mut_data)
}

save_heatmap_file <- function(path, width, height, df, ...) {
  pdf(path, width=width, height=height)
  h = Heatmap(df, show_column_dend = FALSE, show_row_dend = FALSE, ...)
  print(h)
  dev.off()
}

create_and_save_mutation_df_all_cancers <- function() {
  if (!file.exists("processed_data/mut_count_data.csv")) {
    cancer_types_mut_data = get_per_cancer_mut_data(mut, interval.ranges)
    mut_count_data = as.data.frame(do.call(cbind, cancer_types_mut_data))
    colnames(mut_count_data) = colnames(mut)
    rownames(mut_count_data) = names(interval.ranges)
    write.csv(mut_count_data, "processed_data/mut_count_data.csv")
  }
  else {
    mut_count_data = read.csv("processed_data/mut_count_data.csv")
  }
  return(mut_count_data)
}

mut = load_mutation_data()
mut_count_data = create_and_save_mutation_df_all_cancers()
rownames(mut_count_data) = mut_count_data[, 1]
mut_count_data = mut_count_data[, 2:length(colnames(mut_count_data))]
  
combined_data_path = "processed_data/count_overlap_data/combined_count_overlaps/"

combined_filepath = paste(combined_data_path,
                          "count_filter_",
                          CELL_NUMBER_FILTER,
                          "_combined_count_overlaps.rds", sep="")
combined_count_overlaps = t(readRDS(combined_filepath))

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

##### Cell type Cell type correlation ####
