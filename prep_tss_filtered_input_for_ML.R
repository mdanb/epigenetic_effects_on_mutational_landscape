library(optparse)
library(dplyr)

option_list <- list( 
  make_option("--fragment_count_range", type="character"),
  make_option("--bing_ren", action="store_true", type="logical", default=F),
  make_option("--shendure", action="store_true", type="logical", default=F),                         
  make_option("--tsankov", action="store_true", type="logical", default=F)
)

args = parse_args(OptionParser(option_list=option_list))
# args = parse_args(OptionParser(option_list=option_list), args=
#                   c("--fragment_count_range=100000,200000"))

fragment_count_range = unlist(strsplit(args$fragment_count_range, ","))
bing_ren = args$bing_ren
shendure = args$shendure
tsankov = args$tsankov

root_tss_filtered = "/ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/processed_data/count_overlap_data/tsse_filtered"
combined_count_overlaps_path = "/ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/processed_data/count_overlap_data/combined_count_overlaps"
bingren_path = paste(root_tss_filtered, "bing_ren", sep="/")
shendure_path = paste(root_tss_filtered, "shendure", sep="/")
tsankov_path = paste(root_tss_filtered, "tsankov", sep="/")

combined_tss_filtered_cells <- function(fragment_count_range, tss_path, 
                                        combined_overlaps_filepath) {
  combined_overlaps_unfiltered = readRDS(combined_overlaps_filepath)
  combined_path = paste(tss_path, "combined", sep="/")
  dir.create(combined_path)
  for (fragment_count in fragment_count_range) {
    combined_df_list = list()
    for (dir in setdiff(list.dirs(tss_path, recursive = F, full.names = F),
                        "combined")) {
      df = readRDS(paste(paste(tss_path, dir, sep="/"), 
                         paste0("count_overlaps_frag_count_filter_",
                                fragment_count, 
                                ".rds"), 
                   sep="/"))
      if (nrow(df) == 0) {
        print(paste(dir, "has no count overlaps for", fragment_count, 
                    "fragments"), sep=" ")
        next
      }
      
      colnames(df) <- paste(dir, colnames(df))
      combined_df_list = append(combined_df_list, list(df))
    }
    combined_df = bind_cols(combined_df_list)
    for (cell_type in rownames(combined_overlaps_unfiltered)) {
      if (!(cell_type %in% colnames(combined_df))) {
        combined_df = cbind(unlist(combined_overlaps_unfiltered[cell_type, ]), 
                            combined_df)
        colnames(combined_df)[1] = cell_type
      }
    }
    filename = paste0("combined_", fragment_count, "_fragments",".rds")
    filepath = paste(combined_path, filename, sep="/")
    saveRDS(combined_df, filepath)
  }
}

if (bing_ren) {
  combined_overlaps_file = "count_filter_1_combined_count_overlaps.rds"
  combined_overlaps_filepath = paste(combined_count_overlaps_path, 
                                     combined_overlaps_file, sep="/")
	combined_tss_filtered_cells(fragment_count_range, bingren_path, 
	                            combined_overlaps_filepath)
}

if (shendure) {
  combined_overlaps_file = "shendure_count_filter_1_combined_count_overlaps.rds"
  combined_overlaps_filepath = paste(combined_count_overlaps_path, 
                                     combined_overlaps_file, sep="/")
	combined_tss_filtered_cells(fragment_count_range, shendure_path, 
	                            combined_overlaps_filepath)
}

if (tsankov) {
  combined_overlaps_file = "tsankov_count_filter_1_combined_count_overlaps.rds"
  combined_overlaps_filepath = paste(combined_count_overlaps_path, 
                                     combined_overlaps_file, sep="/")
	combined_tss_filtered_cells(fragment_count_range, tsankov_path, 
	                            combined_overlaps_filepath)
}
