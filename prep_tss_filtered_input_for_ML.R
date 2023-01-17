library(optparse)
library(dplyr)

option_list <- list( 
  make_option("--fragment_count_range", type="character"),
  make_option("--bing_ren", action="store_true", type="logical", default=F),
  make_option("--shendure", action="store_true", type="logical", default=F),                         
  make_option("--tsankov", action="store_true", type="logical", default=F)
)

args = parse_args(OptionParser(option_list=option_list))
#args = parse_args(OptionParser(option_list=option_list), args=
#                   c("--fragment_count_range=50000,100000,300000"))

fragment_count_range = unlist(strsplit(args$fragment_count_range, ","))
bing_ren = args$bing_ren
shendure = args$shendure
tsankov = args$tsankov

root = "/ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/processed_data/count_overlap_data/tsse_filtered"
bingren_path = paste(root, "bing_ren", sep="/")
shendure_path = paste(root, "shendure", sep="/")
tsankov_path = paste(root, "tsankov", sep="/")

combined_tss_filtered_cells <- function(fragment_count_range, tss_path) {
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
    filename = paste0("combined_", fragment_count, "_fragments",".rds")
    filepath = paste(combined_path, filename, sep="/")
    saveRDS(combined_df, filepath)
  }
}

if (bing_ren) {
	combined_tss_filtered_cells(fragment_count_range, bingren_path)
}

if (shendure) {
	combined_tss_filtered_cells(fragment_count_range, shendure_path)
}

if (tsankov) {
	combined_tss_filtered_cells(fragment_count_range, tsankov_path)
}
