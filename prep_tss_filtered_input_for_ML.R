library(optparse)
library(dplyr)

option_list <- list( 
  make_option("--fragment_count_range", type="character")
)

args = parse_args(OptionParser(option_list=option_list))
# args = parse_args(OptionParser(option_list=option_list), args=
#                    c("--fragment_count_range=100000,150000,200000"))

fragment_count_range = unlist(strsplit(args$fragment_count_range, ","))

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
      df = readRDS(list.files(paste(tss_path, dir, sep="/"), 
                              pattern=paste0(fragment_count, ".rds"), 
                              full.names = T))
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

combined_tss_filtered_cells(fragment_count_range, bingren_path)
combined_tss_filtered_cells(fragment_count_range, shendure_path)
combined_tss_filtered_cells(fragment_count_range, tsankov_path)

