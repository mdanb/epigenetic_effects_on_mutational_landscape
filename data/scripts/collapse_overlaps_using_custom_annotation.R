library(optparse)
####### CUSTOM ANNOTATION OPTIONS ####### 
### Greenleaf pbmc bm
# CD14-mono_CDlike-T_preB+B-B_late+early-no+distinction_Unk-rm

option_list <- list( 
  make_option("--dataset", type="character"),
  make_option("--annotation", type="character"),
  make_option("--cell_number_filter", type="integer")
)

args = parse_args(OptionParser(option_list=option_list))

cell_number_filter = args$cell_number_filter
dataset = args$dataset
annotation = args$annotation

collapse_using_mapping <- function(mapping, df) {
  for (pattern_replacement in mapping) {
    idxs = grep(pattern_replacement[[1]], rownames(df))
    collapsed_co = rep(0, length(colnames(df)))
    for (idx in idxs) {
      row_to_add = df[idx, ]
      collapsed_co = collapsed_co + row_to_add
    }
    df = df[-idxs, ]
    df[pattern_replacement[[2]], ] = collapsed_co
  }
  return(df)
}

root = "../processed_data/count_overlap_data/combined_count_overlaps"
if (dataset == "Greenleaf_pbmc_bm") {
  default_annotation_fn = paste("Greenleaf_pbmc_bm_count_filter", 
                                cell_number_filter, 
                                "combined_count_overlaps.rds", sep="_")
  default_annotation_fp = paste(root, "default_annotation", 
                                default_annotation_fn, sep="/")
  default_combined_combine_count_ovs = readRDS(default_annotation_fp)
  
  if (annotation == "CD14-mono_CDlike-T_preB+B-B_late+early-no+distinction_Unk-rm") {
    df = default_combined_combine_count_ovs
    # pattern --> replacement (collapsing)
    mapping = list(
                c("bonemarrow CD14.Mono", "bonemarrow Monocytes"),
                c("blood CD14\\.Mono", "blood Monocytes"), 
                c("bonemarrow CD", "bonemarrow T"),
                c("blood CD", "blood T"), 
                c("bonemarrow B|bonemarrow Pre\\.B", "bonemarrow B"),
                c("bonemarrow .*Eryth", "bonemarrow Eryth"))
    
    df = collapse_using_mapping(mapping, default_combined_combine_count_ovs)
    df = df[-grep("Unk", rownames(df)), ]
    save_dir = paste("../processed_data/count_overlap_data/combined_count_overlaps",
                annotation, sep="/")
    save_file = paste(dataset, "count_filter", cell_number_filter, 
                      "combined_count_overlaps.rds", sep = "_")
    dir.create(save_dir)
    saveRDS(df, paste(save_dir, save_file, sep="/"))
  } 
}