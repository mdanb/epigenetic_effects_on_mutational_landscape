library(optparse)
####### CUSTOM ANNOTATION OPTIONS ####### 
### Greenleaf pbmc bm
# Greenleaf_pbmc_bm_CD14-mono_CDlike-T_preB+B-B_late+early-no+distinction_Unk-rm

### Greenleaf brain
# Greenleaf_brain_same+as+paper+but+Early+RG+Late+RG-RG_Unk-rm

option_list <- list( 
  make_option("--dataset", type="character"),
  make_option("--annotation", type="character"),
  make_option("--cell_number_filter", type="integer")
)

args = parse_args(OptionParser(option_list=option_list))

cell_number_filter = args$cell_number_filter
dataset = args$dataset
annotation = args$annotation

save_collapsed_df <- function(df, dataset, annotation, cell_number_filter) {
  save_dir = paste("../processed_data/count_overlap_data/combined_count_overlaps",
                   annotation, sep="/")
  save_file = paste(dataset, "count_filter", cell_number_filter, 
                    "combined_count_overlaps.rds", sep = "_")
  dir.create(save_dir)
  saveRDS(df, paste(save_dir, save_file, sep="/"))
}

collapse_using_mapping <- function(mapping, df) {
  for (pattern_replacement in mapping) {
    idxs = matching_func(pattern_replacement[1], rownames(df))
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
  default_combined_count_ovs = readRDS(default_annotation_fp)
  
  if (annotation == "Greenleaf_pbmc_bm_CD14-mono_CDlike-T_preB+B-B_late+early-no+distinction_Unk-rm") {
    df = default_combined_count_ovs
    # pattern --> replacement (collapsing)
    mapping = list(
                c("bonemarrow CD14.Mono", "bonemarrow Monocytes"),
                c("blood CD14\\.Mono", "blood Monocytes"), 
                c("bonemarrow CD", "bonemarrow T"),
                c("blood CD", "blood T"), 
                c("bonemarrow B|bonemarrow Pre\\.B", "bonemarrow B"),
                c("bonemarrow .*Eryth", "bonemarrow Eryth"))
    
    df = collapse_using_mapping(mapping, default_combined_count_ovs,
                                grep)
    df = df[-grep("Unk", rownames(df)), ]
    save_collapsed_df(df, dataset, annotation, cell_number_filter)
  } 
} else if (dataset == "Greenleaf_brain") {
  lowest_level_annotation_fn = paste("Greenleaf_brain_count_filter", 
                                cell_number_filter, 
                                "combined_count_overlaps.rds", sep="_")
  lowest_level_annotation_fp = paste(root, "Greenleaf_brain_lowest_level_annotation", 
                                     lowest_level_annotation_fn, sep="/")
  lowest_level_combined_count_ovs = readRDS(lowest_level_annotation_fp)
  if (annotation == "Greenleaf_brain_same+as+paper+but+Early+RG+Late+RG-RG_Unk-rm") {
    df = lowest_level_combined_count_ovs
    # pattern --> replacement (collapsing)
    mapping = list(
      c("brain (c0$|c1$|c2$|c5|c6|c7|c13|c14)", "brain GluN"),
      c("brain (c3|c4|c12|c16|c20)", "brain IN"), 
      c("brain c8", "brain nIPC"),
      c("brain (c9|c11)", "brain RG"),
      c("brain c10", "brain mGPC"),
      c("brain c15", "brain OPC/Oligo"),
      c("brain c17", "brain Peric"),
      c("brain c19", "brain MG"),
      c("brain c21", "brain EC"))

    df = collapse_using_mapping(mapping, lowest_level_combined_count_ovs)
    df = df[-grep("c18", rownames(df)), ]
    save_collapsed_df(df, dataset, annotation, cell_number_filter)
  }
}