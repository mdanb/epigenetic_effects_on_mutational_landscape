library(optparse)

####### CUSTOM ANNOTATION OPTIONS ####### 
### Greenleaf pbmc bm
# Greenleaf_pbmc_bm_CD14-mono_CDlike-T_preB+B-B_late+early-no+distinction_Unk-rm

### Greenleaf brain
# Greenleaf_brain_same+as+paper+but+Early+RG+Late+RG-RG_Unk-rm

### Yang Kidney
# Yang_kidney_remove_PT_distinctions

option_list <- list( 
  make_option("--dataset", type="character"),
  make_option("--annotation", type="character")
)

args = parse_args(OptionParser(option_list=option_list))

dataset = args$dataset
annotation = args$annotation

save_collapsed_df <- function(df, df_metadata, dataset, annotation) {
  save_dir = paste("../processed_data/count_overlap_data/combined_count_overlaps",
                   annotation, sep="/")
  save_file = paste(dataset, "combined_count_overlaps.rds", sep = "_")
  save_file_metadata = paste(dataset, "combined_count_overlaps_metadata.rds", 
                             sep = "_")
  
  dir.create(save_dir)
  saveRDS(df, paste(save_dir, save_file, sep="/"))
  saveRDS(df_metadata, paste(save_dir, save_file_metadata, sep="/"))
}

collapse_using_mapping <- function(mapping, df, df_metadata) {
  for (pattern_replacement in mapping) {
    print(paste("Replacing", pattern_replacement[1], 
                "with", pattern_replacement[2]))
    idxs = grep(pattern_replacement[1], rownames(df))
    idxs_metadata = grep(pattern_replacement[1],
                         paste(df_metadata[["tissue_name"]],
                               df_metadata[["cell_type"]]))
    collapsed_co = rep(0, length(colnames(df)))
    collapsed_counts = sum(as.numeric(df_metadata[idxs_metadata, "num_cells"]))
    
    for (idx in idxs) {
      row_to_add = df[idx, ]
      collapsed_co = collapsed_co + row_to_add
    }

    df = df[-idxs, ]
    df_metadata = df_metadata[-idxs_metadata, ]
    if (pattern_replacement[[2]] %in% rownames(df)) {
      df[pattern_replacement[[2]], ] = collapsed_co + df[pattern_replacement[[2]], ]
    } else {
      df[pattern_replacement[[2]], ] = collapsed_co
    }
    temp = unlist(strsplit(pattern_replacement[[2]], split=" "))
    tissue = temp[1]
    cell_type = paste(temp[2:length(temp)], collapse=" ")
    # if (cell_type %in% df_metadata[["cell_type"]]) {
    #   df_metadata[]
    # } else {
    #   df_metadata = rbind(df_metadata, c(tissue, cell_type, collapsed_counts))
    # }
    df_metadata = rbind(df_metadata, c(tissue, cell_type, collapsed_counts))
    # df_metadata[pattern_replacement[[2]], ] = collapsed_counts
  }
  return(list(df, df_metadata))
  # return(list(df, df_metadata))
}

root = "../processed_data/count_overlap_data/combined_count_overlaps"
if (dataset == "Greenleaf_pbmc_bm") {
  default_annotation_fn = "Greenleaf_pbmc_bm_combined_count_overlaps.rds"
  default_annotation_fp = paste(root, "default_annotation", 
                                default_annotation_fn, sep="/")
  default_combined_count_ovs = readRDS(default_annotation_fp)
  
  if (annotation == "Greenleaf_pbmc_bm_CD14-mono_CDlike-T_preB+B-B_late+early-no+distinction_Unk-rm") {
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
    save_collapsed_df(df, dataset, annotation)
  } 
} else if (dataset == "Greenleaf_brain") {
  lowest_level_annotation_fn = "Greenleaf_brain_count_filter_combined_count_overlaps.rds"
  lowest_level_annotation_fp = paste(root, "Greenleaf_brain_lowest_level_annotation", 
                                     lowest_level_annotation_fn, sep="/")
  lowest_level_combined_count_ovs = readRDS(lowest_level_annotation_fp)
  if (annotation == "Greenleaf_brain_same+as+paper+but+Early+RG+Late+RG-RG_Unk-rm") {
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
    save_collapsed_df(df, dataset, annotation)
  }
} else if (dataset == "Yang_kidney") {
  default_annotation_fn = "Yang_kidney_combined_count_overlaps.rds"
  default_annotation_fp = paste(root, "default_annotation", 
                                default_annotation_fn, sep="/")
  default_combined_count_ovs = readRDS(default_annotation_fp)
  default_annotation_metadata_fn = "Yang_kidney_combined_count_overlaps_metadata.rds"
  default_annotation_metadata_fp = paste(root, "default_annotation", 
                                         default_annotation_metadata_fn, sep="/")
  
  default_combined_metadata = readRDS(default_annotation_metadata_fp)
  if (annotation == "Yang_kidney_remove_cell_number_distinctions") {
    # pattern --> replacement (collapsing)
    mapping = list(
      c("kidney (PT1|PT2|PT3|PT4|PT5|PT6|PT7)", "kidney PT"),
      c("kidney (CD_ICA1|CD_ICA2)", "kidney ICA"),
      c("kidney (DCT1|DCT2)", "kidney DCT"),
      c("kidney (CD_PC1|CD_PC2)", "kidney CD_PC"))
    l = collapse_using_mapping(mapping, default_combined_count_ovs, 
                                default_combined_metadata)
    df = l[[1]]
    df_metadata = l[[2]]
    save_collapsed_df(df, df_metadata, dataset, annotation)
  }
} else if (dataset == "Bingren") {
  default_annotation_fn = "Bingren_combined_count_overlaps.rds"
  default_annotation_fp = paste(root, "default_annotation", 
                                default_annotation_fn, sep="/")
  default_combined_count_ovs = readRDS(default_annotation_fp)
  cell_types = rownames(default_combined_count_ovs)
  last_chars <- substr(cell_types, nchar(cell_types), nchar(cell_types))
  is_numbered = grepl("\\d+", last_chars)
  original_string_numbered = cell_types[is_numbered]
  no_numbers = gsub(" \\d+", "", original_string_numbered)
  mapping = unname(mapply(c, original_string_numbered, no_numbers, SIMPLIFY = FALSE))
  default_annotation_metadata_fn = "Bingren_combined_count_overlaps_metadata.rds"
  default_annotation_metadata_fp = paste(root, "default_annotation", 
                                         default_annotation_metadata_fn, sep="/")
  default_combined_metadata = readRDS(default_annotation_metadata_fp)
  l = collapse_using_mapping(mapping, default_combined_count_ovs, 
                             default_combined_metadata)
  df = l[[1]]
  df_metadata = l[[2]]
  save_collapsed_df(df, df_metadata, dataset, annotation)
}