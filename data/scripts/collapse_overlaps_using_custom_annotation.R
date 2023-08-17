library(optparse)
library(dplyr)
library(magrittr)
library(tibble)
####### CUSTOM ANNOTATION OPTIONS ####### 
### Greenleaf pbmc bm
# Greenleaf_pbmc_bm_CD14-mono_CDlike-T_preB+B-B_late+early-no+distinction_Unk-rm

### Greenleaf brain
# Greenleaf_brain_same+as+paper+but+Early+RG+Late+RG-RG_Unk-rm

### Yang Kidney
# Yang_kidney_remove_PT_distinctions

option_list <- list( 
  make_option("--dataset", type="character"),
  make_option("--annotation", type="character"),
  make_option("--which_interval_ranges", type="character")
)

args = parse_args(OptionParser(option_list=option_list))

dataset = args$dataset
annotation = args$annotation
which_interval_ranges = args$which_interval_ranges

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
    idxs = grep(pattern_replacement[1], rownames(df))
    idxs_metadata = grep(pattern_replacement[1], paste(df_metadata[["tissue_name"]], 
                                                       df_metadata[["cell_type"]]))
    collapsed_co = rep(0, length(colnames(df)))
    collapsed_counts = sum(as.numeric(df_metadata[idxs_metadata, "num_cells"]))
    
    for (idx in idxs) {
      row_to_add = df[idx, ]
      collapsed_co = collapsed_co + row_to_add
    }
    
    df = df[-idxs, ]
    df_metadata = df_metadata[-idxs_metadata, ]
    df[pattern_replacement[[2]], ] = collapsed_co
    temp = unlist(strsplit(pattern_replacement[[2]], split=" "))
    tissue = temp[1]
    cell_type = paste(temp[2:length(temp)], collapse=" ")
    df_metadata = rbind(df_metadata, c(tissue, cell_type, collapsed_counts))
  }
  return(list(df, df_metadata))
}

root = "../processed_data/count_overlap_data/combined_count_overlaps"
if (dataset == "Greenleaf_pbmc_bm") {
  default_annotation_fn = "Greenleaf_pbmc_bm_combined_count_overlaps.rds"
  if (which_interval_ranges != "polak") {
    default_annotation_fn = paste("interval_ranges", which_interval_ranges,
                                  default_annotation_fn, sep="_")
  }
  default_annotation_fp = paste(root, "default_annotation", 
                                default_annotation_fn, sep="/")
  default_combined_count_ovs = readRDS(default_annotation_fp)
  default_annotation_metadata_fn = "Greenleaf_pbmc_bm_combined_count_overlaps_metadata.rds"
  default_annotation_metadata_fp = paste(root, "default_annotation", 
                                         default_annotation_metadata_fn, sep="/")
  
  default_combined_metadata = readRDS(default_annotation_metadata_fp)
  
  if (annotation == "Greenleaf_pbmc_bm_CD14-mono_CDlike-T_preB+B-B_late+early-no+distinction_Unk-rm") {
    # pattern --> replacement (collapsing)
    mapping = list(
      c("bonemarrow CD14.Mono", "bonemarrow Monocytes"),
      c("blood CD14\\.Mono", "blood Monocytes"), 
      c("bonemarrow CD", "bonemarrow T"),
      c("blood CD", "blood T"), 
      c("bonemarrow B|bonemarrow Pre\\.B", "bonemarrow B"),
      c("bonemarrow .*Eryth", "bonemarrow Eryth"))
    
    l = collapse_using_mapping(mapping, default_combined_count_ovs,
                                default_combined_metadata)
    df = l[[1]]
    df = df[-grep("Unk", rownames(df)), ]
    df_metadata = l[[2]]
    df_metadata = df_metadata[-grep("Unk", df_metadata[["cell_type"]]), ]
    
    save_collapsed_df(df, df_metadata, dataset, annotation)
  } 
} else if (dataset == "Greenleaf_brain") {
  lowest_level_annotation_fn = "Greenleaf_brain_combined_count_overlaps.rds"
  if (which_interval_ranges != "polak") {
    lowest_level_annotation_fn = paste("interval_ranges", which_interval_ranges,
                                       lowest_level_annotation_fn, sep="_")
  }
  lowest_level_annotation_fp = paste(root, "Greenleaf_brain_lowest_level_annotation", 
                                     lowest_level_annotation_fn, sep="/")
  lowest_level_combined_count_ovs = readRDS(lowest_level_annotation_fp)
  
  lowest_level_annotation_metadata_fn = "Greenleaf_brain_combined_count_overlaps_metadata.rds"
  lowest_level_annotation_metadata_fp = paste(root, "Greenleaf_brain_lowest_level_annotation", 
                                              lowest_level_annotation_metadata_fn, sep="/")
  
  lowest_level_annotation_combined_metadata = readRDS(lowest_level_annotation_metadata_fp)
  
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
    l = collapse_using_mapping(mapping, lowest_level_combined_count_ovs, 
                               lowest_level_annotation_combined_metadata)
    
    df = l[[1]]
    df = df[-grep("c18", rownames(df)), ]
    df_metadata = l[[2]]
    df_metadata = df_metadata[-grep("c18", df_metadata[["cell_type"]]), ]
    save_collapsed_df(df, df_metadata, dataset, annotation)
  }
} else if (dataset == "Yang_kidney") {
  default_annotation_fn = "Yang_kidney_combined_count_overlaps.rds"
  if (which_interval_ranges != "polak") {
    default_annotation_fn = paste("interval_ranges", which_interval_ranges,
                                  default_annotation_fn, sep="_")
  }

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
  if (annotation == "Bingren_remove_same_celltype_indexing") {
    default_annotation_fn = "Bingren_combined_count_overlaps.rds"
    if (which_interval_ranges != "polak") {
      default_annotation_fn = paste("interval_ranges", which_interval_ranges,
                                    default_annotation_fn, sep="_")
    }
    default_annotation_fp = paste(root, "default_annotation", 
                                  default_annotation_fn, sep="/")
    default_combined_count_ovs = readRDS(default_annotation_fp)
    cell_types = gsub(" \\d+", "", rownames(default_combined_count_ovs))
    default_combined_count_ovs = as_tibble(default_combined_count_ovs) %>%
      add_column(cell_types, .before=1)
    default_combined_count_ovs = default_combined_count_ovs %>% 
      group_by(cell_types) %>%
      summarise_all(sum)
    cell_types = default_combined_count_ovs[["cell_types"]]
    default_combined_count_ovs = as.data.frame(default_combined_count_ovs)[, 
                                                                           2:ncol(default_combined_count_ovs)]
    rownames(default_combined_count_ovs) = cell_types
    default_annotation_metadata_fn = "Bingren_combined_count_overlaps_metadata.rds"
    default_annotation_metadata_fp = paste(root, "default_annotation", 
                                           default_annotation_metadata_fn, sep="/")
    default_combined_metadata = readRDS(default_annotation_metadata_fp)
    prev_tissue = default_combined_metadata[["tissue_name"]]
    prev_cell_type = default_combined_metadata[["cell_type"]]
    df_metadata = default_combined_metadata
    df_metadata["cell_type"] = paste(prev_tissue,
                                     prev_cell_type)
    df_metadata = df_metadata[, 2:3]
    df_metadata["cell_type"] = gsub(" \\d+", "", 
                                    df_metadata[["cell_type"]])
    df_metadata = df_metadata %>% 
      group_by(cell_type) %>%
      summarise_all(sum)
    temp = strsplit(df_metadata[["cell_type"]], split = " ")
    tissue = unlist(lapply(temp, "[", 1))
    cell_type = lapply(temp, function(x) x[2:length(x)])
    cell_type = unlist(lapply(cell_type, paste, collapse = " "))
    df_metadata["cell_type"] = cell_type
    df_metadata["tissue"] = tissue
    save_collapsed_df(default_combined_count_ovs, df_metadata, dataset, 
                      annotation)
  }
} else if (dataset == "Shendure") {
  if (annotation == "Shendure_remove_unknown_unsure") {
    default_annotation_fn = "Shendure_combined_count_overlaps.rds"
    if (which_interval_ranges != "polak") {
      default_annotation_fn = paste("interval_ranges", which_interval_ranges,
                                    default_annotation_fn, sep="_")
    }
    default_annotation_fp = paste(root, "default_annotation", 
                                  default_annotation_fn, sep="/")
    default_combined_count_ovs = readRDS(default_annotation_fp)
    combined_count_ovs = default_combined_count_ovs[!grepl("Unknown", 
                                                           rownames(default_combined_count_ovs)), ]
    combined_count_ovs = combined_count_ovs[!grepl("\\?", 
                                                   rownames(combined_count_ovs)), ]
    
    default_annotation_metadata_fn = "Shendure_combined_count_overlaps_metadata.rds"
    default_annotation_metadata_fp = paste(root, "default_annotation", 
                                           default_annotation_metadata_fn, sep="/")
    default_combined_metadata = readRDS(default_annotation_metadata_fp)
    combined_metadata = default_combined_metadata[!grepl("Unknown", 
                                                         default_combined_metadata[["cell_type"]]), ]
    combined_metadata = combined_metadata[!grepl("\\?", 
                                                 combined_metadata[["cell_type"]]), ]
    save_collapsed_df(combined_count_ovs, combined_metadata, 
                      dataset, annotation)
  }
} else if (dataset == "Rawlins_fetal_lung") {
  default_annotation_fn = "Rawlins_fetal_lung_combined_count_overlaps.rds"
  if (which_interval_ranges != "polak") {
    default_annotation_fn = paste("interval_ranges", which_interval_ranges,
                                  default_annotation_fn, sep="_")
  }
  default_annotation_fp = paste(root, "default_annotation", 
                                default_annotation_fn, sep="/")
  default_combined_count_ovs = readRDS(default_annotation_fp)
  default_annotation_metadata_fn = "Rawlins_fetal_lung_combined_count_overlaps_metadata.rds"
  default_annotation_metadata_fp = paste(root, "default_annotation", 
                                         default_annotation_metadata_fn, sep="/")
  
  default_combined_metadata = readRDS(default_annotation_metadata_fp)

  if (annotation == "triangle_cells_combine_meso_and_neuroendocrine_no_schwann") {
    mapping = list(
      c("fetal_lung (Pulmonary NE|GHRL\\+ NE)", "fetal_lung Neuroendocrine"),
      c("fetal_lung (Mid/late meso|Early meso)", "fetal_lung Mesothelium"))
      l = collapse_using_mapping(mapping, default_combined_count_ovs, 
                               default_combined_metadata)
    df = l[[1]]
    df_metadata = l[[2]]
    keep = c("Neuroendocrine", "Neuron", "Early fibro", "Mesothelium",
             "Mid fibro", "Myofibro", "Airway SMC", "Airway fibro",
             "Alveolar fibro", "Interm fibro", "Adventitial fibro", 
             "Vascular SMC", "Pericyte")
    df_metadata = df_metadata %>% filter(cell_type %in% keep)
    df = df[rownames(df) %in% paste("fetal_lung", keep), ]
    save_collapsed_df(df, df_metadata, dataset, annotation)
  } else if (annotation == "triangle_cells_combine_meso_no_schwann") {
    mapping = list(
      c("fetal_lung (Mid/late meso|Early meso)", "fetal_lung Mesothelium"))
    l = collapse_using_mapping(mapping, default_combined_count_ovs, 
                               default_combined_metadata)
    df = l[[1]]
    df_metadata = l[[2]]
    keep = c("Neuroendocrine", "Neuron", "Early fibro", "Mesothelium",
             "Mid fibro", "Myofibro", "Airway SMC", "Airway fibro",
             "Alveolar fibro", "Interm fibro", "Adventitial fibro", 
             "Vascular SMC", "Pericyte")
    df_metadata = df_metadata %>% filter(cell_type %in% keep)
    df = df[rownames(df) %in% paste("fetal_lung", keep), ]
    save_collapsed_df(df, df_metadata, dataset, annotation)
  }
} else if (dataset == "Tsankov") {
    refined_annotation_fn = "Tsankov_combined_count_overlaps.rds"
    if (which_interval_ranges != "polak") {
      refined_annotation_fn = paste("interval_ranges", which_interval_ranges,
                                    refined_annotation_fn, sep="_")
    }
    refined_annotation_fp = paste(root, "Tsankov_refined", 
                                  refined_annotation_fn, sep="/")
    refined_combined_count_ovs = readRDS(refined_annotation_fp)
    refined_annotation_metadata_fn = "Tsankov_combined_count_overlaps_metadata.rds"
    refined_annotation_metadata_fp = paste(root, "Tsankov_refined", 
                                           refined_annotation_metadata_fn, sep="/")
    
    refined_combined_metadata = readRDS(refined_annotation_metadata_fp)
    if (annotation == "remove_basal") {
      df = refined_combined_count_ovs
      df_metadata = refined_combined_metadata
      df = df[-grep("Basal", rownames(df)), ]
      df_metadata = df_metadata[-grep("Basal", df_metadata[["cell_type"]]), ]
      save_collapsed_df(df, df_metadata, dataset, annotation)
    } 
} else if (dataset == "Wang_lung") {
  default_annotation_fn = "Wang_lung_combined_count_overlaps.rds"
  default_annotation_fp = paste(root, "default_annotation", 
                                default_annotation_fn, sep="/")
  default_combined_count_ovs = readRDS(default_annotation_fp)
  default_annotation_metadata_fn = "Wang_lung_combined_count_overlaps_metadata.rds"
  default_annotation_metadata_fp = paste(root, "default_annotation", 
                                         default_annotation_metadata_fn, sep="/")
  
  default_combined_metadata = readRDS(default_annotation_metadata_fp)
  if (annotation == "remove_basal") {
    df = default_combined_count_ovs
    df_metadata = default_combined_metadata
    df = df[-grep("basal", rownames(df)), ]
    df_metadata = df_metadata[-grep("basal", df_metadata[["cell_type"]]), ]
    save_collapsed_df(df, df_metadata, dataset, annotation)
  }
}
