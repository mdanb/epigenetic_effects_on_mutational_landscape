library(exomeCopy)
library(data.table)
library(stringr)
library(rtracklayer)                                                            

add_cell_barcodes_to_metadata <- function(i, metadatas, barcodes) {
  metadatas[[i]]["cell_barcode"] = barcodes[i]
  return(metadatas[[i]])
}

count_fragments_per_cell <- function(fragments) {
  return(table(fragments$name))
}

filter_metadata_by_fragment_count <- function(metadata, min_frag_counts) {
  if ((tail(metadata["frag_counts_cumsum"], n=1) %>% pull) < min_frag_counts) {
    return(NULL)
  }
  metadata_cut_index_minus_one = sum(metadata["frag_counts_cumsum"] < 
                                       min_frag_counts)
  metadata_cut_index = metadata_cut_index_minus_one + 1
  metadata <- metadata[1:metadata_cut_index, ]
  return(metadata)
}

filter_metadata_by_sample_name <- function(sample_name, metadata) {
  return(metadata[metadata$sample %like% sample_name, ])
}

filter_samples_to_contain_only_cells_in_metadata <- function(i,
                                                             samples,
                                                             sample_barcodes_in_metadatas) {
  samples[[i]] = samples[[i]][samples[[i]]$name %in% 
                                sample_barcodes_in_metadatas[[i]]]
  return(samples[[i]])
}

import_sample <- function(file, dataset) {
  if (dataset == "Bingren") {
    sample = import(paste("/broad", "hptmp", "bgiotti", 
                          "BingRen_scATAC_atlas", 
                          "data", "bed_files", "bingren_scATAC", 
                          "migrated_to_hg19", file, sep="/"), 
                    format="bed")
  }
  else if (dataset == "Shendure") {
    sample = import(paste("/broad", "hptmp", "bgiotti", 
                          "BingRen_scATAC_atlas", 
                          "data", "bed_files", "JShendure_scATAC", 
                          file, sep="/"), format="bed")
  }
  else if (dataset == "Tsankov") {
    sample = import(paste("/broad", "hptmp", "bgiotti", 
                          "BingRen_scATAC_atlas", "data", "bed_files",
                          "Tsankov_scATAC", "migrated_to_hg19",
                          file, sep="/"), format="bed")
  }
  else if (dataset == "Greenleaf_brain") {
    sample = import(paste("/broad", "hptmp", "bgiotti", 
                          "BingRen_scATAC_atlas", "data", "bed_files",
                          "greenleaf_brain_scATAC", "migrated_to_hg19",
                          file, sep="/"), format="bed")
  }
  else if (dataset == "Greenleaf_pbmc_bm") {
    sample = import(paste("/broad", "hptmp", "bgiotti", 
                          "BingRen_scATAC_atlas", "data", "bed_files",
                          "greenleaf_pbmc_bm_scATAC", "migrated_to_hg19",
                          file, sep="/"), format="bed")
  }
  else if (dataset == "Yang_kidney") {
    sample = import(paste("/broad", "hptmp", "bgiotti", 
                          "BingRen_scATAC_atlas", "data", "bed_files",
                          "yang_kidney_scATAC", 
                          file, sep="/"), format="bed")
  }
  return(sample)
}

remove_extension <- function(path) {
  return(gsub("\\.(tsv|txt)\\.b*gz$", "", path))
}

get_sample_filename <- function(file, dataset) {
  if (dataset == "Bingren") {
    filename = paste(paste("Bingren_count_overlaps", remove_extension(file),
                           sep="_"),
                     "rds", sep=".")
  }
  else if (dataset == "Shendure") {
    filename = paste(paste("Shendure_count_overlaps", remove_extension(file),
                           sep="_"),
                     "rds", sep=".")
  }
  else if (dataset == "Tsankov") {
    filename = paste(remove_extension(file))
    filename = unlist(strsplit(filename, split = "_"))
    filename = paste(filename[1:length(filename) - 1], collapse="_")
    filename = paste("Tsankov_count_overlaps", filename, sep="_")
    filename = paste(filename, "rds", sep=".")
  }
  else if (dataset == "Greenleaf_brain") {
    filename = paste(paste("Greenleaf_brain_count_overlaps", 
                           remove_extension(file), sep="_"), "rds", sep = ".")
  }
  else if (dataset == "Greenleaf_pbmc_bm") {
    filename = paste("Greenleaf_pbmc_bm_count_overlaps", 
                     paste(remove_extension(file),
                     "rds", sep="."), sep="_")
  }
  else if (dataset == "Yang_kidney") {
    filename = paste("Yang_kidney_count_overlaps",
                     paste(remove_extension(file),
                     "rds", sep="."), sep="_")
  }
  return(filename)
}

get_fragments_by_cell_barcode <- function(i, sample_idx, fragments, 
                                          top_barcodes) {
  return(fragments[[sample_idx[i]]][fragments[[sample_idx[i]]]$name %in% 
                                      top_barcodes[i]])
}

get_sample_barcodes_in_metadata <- function(filtered_metadata, dataset) {
  if (dataset == "Bingren") {
    sample_barcodes_in_metadata = get_sample_barcodes_in_metadata_helper(filtered_metadata,
                                                                         "cellID",
                                                                         "\\+")
  }
  else if (dataset == "Shendure") {
    sample_barcodes_in_metadata = filtered_metadata[["cell"]]
  }
  else if (dataset == "Tsankov") {
    sample_barcodes_in_metadata = get_sample_barcodes_in_metadata_helper(filtered_metadata, 
                                                                         "sample",
                                                                         "#")
    sample_barcodes_in_metadata = substr(sample_barcodes_in_metadata, 1, 16)
  }
  else if (dataset == "Greenleaf_pbmc_bm" || dataset == "Yang_kidney") {
    sample_barcodes_in_metadata = 
      substr(filtered_metadata[["barcode"]], 1, 16)
  }
  else if (dataset == "Greenleaf_brain") {
    sample_barcodes_in_metadata = 
      substr(filtered_metadata[["Cell.Barcode"]], 1, 16)
  }
  return(sample_barcodes_in_metadata)
}

get_sample_barcodes_in_metadata_helper <- function(metadata, cellID_col_name,
                                                   separator_char) {
  sample_barcodes_in_metadata = unlist(lapply(metadata[[cellID_col_name]],
                                              function(x) unlist(strsplit(x,
                                                                          separator_char))[2]))
  return(sample_barcodes_in_metadata)
}

get_sample_name <- function(file, dataset) {
  if (dataset == "Bingren") {
    sample_name = get_sample_name_bingren(file)
  }
  else if (dataset == "Shendure") {
    sample_name = get_sample_name_shendure(file)
  }
  else if (dataset == "Tsankov") {
    sample_name = get_sample_name_tsankov(file)
  }
  else if (dataset == "Greenleaf_brain") {
    sample_name = get_sample_name_greenleaf_brain(file)
  }
  else if (dataset == "Greenleaf_pbmc_bm") {
    sample_name =  get_sample_name_greenleaf_pbmc_bm(file)
  }
  else if (dataset == "Yang_kidney") {
    sample_name =  get_sample_name_yang(file)
  }
  return(sample_name)
}

get_sample_name_bingren <- function(file) {
  if (grepl("UMB4540_snATAC_frontal_cortex_rep1", file)) {
    return("Human_brain_1")
  }
  else if (grepl("UMB4540_snATAC_frontal_cortex_rep2", file)) {
    return("Human_brain_2")
  }
  sample_name = str_remove(file, "rep")
  sample_name = str_remove(sample_name, "_fragments.bed.bgz")
  sample_name = unlist(strsplit(sample_name, split="_"))
  sample_name = paste(sample_name[2:length(sample_name)], collapse="_")
  return(sample_name)
}

get_sample_name_shendure <- function(file) {
  sample_name = unlist(strsplit(file, split="\\."))[1]
  if (grepl("cerebrum", file)) {
    sample_name = gsub("cerebrum", "brain", sample_name)
  }
  return(sample_name)
}

get_sample_name_tsankov <- function(file) {
  sample_name = str_remove(file, "_fragments.tsv.bgz")
  return(sample_name)
}

get_sample_name_greenleaf_brain <- function(file) {
  sample_name = str_remove(file, "_fragments.tsv.bgz")
  return(sample_name)
}

get_sample_name_greenleaf_pbmc_bm <- function(file) {
  part_one = unlist(strsplit(file, split="_"))[3]
  part_two = unlist(strsplit(unlist(strsplit(file, split="_"))[4], split="\\."))[1]
  sample_name = paste(part_one, part_two, sep="_")
  return(sample_name)
}

get_sample_name_yang <- function(file) {
  sample_name = str_remove(file, ".fragments.tsv.gz")
  sample_name = str_remove(sample_name, ".fragments.tsv.bgz")
  return(sample_name)
}

get_tissue_name <- function(file, dataset, annotation) {
  if (dataset == "Bingren") {
    tissue_name = get_tissue_name_bingren(file)
  }
  else if (dataset == "Shendure") {
    tissue_name = get_tissue_name_shendure(file)
  }
  else if (dataset == "Tsankov") {
    tissue_name = get_tissue_name_tsankov(file, annotation)
  }
  else if (dataset == "Greenleaf_brain") {
    tissue_name = "brain"
  }
  else if (dataset == "Greenleaf_pbmc_bm") {
    tissue_name = get_tissue_name_greenleaf_pbmc_bm(file)
  }
  else if (dataset == "Yang_kidney") {
    tissue_name = "kidney"
  }
  return(tissue_name)
}

get_tissue_name_bingren <- function(filename) {
  if (grepl("frontal_cortex", filename)) {
    tissue_name = "frontal_cortex"
  }
  else {
    tissue_name = unlist(strsplit(str_extract(filename, 
                                              "count_overlaps_.*[0-9]_.*_SM"), 
                                  "_"))
    tissue_name = paste(tissue_name[4:(length(tissue_name)-1)], 
                        collapse="_")
    # tissue_name = str_to_title(paste(tissue_name[4:(length(tissue_name)-1)], 
    #                                  collapse=" "))
  }
  return(tissue_name)
}

# get_tissue_name_bingren <- function(file) {
#   if (grepl("frontal_cortex", file)) {
#     tissue_name = "frontal_cortex"
#   } 
#   else {
#     tissue_name = str_extract(file, pattern="_.*_SM")              
#     tissue_name = substr(tissue_name, 2, nchar(tissue_name) - 3)    
#   }
#   return(tissue_name)
# }

get_tissue_name_shendure <- function(filename) {
  tissue_name = unlist(strsplit(str_extract(filename, 
                                            "sample_.*[.]fragments"), "_"))[3]
  tissue_name = unlist(strsplit(tissue_name, "[.]"))[1]
  tissue_name = paste(tissue_name, collapse="_")
  # tissue_name = str_to_title(tissue_name)
  return(tissue_name)
}

get_tissue_name_greenleaf_pbmc_bm <- function(filename) {
  if (grepl("PBMC", filename)) {
    tissue_name = "blood"
  }
  else {
    tissue_name = "bonemarrow"
  }
  return(tissue_name)
}

get_tissue_name_tsankov <- function(filename) {
  if (annotation == "default_annotation") {
    if (grepl("RPL", filename)) {
      tissue_name = "distal lung"
    }
    else {
      tissue_name = "proximal lung"
    }
  } else if (annotation == "Tsankov_refined") {
    tissue_name = "lung"
  }
  return(tissue_name)
}

# get_tissue_name_shendure <- function(file) {
#   sample_name = get_sample_name_shendure(file)
#   tissue_name = unlist(strsplit(sample_name, split="_"))[3]
#   return(tissue_name)
# }

migrate_file <- function(sample, chain) {
  seqlevelsStyle(sample) = "UCSC"  # necessary
  sample = unlist(liftOver(sample, ch))
  return(sample)
}
