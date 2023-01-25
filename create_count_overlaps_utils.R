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
  if (dataset == "bingren") {
    sample = import(paste("/broad", "hptmp", "bgiotti", 
                          "BingRen_scATAC_atlas", 
                          "raw_dir", "bed_files", file, sep="/"), 
                    format="bed")
  }
  else if (dataset == "shendure") {
    sample = import(paste("/broad", "hptmp", "bgiotti", 
                          "BingRen_scATAC_atlas", 
                          "raw_dir", "bed_files", "JShendure_scATAC", file, 
                          sep="/"), format="bed")
  }
  else if (dataset == "tsankov") {
    sample = import(paste("/broad", "hptmp", "bgiotti", 
                                 "BingRen_scATAC_atlas", "raw_dir", "bed_files",
                                 "Tsankov_scATAC", file, sep="/"), format="bed")
  }
  else if (dataset == "greenleaf_brain") {
    sample = import(paste("/broad", "hptmp", "bgiotti", 
                                 "BingRen_scATAC_atlas", "raw_dir", "bed_files",
                                 "greenleaf_brain_scATAC", 
                                 file, sep="/"), format="bed")
  }
  return(sample)
}

get_sample_filename <- function(file, dataset) {
  if (dataset == "bingren") {
    filename = paste("Bingren_count_filter", cell_number_filter,  
                     "count_overlaps", paste(file_path_sans_ext(file, TRUE),
                                             "rds", sep="."), sep="_")
  }
  else if (dataset == "shendure") {
    filename = paste("Shendure_count_filter", cell_number_filter,  
                     "count_overlaps", paste(file_path_sans_ext(file, TRUE),
                                             "rds", sep="."), sep="_")
  }
  else if (dataset == "tsankov") {
    filename = paste(file_path_sans_ext(file, TRUE))
    filename = unlist(strsplit(filename, split = "_"))
    filename = paste(filename[1:length(filename) - 1], collapse="_")
    filename = paste("Tsankov_count_filter", cell_number_filter,  
                     "count_overlaps", filename, sep="_")
    filename = paste(filename, "rds", sep=".")
  }
  else if (dataset == "greenleaf_brain") {
    filename = paste("Greenleaf_brain_count_filter", cell_number_filter,  
                     "count_overlaps", paste(file_path_sans_ext(file, TRUE),
                                             "rds", sep="."), sep="_")
  }
  else if (dataset == "greenleaf_pbmc_bm") {
    filename = paste("Greenleaf_pbmc_bm_count_filter", cell_number_filter,  
                     "count_overlaps", paste(file_path_sans_ext(file, TRUE),
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
  if (dataset == "bingren") {
    sample_barcodes_in_metadata = get_sample_barcodes_in_metadata_helper(filtered_metadata,
                                                                  "cellID",
                                                                  "\\+")
  }
  else if (dataset == "shendure") {
    sample_barcodes_in_metadata = filtered_metadata[["cell"]]
  }
  else if (dataset == "tsankov") {
    sample_barcodes_in_metadata = get_sample_barcodes_in_metadata_helper(filtered_metadata, 
                                                                  "X",
                                                                  "#")
    sample_barcodes_in_metadata = substr(sample_barcodes_in_metadata, 1, 16)
  }
  else if (dataset == "greenleaf_blood_bm") {
    print("todo")
  }
  else if (dataset == "greenleaf_brain") {
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
  if (dataset == "bingren") {
    sample_name = get_sample_name_bingren(file)
  }
  else if (dataset == "shendure") {
    sample_name = get_sample_name_shendure(file)
  }
  else if (dataset == "tsankov") {
    sample_name = get_sample_name_tsankov(file)
  }
  else if (dataset == "greenleaf_blood_bm") {
    print("todo")
  }
  else if (dataset == "greenleaf_brain") {
    sample_name = get_sample_name_greenleaf_brain(file)
  }
  return(sample_name)
}

get_sample_name_bingren <- function(file) {
  sample_name = str_remove(file, "rep")
  sample_name = str_remove(sample_name, "_fragments.bed.gz")
  sample_name = unlist(strsplit(sample_name, split="_"))
  sample_name = paste(sample_name[2:length(sample_name)], collapse="_")
  return(sample_name)
}

get_sample_name_shendure <- function(file) {
  sample_name = unlist(strsplit(file, split="\\."))[1]
  return(sample_name)
}

get_sample_name_tsankov <- function(file) {
  sample_name = str_remove(file, "_fragments.tsv")
  return(sample_name)
}

get_sample_name_greenleaf_brain <- function(file) {
  sample_name = str_remove(file, "_fragments.tsv.gz")
  return(sample_name)
}

get_tissue_name <- function(file) {
  if (grepl("frontal_cortex", file)) {
    tissue_name = "frontal_cortex"
  } 
  else {
    tissue_name = str_extract(file, pattern="_.*_SM")              
    tissue_name = substr(tissue_name, 2, nchar(tissue_name) - 3)    
  }
  return(tissue_name)
}

get_tissue_name_shendure <- function(file) {
  sample_name = get_sample_name_shendure(file)
  tissue_name = unlist(strsplit(sample_name, split="_"))[3]
  return(tissue_name)
}

migrate_bed_file_to_hg37 <- function(bed_sample, chain) {
  seqlevelsStyle(bed_sample) = "UCSC"  # necessary
  bed_sample = unlist (liftOver(bed_sample, ch))
  return(bed_sample)
}