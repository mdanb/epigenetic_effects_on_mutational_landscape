library(exomeCopy)
library(data.table)
library(stringr)

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
  samples[i] = samples[[i]][samples[[i]]$name %in% 
                              sample_barcodes_in_metadatas[[i]]]
  return(samples[i])
}

get_fragments_by_cell_barcode <- function(i, sample_idx, fragments, 
                                          top_barcodes) {
  return(fragments[[sample_idx[i]]][fragments[[sample_idx[i]]]$name %in% 
                                      top_barcodes[i]])
}

get_sample_barcodes_in_metadata <- function(metadata, cellID_col_name,
                                            separator_char) {
  sample_barcodes_in_metadata = unlist(lapply(metadata[[cellID_col_name]],
                                              function(x) unlist(strsplit(x,
                                                           separator_char))[2]))
  return(sample_barcodes_in_metadata)
}

get_sample_name <- function(file) {
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