library(RColorBrewer)

get_tissue_name <- function(filename) {
  if (grepl("frontal_cortex", filename)) {
    tissue_name = "frontal_cortex"
  }
  else {
    tissue_name = unlist(strsplit(str_extract(filename, 
                                              "count_overlaps_.*[0-9]_.*_SM"), 
                                  "_"))
    tissue_name = str_to_title(paste(tissue_name[4:(length(tissue_name)-1)], 
                                     collapse=" "))
  }
  return(tissue_name)
}

get_tissue_name_shendure <- function(filename) {
  tissue_name = unlist(strsplit(str_extract(filename, 
                                            "sample_.*[.]fragments"), "_"))[3]
  tissue_name = unlist(strsplit(tissue_name, "[.]"))[1]
  tissue_name = str_to_title(tissue_name)
  return(tissue_name)
}

get_tissue_name_tsankov <- function(filename) {
  tissue_name = unlist(strsplit(str_extract(filename, 
                                            "count_overlaps_.*_patient"), 
                                "_"))
  tissue_name = tissue_name[3:(length(tissue_name) - 1)]
  tissue_name = paste(tissue_name, collapse=" ")
  tissue_name = str_to_title(tissue_name)
  return(tissue_name)
}

get_n_colors <- function(n, seed) {
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, 
                             rownames(qual_col_pals)))
  set.seed(seed)
  cols <- sample(col_vector, n)
  return(cols)
}