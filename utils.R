library(RColorBrewer)

get_per_cancer_mut_data <- function(all_mutations, interval_ranges) {
  cancer_types_mut_data <- list()
  for (cancer_type in colnames(all_mutations)) {
    idx_cancer_type = match(rownames(as.data.frame(interval_ranges)), 
                            rownames(all_mutations[cancer_type]))
    cancer_types_mut_data = append(cancer_types_mut_data,
                                   list(mut[cancer_type][idx_cancer_type, ]))
  }
  return(cancer_types_mut_data)
}

get_mutation_df_all_cancers <- function(mut, interval_ranges) {
  if (!file.exists("/ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/data/processed_data/mut_count_data.csv")) {
    cancer_types_mut_data = get_per_cancer_mut_data(mut, interval_ranges)
    mut_count_data = as.data.frame(do.call(cbind, cancer_types_mut_data))
    colnames(mut_count_data) = colnames(mut)
    rownames(mut_count_data) = names(interval_ranges)
    write.csv(mut_count_data, "/ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/data/processed_data/mut_count_data.csv")
  }
  else {
    mut_count_data = read.csv("/ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/data/processed_data/mut_count_data.csv")
  }
  return(mut_count_data)
}


# get_tissue_name_tsankov <- function(filename) {
#   tissue_name = unlist(strsplit(str_extract(filename, 
#                                             "count_overlaps_.*_patient"), 
#                                 "_"))
#   tissue_name = tissue_name[3:(length(tissue_name) - 1)]
#   tissue_name = paste(tissue_name, collapse="_")
#   # tissue_name = paste(tissue_name, collapse=" ")
#   # tissue_name = str_to_title(tissue_name)
#   return(tissue_name)
# }

get_n_colors <- function(n, seed) {
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, 
                             rownames(qual_col_pals)))
  set.seed(seed)
  cols <- sample(col_vector, n)
  return(cols)
}

load_mutation_data <- function() {
  mut = readRDS('/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/mutation_data/mutations.aggregated.PCAWG.RK.20181215.Rds')
  # get all mutations, without distinctions e.g clonal vs subclonal 
  mut = mut[, 1:37]
  return(mut)
}

load_meso_mutation_data <- function(waddell_sarc_biph, waddell_sarc,
                                    waddell_sarc_tsankov_sarc, 
                                    waddell_sarc_biph_tsankov_sarc_biph) {
  if (waddell_sarc_biph) {
    return(read.csv('/ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/data/processed_data/meso_mut_count_data.csv'))
  } else if (waddell_sarc) {
    return(read.csv('/ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/data/processed_data/mesothelioma_epithelioid_sarcomatoid2_WGS_Waddell.csv'))
  }
  else if (waddell_sarc_tsankov_sarc) {
    return(read.csv('/ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/data/processed_data/mesothelioma_p786_waddell_sarco.csv'))
  }
  else if (waddell_sarc_biph_tsankov_sarc_biph) {
    return(read.csv('/ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/data/processed_data/mesothelioma_epithelioid_sarcomatoid_biphasic_WGS_Waddell_786_846.csv'))
  }
}