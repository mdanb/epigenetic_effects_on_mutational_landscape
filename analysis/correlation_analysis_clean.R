library(RColorBrewer)
library(tidyverse)
library(preprocessCore)
library(gtools)
source("color.R")
library(parallel)
library(gtools)
library(optparse)

option_list <- list( 
  make_option("--fig2_mss", action="store_true", default=FALSE),
  make_option("--fig2_cll", action="store_true", default=FALSE),
  make_option("--fig2_aml", action="store_true", default=FALSE),
  make_option("--fig3_adeno", action="store_true", default=FALSE),
  make_option("--fig3_neuro", action="store_true", default=FALSE),
  make_option("--fig4_oligo", action="store_true", default=FALSE),
  make_option("--fig4_astro", action="store_true", default=FALSE),
  make_option("--fig4_gbm", action="store_true", default=FALSE)
)

helper <- function(metacell, agg_cancer, scatac_df) {
  corrs = mclapply(metacell, function(x) {
    cor(
      colSums(scatac_df[x, ]),
      agg_cancer)}, mc.cores=8
  )
  return(corrs)
}

perform_and_plot_metacell_correlation <- function(metacells,
                                                  agg_df, 
                                                  scatac_df,
                                                  metacell_correlations_fname,
                                                  cells_to_metacorrelation_fname, 
                                                  embedding_fname, save_fig_fname) {
  print("Getting metacell correlations 1...")
  helper_function <- function(metacell, agg_df, scatac_df) {
    helper(metacell, agg_df, scatac_df)
  }
  if (!file.exists(metacell_correlations_path)) {
    metacell_correlations <- unlist(lapply(metacells, helper_function, 
                                           agg_df, scatac_df))
    metacell_correlations_path = paste("../data/processed_data", metacell_correlations_fname,
                                       sep="/")
    saveRDS(metacell_correlations, metacell_correlations_path)
  }
  
  if (!file.exists(cells_to_metacorrelation_path)) {
    metacell_correlations = readRDS(metacell_correlations_path)
    unique_cells = lapply(lapply(metacells, unlist),
                          unique)
    print("Getting metacell correlations 2...")
    
    idxs = mclapply(unlist(unique_cells), function(x) {
      which(unlist(lapply(metacells[[1]], function(l) {x %in% l})))
    }, mc.cores=8)
    cell_metacorrelations = unlist(lapply(idxs, function(idx_list) 
      mean(metacell_correlations[idx_list])))
    
    # cell_metacorrelations = mapply(compute_cell_metacorrelation, 
    #                                unique_cells_per_cell_type[[1]],
    #                                metacells[[1]],
    #                                metacell_correlations)
    cells_to_metacorrelation = data.frame(cell_barcode=unname(unlist(unique_cells)),
                                          cell_metacorrelation=unname(unlist(cell_metacorrelations)))
    cells_to_metacorrelation_path = paste("../data/processed_data", cells_to_metacorrelation_fname,
                                          sep = "/")
    write.csv(cells_to_metacorrelation, cells_to_metacorrelation_path)
  }
  
  cells_to_metacorrelation = read.csv(cells_to_metacorrelation_path,
                                      row.names = 1)
  embedding_path = paste("../data/processed_data", 
                         embedding_fname, sep = "/")
  embedding = read.csv(embedding_path)
  embedding = as_tibble(embedding)
  colnames(embedding) = c("id", "umap1", "umap2")
  colnames(cells_to_metacorrelation)[1] = "id"
  df = inner_join(embedding, cells_to_metacorrelation)
  
  df = df %>% filter(!is.na(cell_metacorrelation))
  colors = material.heat(3)
  p = ggplot(df) +
    geom_point(aes(x = umap1, y = umap2, color = cell_metacorrelation)) +
    scale_color_gradient2(
      low = colors[3],
      mid = colors[2],  # Specify your desired midpoint color here
      high = colors[1],
      midpoint = ((min(df$cell_metacorrelation, na.rm = TRUE) +  
                     max(df$cell_metacorrelation, na.rm = TRUE)) / 2),  # Set the midpoint value
      limits = c(min(df$cell_metacorrelation, na.rm = TRUE), 
                 max(df$cell_metacorrelation, na.rm = TRUE))
    ) +
    theme_minimal() +  # Use a minimal theme as a starting point
    theme(
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      panel.background = element_blank(),  # Remove panel background
      axis.line = element_line(colour = "black"),  # Add axis lines
      plot.background = element_blank()  # Remove plot background if desired
    )
  # ggsave(filename=paste("../figures", paste0(save_fig_fname, ".png"), sep="/"), 
  #        width = 20, height = 18)
  ggsave(filename=paste("../figures", paste0(save_fig_fname, ".pdf"), sep="/"), 
         width = 20, height = 18)
  
}

args = parse_args(OptionParser(option_list=option_list))
fig2_mss = args$fig2_mss
fig2_cll = args$fig2_cll
fig2_aml = args$fig2_aml
fig3_adeno = args$fig3_adeno
fig3_neuro = args$fig3_neuro
fig4_oligo = args$fig4_oligo
fig4_astro = args$fig4_astro
fig4_gbm = args$fig4_gbm

chr_keep = read.csv("../data/processed_data/chr_keep.csv")[["chr"]]
chr_ranges = unlist(read.csv("../data/processed_data/chr_ranges.csv"))

brain = read.csv("../data/processed_data/mutations_with_subtypes/brain.csv", 
                 row.names=1)

# cell_types = names(metacells_per_cell_type)

# compute_cell_metacorrelation <- function(unique_cells, metacells, 
#                                          correlations_per_metacell) {
#   idxs = mclapply(unique_cells, function(x) {
#     which(unlist(lapply(metacells, function(l) {x %in% l})))
#   }, mc.cores=8)
#   return(unlist(lapply(idxs, function(idx_list) 
#     mean(correlations_per_metacell[idx_list]))))
# }

# cors = cor(t(scatac_df_GL_brain), agg_astro)
# colnames(cors) = c("correlation")
# write.csv(cors, "../data/processed_data/brain_per_cell_correlations.csv")


if (fig4_gbm || fig4_astro || fig4_oligo) {
  scatac_df_GL_brain = readRDS("../data/processed_data/count_overlap_data/combined_count_overlaps/Greenleaf_brain_lowest_level_annotation/per_cell_Greenleaf_brain_combined_count_overlaps.rds")
  scatac_df_GL_brain = scatac_df_GL_brain[, chr_keep]
  scatac_df_GL_brain = scatac_df_GL_brain[, mixedsort(chr_keep)]
  rownames(scatac_df_GL_brain) = paste0(rownames(scatac_df_GL_brain), "-1")
  load("../data/processed_data/Greenleaf_brain_cell_type_independent_nfrags_filter_1_k_500_knnIteration_10000_metacells.Rdata")
  
  metacells = KNN
  
  if (fig4_astro) {
    print("Astrocytoma")
    astro = brain[brain[["subtype"]] == "Astrocytoma", ]
    astro = astro[, 2:2129]
    agg_astro=colSums(astro)
    agg_astro=data.frame(agg_astro[mixedsort(names(agg_astro))])
    
    perform_and_plot_metacell_correlation(metacells, agg_astro, 
                                          scatac_df_GL_brain,
                                          metacell_correlations_fname="astro_nfrags_1_500k_n_100_metacell_correlations_per_cell_type.rds",
                                          cells_to_metacorrelation_fname="astro_nfrags_1_500k_cell_metacorrelations.csv", 
                                          embedding_fname="Greenleaf_brain_nfrags_filter_1_embedding.csv", 
                                          save_fig_fname="astro")
  }
  
  if (fig4_gbm) {
    print("GBM")
    gbm = brain[brain[["subtype"]] == "GBM", ]
    gbm = gbm[, 2:2129]
    agg_gbm=colSums(gbm)
    agg_gbm=data.frame(agg_gbm[mixedsort(names(agg_gbm))])
    
    perform_and_plot_metacell_correlation(metacells, agg_gbm, 
                                          scatac_df_GL_brain,
                                          metacell_correlations_fname="gbm_nfrags_1_500k_n_100_metacell_correlations_per_cell_type.rds",
                                          cells_to_metacorrelation_fname="gbm_nfrags_1_500k_cell_metacorrelations.csv", 
                                          embedding_fname="Greenleaf_brain_nfrags_filter_1_embedding.csv", 
                                          save_fig_fname="gbm")
  }
  
  if (fig4_oligo) {
    print("Oligo")
    oligo = brain[brain[["subtype"]] == "Oligo", ]
    oligo = oligo[, 2:2129]
    agg_oligo=colSums(oligo)
    agg_oligo=data.frame(agg_oligo[mixedsort(names(agg_oligo))])
    
    perform_and_plot_metacell_correlation(metacells, 
                                          agg_oligo, 
                                          scatac_df_GL_brain,
                                          metacell_correlations_fname="oligo_nfrags_1_500k_n_100_metacell_correlations_per_cell_type.rds",
                                          cells_to_metacorrelation_fname="oligo_nfrags_1_500k_cell_metacorrelations.csv", 
                                          embedding_fname="Greenleaf_brain_nfrags_filter_1_embedding.csv", 
                                          save_fig_fname="oligo")
  }
}

if (fig2_mss) {
    print("MSS")
    colon = read.csv("../data/processed_data/mutations_with_subtypes/all_colorectal.csv")
    colon = colon[, chr_keep]
    agg_colon=colSums(colon)
    agg_colon=data.frame(agg_colon[mixedsort(names(agg_colon))])
    
    scatac_df_GL_colon = readRDS("../data/processed_data/count_overlap_data/combined_count_overlaps/default_annotation/per_cell_Greenleaf_colon_combined_count_overlaps.rds")
    scatac_df_GL_colon = scatac_df_GL_colon[, chr_keep]
    scatac_df_GL_colon = scatac_df_GL_colon[, mixedsort(chr_keep)]
    
    # colon = colon[, 2:2129]
    # agg_oligo=colSums(oligo)
    # agg_oligo=data.frame(agg_oligo[mixedsort(names(agg_oligo))])
    load("../data/processed_data/Greenleaf_colon_cell_type_independent_nfrags_filter_10000_k_500_knnIteration_10000_metacells.Rdata")
    metacells = KNN
    
    perform_and_plot_metacell_correlation(metacells, 
                                          agg_colon, 
                                          scatac_df_GL_colon,
                                          metacell_correlations_fname="mss_nfrags_10000_500k_n_100_metacell_correlations_per_cell_type.rds",
                                          cells_to_metacorrelation_fname="mss_nfrags_10000_500k_cell_metacorrelations.csv", 
                                          embedding_fname="Greenleaf_colon_nfrags_filter_10000_embedding.csv", 
                                          save_fig_fname="mss")
}

if (fig3_adeno || fig3_neuro) {
  load("../data/processed_data/Shendure_cell_type_independent_nfrags_filter_1_k_500_knnIteration_10000_metacells.Rdata")
  metacells = KNN
  
  scatac_df_shendure = readRDS("../data/processed_data/count_overlap_data/combined_count_overlaps/default_annotation/per_cell_Shendure_combined_count_overlaps.rds")
  scatac_df_shendure = scatac_df_shendure[, chr_keep]
  scatac_df_shendure = scatac_df_shendure[, mixedsort(chr_keep)]
  
  pancreas = read.csv("../data/processed_data/mutations_with_subtypes/pancreas_all.csv")
  if (fig3_neuro) {
    print("Neuroendocrine")
    neuroendocrine = pancreas[pancreas[["subtype"]] == "Neoroendocrine carcinoma", chr_keep]
    agg_neuroendocrine=colSums(neuroendocrine)
    agg_neuroendocrine=data.frame(agg_neuroendocrine[mixedsort(names(agg_neuroendocrine))])
    
    perform_and_plot_metacell_correlation(metacells, 
                                          agg_neuroendocrine, 
                                          scatac_df_shendure,
                                          metacell_correlations_fname="neuroendocrine_nfrags_1_500k_n_100_metacell_correlations_per_cell_type.rds",
                                          cells_to_metacorrelation_fname="neuroendocrine_nfrags_1_500k_cell_metacorrelations.csv", 
                                          embedding_fname="Shendure_nfrags_filter_1_embedding.csv", 
                                          save_fig_fname="neuroendocrine")
  }
  if (fig3_adeno) {
    print("Panc Adeno")
    panc_adenoca = pancreas[pancreas[["subtype"]] != "Neoroendocrine carcinoma", chr_keep]
    agg_panc_adenoca=colSums(panc_adenoca)
    agg_panc_adenoca=data.frame(agg_panc_adenoca[mixedsort(names(agg_panc_adenoca))])
    
    perform_and_plot_metacell_correlation(metacells, 
                                          agg_panc_adenoca, 
                                          scatac_df_shendure,
                                          metacell_correlations_fname="panc_adenoca_nfrags_1_500k_n_100_metacell_correlations_per_cell_type.rds",
                                          cells_to_metacorrelation_fname="panc_adenoca_nfrags_1_500k_cell_metacorrelations.csv", 
                                          embedding_fname="Shendure_nfrags_filter_1_embedding.csv", 
                                          save_fig_fname="panc_adenoca")
  }
}



# metacell_correlations = lapply(lapply(metacells_per_cell_type, helper, agg_astro), 
#                                unlist)
# saveRDS(metacell_correlations,
#         "../data/processed_data/astro_nfrags_1_500k_n_100_metacell_correlations_per_cell_type.rds")
# 
# metacell_correlations = readRDS("../data/processed_data/astro_nfrags_1_500k_n_100_metacell_correlations_per_cell_type.rds")
# unique_cells_per_cell_type = lapply(lapply(metacells_per_cell_type, unlist),
#                                     unique)
# cell_metacorrelations = mapply(compute_cell_metacorrelation, 
#                                unique_cells_per_cell_type,
#                                metacells_per_cell_type,
#                                metacell_correlations)
# 
# cells_to_metacorrelation = data.frame(cell_barcode=unname(unlist(unique_cells_per_cell_type)),
#                                       cell_metacorrelation=unname(unlist(cell_metacorrelations)))
# write.csv(cells_to_metacorrelation, "../data/processed_data/astro_nfrags_1_500k_cell_metacorrelations.csv")
# cells_to_metacorrelation = read.csv("../data/processed_data/astro_nfrags_1_500k_cell_metacorrelations.csv",
#                                     row.names = 1)
# 
# embedding = read.csv("../data/processed_data/Greenleaf_brain_nfrags_filter_1_embedding.csv")
# embedding = as_tibble(embedding)
# colnames(embedding) = c("id", "umap1", "umap2")
# colnames(cells_to_metacorrelation)[1] = "id"
# df = inner_join(embedding, cells_to_metacorrelation)
# 
# df = df %>% filter(!is.na(cell_metacorrelation))
# colors = material.heat(3)
# p = ggplot(df) +
#   geom_point(aes(x = umap1, y = umap2, color = cell_metacorrelation)) +
#   scale_color_gradient2(
#     low = colors[3],
#     mid = colors[2],  # Specify your desired midpoint color here
#     high = colors[1],
#     midpoint = ((min(df$cell_metacorrelation, na.rm = TRUE) +  
#                    max(df$cell_metacorrelation, na.rm = TRUE)) / 2),  # Set the midpoint value
#     limits = c(min(df$cell_metacorrelation, na.rm = TRUE), 
#                max(df$cell_metacorrelation, na.rm = TRUE))
#   ) +
#   theme_minimal() +  # Use a minimal theme as a starting point
#   theme(
#     panel.grid.major = element_blank(),  # Remove major grid lines
#     panel.grid.minor = element_blank(),  # Remove minor grid lines
#     panel.background = element_blank(),  # Remove panel background
#     axis.line = element_line(colour = "black"),  # Add axis lines
#     plot.background = element_blank()  # Remove plot background if desired
#   )
# 
# ggsave(filename="astro.png", 
#        width = 20, height = 18)

