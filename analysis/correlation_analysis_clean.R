library(RColorBrewer)
library(tidyverse)
library(preprocessCore)
library(gtools)
source("color.R")
library(parallel)

chr_keep = read.csv("../data/processed_data/chr_keep.csv")[["chr"]]
chr_ranges = unlist(read.csv("../data/processed_data/chr_ranges.csv"))

brain = read.csv("../data/processed_data/mutations_with_subtypes/brain.csv", row.names=1)
astro = brain[brain[["subtype"]] == "Astrocytoma", ]
astro = astro[, 2:2129]
agg_astro=colSums(astro)
agg_astro=data.frame(agg_astro[mixedsort(names(agg_astro))])

scatac_df_GL_brain = readRDS("../data/processed_data/count_overlap_data/combined_count_overlaps/Greenleaf_brain_lowest_level_annotation/per_cell_Greenleaf_brain_combined_count_overlaps.rds")
scatac_df_GL_brain = scatac_df_GL_brain[, chr_keep]
scatac_df_GL_brain = scatac_df_GL_brain[, mixedsort(chr_keep)]
rownames(scatac_df_GL_brain) = paste0(rownames(scatac_df_GL_brain), "-1")
cors = cor(t(scatac_df_GL_brain), agg_astro)
colnames(cors) = c("correlation")
write.csv(cors, "../data/processed_data/brain_per_cell_correlations.csv")
load("../data/processed_data/Greenleaf_brain_cell_type_independent_nfrags_filter_1_k_500_knnIteration_10000_metacells.Rdata")

metacells_per_cell_type = KNN
cell_types = names(metacells_per_cell_type)

helper <- function(metacells, agg_cancer, scatac_df) {
  corrs = mclapply(metacells, function(x) {
    cor(
      colSums(scatac_df[x, ]),
      agg_cancer)}, mc.cores=8
  )
  return(corrs)
}

compute_cell_metacorrelation <- function(unique_cells, metacells, 
                                         correlations_per_metacell) {
  idxs = mclapply(unique_cells, function(x) {
    which(unlist(lapply(metacells, function(l) {x %in% l})))
  }, mc.cores=8)
  return(unlist(lapply(idxs, function(idx_list) 
    mean(correlations_per_metacell[idx_list]))))
}

perform_and_plot_metacell_correlation <- function(metacells_per_cell_type,
                                                  agg_df, 
                                                  scatac_df,
                                                  metacell_correlations_fname,
                                                  cells_to_metacorrelation_fname, 
                                                  embedding_fname, save_fig_fname) {
  print("Getting metacell correlations 1...")
  metacell_correlations = lapply(mapply(helper,
                                        metacells_per_cell_type,
                                        agg_df,
                                        scatac_df), 
                                 unlist)
  metacell_correlations_path = paste("../data/processed_data", metacell_correlations_fname,
                                    sep="/")
  saveRDS(metacell_correlations, metacell_correlations_path)
  metacell_correlations = readRDS(metacell_correlations_path)
  unique_cells_per_cell_type = lapply(lapply(metacells_per_cell_type, unlist),
                                      unique)
  print("Getting metacell correlations 2...")
  cell_metacorrelations = mapply(compute_cell_metacorrelation, 
                                 unique_cells_per_cell_type,
                                 metacells_per_cell_type,
                                 metacell_correlations)
  cells_to_metacorrelation = data.frame(cell_barcode=unname(unlist(unique_cells_per_cell_type)),
                                        cell_metacorrelation=unname(unlist(cell_metacorrelations)))
  cells_to_metacorrelation_path = paste("../data/processed_data", cells_to_metacorrelation_fname,
                                        sep = "/")
  write.csv(cells_to_metacorrelation, cells_to_metacorrelation_path)
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
  ggsave(filename=paste("../figures", save_fig_fname, sep="/"), 
         width = 20, height = 18)
}

perform_and_plot_metacell_correlation(metacells_per_cell_type, agg_astro, 
                                      scatac_df_GL_brain,
                                      metacell_correlations_fname="astro_nfrags_1_500k_n_100_metacell_correlations_per_cell_type.rds",
                                      cells_to_metacorrelation_fname="astro_nfrags_1_500k_cell_metacorrelations.csv", 
                                      embedding_fname="Greenleaf_brain_nfrags_filter_1_embedding.csv", 
                                      save_fig_fname="astro.png")

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


gbm = brain[brain[["subtype"]] == "GBM", ]
gbm = gbm[, 2:2129]
agg_gbm=colSums(gbm)
agg_gbm=data.frame(agg_gbm[mixedsort(names(agg_gbm))])

perform_and_plot_metacell_correlation(metacells_per_cell_type, agg_gbm, 
                                      scatac_df_GL_brain,
                                      metacell_correlations_fname="gbm_nfrags_1_500k_n_100_metacell_correlations_per_cell_type.rds",
                                      cells_to_metacorrelation_fname="gbm_nfrags_1_500k_cell_metacorrelations.csv", 
                                      embedding_fname="Greenleaf_brain_nfrags_filter_1_embedding.csv", 
                                      save_fig_fname="gbm.png")
