library(ComplexHeatmap)
library(RColorBrewer)
library(tidyverse)
library(circlize)
library(edgeR) # for kidney analysis
library(preprocessCore)
library(Seurat)
library(gtools)
source("color.R")

chr_keep = read.csv("../data/processed_data/chr_keep.csv")[["chr"]]
chr_ranges = unlist(read.csv("../data/processed_data/chr_ranges.csv"))

#################
get_filtered_cells <- function(cell_num_filter, annotation, dataset,
                               tissue_col_name=NULL, tissue=NULL) {
  root = "../data/processed_data/count_overlap_data/combined_count_overlaps"
  fn = paste(dataset, "combined_count_overlaps_metadata.rds", sep = "_")
  fp = paste(root, annotation, fn, sep = "/")
  metadata = readRDS(fp)
  if (!is.null(tissue_col_name)) {
    metadata = metadata[metadata[[tissue_col_name]] == tissue, ]
  }
  metadata = metadata[as.numeric(metadata[["num_cells"]]) >= cell_num_filter, ]
  if (dataset == "Bingren") {
    metadata["cell_type"] = paste(metadata[["cell_type"]], "BR")
  }
  else if (dataset == "Shendure") {
    metadata["cell_type"] = paste(metadata[["cell_type"]], "SH")
  }
  else if (dataset == "Greenleaf_brain") {
    metadata["cell_type"] = paste(metadata[["cell_type"]], "GL_Br")
  }
  return(metadata)
}

create_cell_fun = function(corrs = NULL, fs) {
  function(j, i, x, y, w, h, fill) {
    grid::grid.text(sprintf("%.2f", corrs[i, j]), x, y, 
                    gp = grid::gpar(col = "black", fontsize = fs))
  }
}

plot_count_distribution <- function(mutation_df, save_filename) {
  mutation_df <- mutation_df %>%
                    gather(key = "column", value = "value") 
  count_distribution <- mutation_df %>%
    group_by(column, value) %>%
    summarise(count = n(), .groups="keep")
  
  ggplot(count_distribution, aes(x = value, y = count)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    facet_wrap(~column, scales = "free") +
    labs(x = "Number of Mutations", y = "Count", 
         title = "Count Distribution of Mutations per Bin per Patient")
  ggsave(save_filename, height=30, width=30)
}
#################

plot_scatac_vs_mutation <- function(df, x_name, y_name) {
  cor_label = cor(df[[x_name]], df[[y_name]])
  p <- ggplot(df) +
        geom_point(aes(x=.data[[x_name]], y=.data[[y_name]])) +
        ggtitle(paste("Pearson's r = ", round(cor_label, 2)))
  fn = paste(paste(x_name, "vs", y_name, sep="_"), "png", sep=".")
  fp = paste("../figures", fn, sep="/")
  ggsave(fp, height = 10, width = 15)
  print(p)
}

#### Lung, Mutations only ####
lung = read.csv("../data/processed_data/mutations_with_subtypes/all_lung.csv")
lung["subtype_grouped_meso"] = lung["subtype"]
lung[lung["subtype"] == "Not.Otherwise.Specified" | lung["subtype"] == 
       "Epithelioid" | lung["subtype"] == "Sarcomatoid" | 
       lung["subtype"] == "Biphasic", 
     "subtype_grouped_meso"] = "meso"
# lung = lung[1:10, ]
subtype = lung[["subtype_grouped_meso"]]
rownames(lung) = lung[["X"]]
lung = lung[, chr_keep]
lung["donor_id"] = rownames(lung)

sclc_subtypes = read.csv("../data/processed_data/mutations_with_subtypes/SCLC.csv")
sclc_subtypes = sclc_subtypes[c("donor_id", "Subtype")]
sclc_subtypes[sclc_subtypes["Subtype"] == "", "Subtype"] = "SCLC_nos"
colnames(sclc_subtypes) = c("donor_id", "SCLC_subtype")
rownames_lung = rownames(lung)
sclc_subtypes = (lung %>% left_join(sclc_subtypes))["SCLC_subtype"]
sclc_subtypes[is.na(sclc_subtypes), "SCLC_subtype"] = "not_SCLC"
lung = lung[, colnames(lung) != "donor_id"]
subtype_colors <- c(
  "Adeno" = "red",
  "Squamous" = "green",
  "SCLC" = "blue",
  "meso" = "orange"
)

sclc_subtype_colors <- c(
  "not_SCLC" = "red",
  "SCLC-N" = "green",
  "SCLC-A" = "blue",
  "SCLC_nos" = "orange",
  "SCLC-P" = "purple",
  "SCLC-Y" = "black"
)

lung = t(lung)
column_ha <- HeatmapAnnotation(subtype = subtype, 
                               SCLC_subtype=sclc_subtypes[["SCLC_subtype"]],
                               counts = log10(colSums(lung)),
                               col = list(subtype = 
                                            subtype_colors,
                                          SCLC_subtype = 
                                            sclc_subtype_colors))
cor_matrix <- cor(lung, method = "pearson")
cor_distance <- as.dist(1 - cor_matrix)
hc <- hclust(cor_distance, method = "ward.D")
Heatmap(cor_matrix, 
        top_annotation = column_ha,
        name = "correlation",
        cluster_columns = hc,
        cluster_rows = hc,
        show_column_names = F,
        show_row_names = F,
        show_row_dend = F,
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")))

# Heatmap(lung,
#         top_annotation = column_ha,
#         clustering_distance_columns = "pearson",
#         clustering_method_columns = "ward.D",
#         cluster_rows = F,
#         show_column_names = F,
#         show_row_names = F, 
#         show_row_dend = F,
#         column_dend_reorder=F)
#### Brain ####
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

helper <- function(metacells, agg_cancer) {
  corrs = mclapply(metacells, function(x) {
    cor(
      colSums(scatac_df_GL_brain[x, ]),
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

metacell_correlations = lapply(lapply(metacells_per_cell_type, helper, agg_astro), 
                               unlist)
saveRDS(metacell_correlations,
        "../data/processed_data/astro_nfrags_1_500k_n_100_metacell_correlations_per_cell_type.rds")
metacell_correlations = readRDS("../data/processed_data/astro_nfrags_1_500k_metacell_correlations_per_cell_type.rds")
unique_cells_per_cell_type = lapply(lapply(metacells_per_cell_type, unlist),
                                    unique)
cell_metacorrelations = mapply(compute_cell_metacorrelation, 
                               unique_cells_per_cell_type,
                               metacells_per_cell_type,
                               metacell_correlations)

cells_to_metacorrelation = data.frame(cell_barcode=unname(unlist(unique_cells_per_cell_type)),
                                      cell_metacorrelation=unname(unlist(cell_metacorrelations)))
write.csv(cells_to_metacorrelation, "../data/processed_data/astro_nfrags_1_500k_cell_metacorrelations.csv")
cells_to_metacorrelation = read.csv("../data/processed_data/astro_nfrags_1_500k_cell_metacorrelations.csv",
                                    row.names = 1)

embedding = read.csv("../data/processed_data/Greenleaf_brain_nfrags_filter_1_embedding.csv")
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
ggsave(filename="astro.png", 
       width = 20, height = 18)

gbm = brain[brain[["subtype"]] == "GBM", ]
gbm = gbm[, 2:2129]
agg_gbm=colSums(gbm)
agg_gbm=data.frame(agg_gbm[mixedsort(names(agg_gbm))])
metacell_correlations = lapply(lapply(metacells_per_cell_type, helper, agg_gbm), 
                               unlist)
cors = cor(t(scatac_df_GL_brain), agg_gbm)

saveRDS(metacell_correlations,
        "../data/processed_data/gbm_nfrags_1_500k_n_100_metacell_correlations_per_cell_type.rds")
unique_cells_per_cell_type = lapply(lapply(metacells_per_cell_type, unlist),
                                    unique)
cell_metacorrelations = mapply(compute_cell_metacorrelation, 
                               unique_cells_per_cell_type,
                               metacells_per_cell_type,
                               metacell_correlations)
cells_to_metacorrelation = data.frame(cell_barcode=unname(unlist(unique_cells_per_cell_type)),
                                      cell_metacorrelation=unname(unlist(cell_metacorrelations)))
write.csv(cells_to_metacorrelation, "../data/processed_data/gbm_nfrags_1_500k_cell_metacorrelations.csv")
cells_to_metacorrelation = read.csv("../data/processed_data/gbm_nfrags_1_500k_cell_metacorrelations.csv",
                                    row.names = 1)


oligo = brain[brain[["subtype"]] == "Oligo", ]
oligo = oligo[, 2:2129]
agg_oligo=colSums(oligo)
agg_oligo=data.frame(agg_oligo[mixedsort(names(agg_oligo))])
cors = cor(t(scatac_df_GL_brain), agg_oligo)

metacell_correlations = lapply(lapply(metacells_per_cell_type, helper, agg_oligo), 
                               unlist)

# p=ggplot(df) +
#   geom_point(aes(umap1, umap2, color=cell_metacorrelation)) +
#   scale_color_gradient(
#     low = "blue", 
#     high = "red",
#     limits = c(min(df$cell_metacorrelation, na.rm = TRUE), 
#                max(df$cell_metacorrelation, na.rm = TRUE)))
# ggsave(filename="umap_10000_frags_variable_k_n_10_knnIteration_100.png", width = 20, height = 20)
ggsave(filename="umap_cell_type_independent_1_frags_k_100_knnIteration_10000.png", 
       width = 20, height = 20)

#### Colon, Mutations only ####
colon = read.csv("../data/processed_data/mutations_with_subtypes/all_colorectal.csv")
# lung = lung[1:10, ]
subtype = colon[["subtype"]]
cbio_subtype = colon[["cbio_subtype"]]
cbio_subtype[cbio_subtype == ""] = "NA"
cancer_type_detailed = colon[["cancer_type_detailed"]]
rownames(colon) = colon[["donor_id"]]
colon = colon[, chr_keep]

agg_colon=colSums(colon)
agg_colon=data.frame(agg_colon[mixedsort(names(agg_colon))])

scatac_df_GL_colon = readRDS("../data/processed_data/count_overlap_data/combined_count_overlaps/default_annotation/per_cell_Greenleaf_colon_combined_count_overlaps.rds")
scatac_df_GL_colon = scatac_df_GL_colon[, chr_keep]
scatac_df_GL_colon = scatac_df_GL_colon[, mixedsort(chr_keep)]

cors = cor(t(scatac_df_GL_colon), agg_colon)
colnames(cors) = c("correlation")
write.csv(cors, "../data/processed_data/colon_per_cell_correlations.csv")

load("../data/processed_data/colorectal_nfrags_filter_10000_500k_metacells.Rdata")
load("../data/processed_data/colorectal_nfrags_filter_10000_100k_metacells.Rdata")
load("../data/processed_data/colorectal_nfrags_filter_5000_variable_k_metacells.Rdata")
load("../data/processed_data/colorectal_nfrags_filter_10000_variable_k_n_10_knnIteration_100_metacells.Rdata")
load("../data/processed_data/colorectal_cell_type_independent_nfrags_filter_10000_k_100_knnIteration_10000_metacells.Rdata")

metacells_per_cell_type = KNN
cell_types = names(metacells_per_cell_type)

helper <- function(metacells) {
  corrs = mclapply(metacells, function(x) {
    cor(
    colSums(scatac_df_GL_colon[x, ]),
    agg_colon)}, mc.cores=8
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

metacell_correlations = lapply(lapply(metacells_per_cell_type, helper), unlist)
saveRDS(metacell_correlations,
    "../data/processed_data/colorectal_nfrags_10000_variable_k_n_100_metacell_correlations_per_cell_type.rds")
metacell_correlations = readRDS("../data/processed_data/colorectal_nfrags_10000_500k_metacell_correlations_per_cell_type.rds")
unique_cells_per_cell_type = lapply(lapply(metacells_per_cell_type, unlist),
                                    unique)
cell_metacorrelations = mapply(compute_cell_metacorrelation, 
                                 unique_cells_per_cell_type,
                                 metacells_per_cell_type,
                                 metacell_correlations)

cells_to_metacorrelation = data.frame(cell_barcode=unname(unlist(unique_cells_per_cell_type)),
                                      cell_metacorrelation=unname(unlist(cell_metacorrelations)))
write.csv(cells_to_metacorrelation, "../data/processed_data/colorectal_nfrags_10000_variable_k_cell_metacorrelations.csv")
cells_to_metacorrelation = read.csv("../data/processed_data/colorectal_nfrags_10000_variable_k_cell_metacorrelations.csv",
                                    row.names = 1)

embedding = read.csv("../data/processed_data/colorectal_nfrags_filter_10000_embedding.csv")
embedding = as_tibble(embedding)
colnames(embedding) = c("id", "umap1", "umap2")
colnames(cells_to_metacorrelation)[1] = "id"
df = inner_join(embedding, cells_to_metacorrelation)

df = df %>% filter(!is.na(cell_metacorrelation))
p=ggplot(df) +
  geom_point(aes(umap1, umap2, color=cell_metacorrelation)) +
  scale_color_gradient(
    low = "blue", 
    high = "red",
    limits = c(min(df$cell_metacorrelation, na.rm = TRUE), 
               max(df$cell_metacorrelation, na.rm = TRUE)))
# ggsave(filename="umap_10000_frags_variable_k_n_10_knnIteration_100.png", width = 20, height = 20)
ggsave(filename="umap_cell_type_independent_10000_frags_k_100_knnIteration_10000.png", width = 20, height = 20)

cbio_subtype_colors <- c(
  "NA" = "red",
  "GS" = "green",
  "CIN" = "blue",
  "POLE" = "orange",
  "MSI" = "purple"
)

subtype_colors <- c(
  "colorectal_adeno_mucinous" = "red",
  "colorectal_adeno_nos" = "green",
  "colorectal_nos" = "blue"
)

cancer_type_detailed_colors = c(
  "Colon Adenocarcinoma" = "red",
  "Mucinous Adenocarcinoma of the Colon and Rectum" = "blue",
  "Rectal Adenocarcinoma" = "green"
)

colon = t(colon)


column_ha <- HeatmapAnnotation(subtype = subtype, 
                               cbio_subtype = cbio_subtype,
                               cancer_type_detailed = cancer_type_detailed,
                               counts = log10(colSums(colon)),
                               col = list(subtype = 
                                            subtype_colors,
                                          cbio_subtype = 
                                            cbio_subtype_colors,
                                          cancer_type_detailed = 
                                            cancer_type_detailed_colors))

# ,
# annotation_legend_param = list(
#   subtype = list(labels = names(subtype_colors), 
#                  title = "Subtype"),
#   cbio_subtype = list(labels = names(cbio_subtype_colors), 
#                       title = "cBio Subtype"),
#   cancer_type_detailed = list(labels = names(cancer_type_detailed_colors), 
#                               title = "Cancer Type Detailed")
# ))


cor_matrix <- cor(colon, method = "pearson")
cor_distance <- as.dist(1 - cor_matrix)
hc <- hclust(cor_distance, method = "ward.D")

# lgd2 = Legend(col_fun = col_fun, title = "legend2", at = c(0, 0.25, 0.5, 0.75, 1))
# lgd3 = Legend(labels = month.name[1:3], legend_gp = gpar(fill = 7:9), title = "legend3")

Heatmap(cor_matrix, 
        top_annotation = column_ha,
        name = "correlation",
        cluster_columns = hc,
        cluster_rows = hc,
        show_column_names = F,
        show_row_names = F,
        show_row_dend = F,
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")))

        

#################
#################

scatac_df_br = t(readRDS("../data/processed_data/count_overlap_data/combined_count_overlaps/Bingren_remove_same_celltype_indexing/Bingren_combined_count_overlaps.rds"))
scatac_df_sh = t(readRDS("../data/processed_data/count_overlap_data/combined_count_overlaps/Shendure_remove_unknown_unsure/Shendure_combined_count_overlaps.rds"))
colnames(scatac_df_sh) = paste(colnames(scatac_df_sh), "SH")
colnames(scatac_df_br) = paste(colnames(scatac_df_br), "BR")
scatac_df = cbind(scatac_df_sh, scatac_df_br)
scatac_df = scatac_df[rownames(scatac_df) %in% chr_keep, ]

#### Meso ####
scatac_df_lung = t(readRDS("../data/processed_data/count_overlap_data/combined_count_overlaps/Tsankov_refined/Tsankov_combined_count_overlaps.rds"))
scatac_df_lung = scatac_df_lung[rownames(scatac_df_lung) %in% chr_keep, ]

mutations_df_meso = read.csv("../data/processed_data/mesothelioma.csv",
                             row.names = 1)

mutations_df_meso = mutations_df_meso[rownames(mutations_df_meso) %in% chr_keep, ]

scatter_plot_data = cbind(scatac_df_lung[, c("lung Fibroblast.WT1+", 
                                         "lung Mesothelium")], 
                          mutations_df_meso[, c("blum_top_10_perc",
                                                "blum_bottom_10_perc",
                                                "blum_50_top_10_perc",
                                                "blum_50_bottom_10_perc",
                                                "nmf_top_10_perc",
                                                "nmf_bottom_10_perc")])
colnames(scatter_plot_data) = gsub("\\s", "_", colnames(scatter_plot_data))
colnames(scatter_plot_data) = gsub("\\+", "pos", colnames(scatter_plot_data))

rankings = c("nmf_bottom_10_perc", "nmf_top_10_perc",
             "blum_bottom_10_perc", "blum_top_10_perc",
             "blum_50_bottom_10_perc", "blum_50_top_10_perc")
cell_types = c("lung_Fibroblast.WT1pos", "lung_Mesothelium")
for (ranking in rankings) {
  for (cell_type in cell_types) {
    plot_scatac_vs_mutation(scatter_plot_data, cell_type, 
                            ranking)
  }
}

# scatac_df_lung = scatac_df_lung[, !(colnames(scatac_df_lung) == "lung Myeloid")]
corrs_pearson = cor(scatac_df_lung, mutations_df_meso)

Heatmap(scale(corrs_pearson), 
        col = RColorBrewer::brewer.pal(9, "RdBu"),
        column_names_gp = grid::gpar(fontsize = 8),
        row_names_gp = grid::gpar(fontsize = 8),
        cell_fun = create_cell_fun(corrs = corrs_pearson, fs=3))

meso_paper = corrs_pearson[, c("mesomics_top_r1_90th_perc", 
                               "mesomics_bottom_r1_10th_perc")]
Heatmap(scale(meso_paper), 
        col = RColorBrewer::brewer.pal(9, "RdBu"),
        column_names_gp = grid::gpar(fontsize = 8),
        row_names_gp = grid::gpar(fontsize = 8),
        cluster_columns = F,
        name = "Pearson's r, \nscaled per column")

df_top_fibroblast = data.frame(muts = mutations_df_meso["mesomics_top_r1_90th_perc"],
                                scatac = scatac_df_lung[, "distal lung Mesothelium"])
df_bottom_mesothelium = data.frame(muts = mutations_df_meso["mesomics_bottom_r1_10th_perc"],
                                scatac = scatac_df_lung[, "distal lung"])
df_top_fibroblast = data.frame(muts = mutations_df_meso["mesomics_top_r1_90th_perc"],
                                scatac = scatac_df_lung)
df_bottom_mesothelium = data.frame(muts = mutations_df_meso["mesomics_bottom_r1_10th_perc"],
                                   scatac = scatac_df_lung)

mutations_df_meso_biphasic = read.csv("../data/processed_data/Mesothelioma_MMB_MutCountperBin.txt",
                                      sep="\t")[4:28]
mutations_df_meso_epithelioid = read.csv("../data/processed_data/Mesothelioma_MME_MutCountperBin.txt",
                                      sep="\t")[4:80]
mutations_df_meso_sarcomatoid = read.csv("../data/processed_data/Mesothelioma_MMS_MutCountperBin.txt",
                                         sep="\t")[4:16]

meso_individuals = cbind(mutations_df_meso_biphasic, 
                         mutations_df_meso_epithelioid,
                         mutations_df_meso_sarcomatoid)

subtypes = c(rep("Biphasic", ncol(mutations_df_meso_biphasic)), 
                    rep("Epithelioid", ncol(mutations_df_meso_epithelioid)), 
                    rep("Sarcomatoid", ncol(mutations_df_meso_sarcomatoid)))
subtypes_df <- data.frame(subtype = subtypes)
rownames(subtypes_df) <- colnames(meso_individuals)

rownames(meso_individuals) = chr_ranges
meso_individuals = meso_individuals[rownames(meso_individuals) %in% chr_keep, ]
scatac_df_lung = t(readRDS("../data/processed_data/count_overlap_data/combined_count_overlaps/Tsankov_separate_fibroblasts/Tsankov_count_filter_1_combined_count_overlaps.rds"))
scatac_df_lung = scatac_df_lung[rownames(scatac_df_lung) %in% chr_keep, ]

subtype_colors <- c(
  "Biphasic" = "red",
  "Epithelioid" = "green",
  "Sarcomatoid" = "blue"
)

ha <- HeatmapAnnotation(subtype = subtypes, col = list(subtype = subtype_colors))

lung_cells_to_keep = c("proximal lung Basal", 
                        "proximal lung Ciliated",
                        "proximal lung Ionocytes",
                        "proximal lung Neuroendocrine",
                        "proximal lung Secretory", 
                        "proximal lung Tuft.like",
                        "distal lung AT1",
                        "distal lung AT2",
                        "distal lung Ciliated",
                        "distal lung Secretory",
                        "distal lung Fibroblasts_C12",
                        "distal lung Fibroblasts_C14",
                        "distal lung Mesothelium")
scatac_df_lung = scatac_df_lung[, colnames(scatac_df_lung) %in% lung_cells_to_keep]
corrs_pearson = cor(scatac_df_lung, meso_individuals)
Heatmap(scale(corrs_pearson), 
        col = RColorBrewer::brewer.pal(9, "RdBu"),
        column_names_gp = grid::gpar(fontsize = 4),
        row_names_gp = grid::gpar(fontsize = 8),
        top_annotation = ha)

top_predictors = c("distal lung Fibroblasts_C14", 
                   "proximal lung Ionocytes",
                   "proximal lung Neuroendocrine",
                   "distal lung AT2",
                   "proximal lung Ciliated",
                   "distal lung Mesothelium",
                   "distal lung SmoothMuscle")

scatac_df_lung = t(readRDS("../data/processed_data/count_overlap_data/combined_count_overlaps/Tsankov_separate_fibroblasts/Tsankov_count_filter_1_combined_count_overlaps.rds"))
scatac_df_lung = scatac_df_lung[rownames(scatac_df_lung) %in% chr_keep, ]
scatac_df_lung = scatac_df_lung[, colnames(scatac_df_lung) %in% top_predictors]
corrs_pearson = cor(scatac_df_lung, meso_individuals)
Heatmap(scale(corrs_pearson), 
        col = RColorBrewer::brewer.pal(9, "RdBu"),
        column_names_gp = grid::gpar(fontsize = 4),
        row_names_gp = grid::gpar(fontsize = 8),
        top_annotation = ha)

plot_count_distribution(meso_individuals, "meso_counts.pdf")
# scatac_df = scatac_df[, !(colnames(scatac_df) == "distal lung Immune")]
# corrs_pearson = cor(scatac_df, meso_individuals)
# 
# Heatmap(scale(corrs_pearson), 
#         col = RColorBrewer::brewer.pal(9, "RdBu"),
#         column_names_gp = grid::gpar(fontsize = 4),
#         row_names_gp = grid::gpar(fontsize = 8),
#         top_annotation = ha)

# ordered_cor_matrix <- cor_matrix[rev(hc$order), rev(hc$order)]
# column_ha <- HeatmapAnnotation(subtype = subtype[rev(hc$order)],
#                                counts = log10(colSums(lung))[rev(hc$order)],
#                                col = list(subtype = subtype_colors))
# 
# # pdf("lung_heatmap.pdf", width = 8, height = 8)
# Heatmap(ordered_cor_matrix, 
#         top_annotation = column_ha,
#         name = "correlation",
#         cluster_columns = F,
#         cluster_rows = F,
#         show_column_names = F,
#         show_row_names = F,
#         col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")))

# dev.off()


#### Lung-AdenoCA ####
lung_cells_to_keep = c("proximal lung Basal", 
                       "proximal lung Ciliated",
                       "proximal lung Ionocytes",
                       "proximal lung Neuroendocrine",
                       "proximal lung Secretory", 
                       "proximal lung Tuft.like",
                       "distal lung AT1",
                       "distal lung AT2",
                       "distal lung Ciliated",
                       "distal lung Secretory")

scatac_df_lung = t(readRDS("../data/processed_data/count_overlap_data/combined_count_overlaps/Tsankov_refined/Tsankov_combined_count_overlaps.rds"))
scatac_df_lung = scatac_df_lung[rownames(scatac_df_lung) %in% chr_keep, ]
lung_adenoca = read.csv("../data/mutation_data/Lung-AdenoCA.txt",
                        sep="\t")
lung_adenoca = lung_adenoca[, 4:ncol(lung_adenoca)]
lung_adenoca = lung_adenoca[chr_ranges %in% chr_keep, ]

lung_squamous = read.csv("../data/mutation_data/Lung-SCC.txt",
                         sep="\t")
lung_squamous = lung_squamous[, 4:ncol(lung_squamous)]
lung_squamous = lung_squamous[chr_ranges %in% chr_keep, ]
lung_subtypes = read.csv("../data/TCGA_ICGCDonor_Subtype.txt", sep = "\t",
                         header = 0)
rownames(lung_subtypes) = lung_subtypes[["V2"]]
adeno_subtypes = lung_subtypes[colnames(lung_adenoca), "V3"]
squamous_subtypes = lung_subtypes[colnames(lung_squamous), "V3"]

subtype_colors <- c(
  "Immune" = "red",
  "S100.MAFK" = "green",
  "Secretory.AT2" = "blue",
  "Mucinous" = "black",
  "Secretory.Ciliated.AT" = "purple",
  "FattyAcid.AT" = "pink",
  "MTORC1.SOX2" = "cyan",
  "Neuroendocrine" = "orange", 
  "Immune.EMT" = "grey",
  "Mesenchymal" = "yellow" 
)

# scatac_df_lung = scatac_df_lung[, !(colnames(scatac_df_lung) == "distal lung Immune")]

scatac_df_lung = scatac_df_lung[, colnames(scatac_df_lung) %in% 
                                  lung_cells_to_keep]
corrs_pearson = cor(scatac_df_lung, lung_adenoca)
ha <- HeatmapAnnotation(subtype = adeno_subtypes, col = list(subtype = subtype_colors))
Heatmap(scale(corrs_pearson), 
        col = RColorBrewer::brewer.pal(9, "RdBu"),
        column_names_gp = grid::gpar(fontsize = 2),
        row_names_gp = grid::gpar(fontsize = 8),
        top_annotation = ha,
        cell_fun = create_cell_fun(corrs = corrs_pearson, fs=3.7))

scatac_df_lung_no_AT1 = as.data.frame(scatac_df_lung) %>% select(-`distal lung AT1`)
corrs_pearson = cor(scatac_df_lung_no_AT1, lung_adenoca)
ha <- HeatmapAnnotation(subtype = adeno_subtypes, col = list(subtype = subtype_colors))
Heatmap(scale(corrs_pearson), 
        col = RColorBrewer::brewer.pal(9, "RdBu"),
        column_names_gp = grid::gpar(fontsize = 2),
        row_names_gp = grid::gpar(fontsize = 8),
        top_annotation = ha,
        cell_fun = create_cell_fun(corrs = corrs_pearson, fs=3.7))

scatac_df_top_5 = scatac_df[, c("lung Ciliated epithelial cells SH",
                                "colon_transverse Small Intestinal Enterocyte BR")]
scatac_df_top_5 = cbind(scatac_df_top_5, scatac_df_lung[, c("distal lung AT2",
                                                            "distal lung Secretory",
                                                            "proximal lung Ciliated")])
corrs_pearson = cor(scatac_df_top_5, lung_adenoca)
ha <- HeatmapAnnotation(subtype = adeno_subtypes, col = list(subtype = subtype_colors))
Heatmap(scale(corrs_pearson), 
        col = RColorBrewer::brewer.pal(9, "RdBu"),
        column_names_gp = grid::gpar(fontsize = 2),
        row_names_gp = grid::gpar(fontsize = 8),
        top_annotation = ha,
        cell_fun = create_cell_fun(corrs = corrs_pearson, fs=3.7))

scatac_df_top_5 = as.data.frame(scatac_df_top_5) %>% select(-`colon_transverse Small Intestinal Enterocyte BR`)
corrs_pearson = cor(scatac_df_top_5, lung_adenoca)
ha <- HeatmapAnnotation(subtype = adeno_subtypes, col = list(subtype = subtype_colors))
Heatmap(scale(corrs_pearson), 
        col = RColorBrewer::brewer.pal(9, "RdBu"),
        column_names_gp = grid::gpar(fontsize = 2),
        row_names_gp = grid::gpar(fontsize = 8),
        top_annotation = ha,
        cell_fun = create_cell_fun(corrs = corrs_pearson, fs=3.7))

plot_count_distribution(lung_adenoca, "lung_adeno_counts.pdf")

subtype_colors <- c(
  "AT.NKX2-1" = "red",
  "mTOR.SOX2" = "green",
  "COL.SOX4" = "blue",
  "KRT13,4.FOXA1" = "black",
  "CellCycle" = "purple",
  "S100A.HES2" = "cyan"
)

corrs_pearson = cor(scatac_df_lung, lung_squamous)
ha <- HeatmapAnnotation(subtype = squamous_subtypes, col = list(subtype = subtype_colors))
Heatmap(scale(corrs_pearson), 
        col = RColorBrewer::brewer.pal(9, "RdBu"),
        column_names_gp = grid::gpar(fontsize = 2),
        row_names_gp = grid::gpar(fontsize = 8),
        top_annotation = ha,
        cell_fun = create_cell_fun(corrs = corrs_pearson, fs=3.7))

scatac_df_lung = as.data.frame(scatac_df_lung) %>% select(-`distal lung AT1`)
corrs_pearson = cor(scatac_df_lung, lung_squamous)
ha <- HeatmapAnnotation(subtype = squamous_subtypes, col = list(subtype = subtype_colors))
Heatmap(scale(corrs_pearson), 
        col = RColorBrewer::brewer.pal(9, "RdBu"),
        column_names_gp = grid::gpar(fontsize = 2),
        row_names_gp = grid::gpar(fontsize = 8),
        top_annotation = ha,
        cell_fun = create_cell_fun(corrs = corrs_pearson, fs=3.7))

scatac_df_top_5 = scatac_df[, c("lung Bronchiolar and alveolar epithelial cells SH",
                                "lung Ciliated epithelial cells SH")]
scatac_df_top_5 = cbind(scatac_df_top_5, scatac_df_lung[, c("proximal lung Basal",
                                                            "proximal lung Tuft.like",
                                                            "proximal lung Ciliated")])
corrs_pearson = cor(scatac_df_top_5, lung_squamous)
ha <- HeatmapAnnotation(subtype = squamous_subtypes, col = list(subtype = subtype_colors))
Heatmap(scale(corrs_pearson), 
        col = RColorBrewer::brewer.pal(9, "RdBu"),
        column_names_gp = grid::gpar(fontsize = 2),
        row_names_gp = grid::gpar(fontsize = 8),
        top_annotation = ha,
        cell_fun = create_cell_fun(corrs = corrs_pearson, fs=3))
plot_count_distribution(lung_squamous, "lung_squamous_counts.pdf")


#### Pancreas-AdenoCA ####
pancreas_adenoca = read.csv("../data/mutation_data/Panc-AdenoCA.txt",
                            sep="\t")
pancreas_adenoca = pancreas_adenoca[chr_ranges %in% chr_keep, ]
pancreas_adenoca = pancreas_adenoca[, 4:ncol(pancreas_adenoca)]

pancreas_endocrine = read.csv("../data/mutation_data/Panc-Endocrine.txt",
                            sep="\t")
pancreas_endocrine = pancreas_endocrine[chr_ranges %in% chr_keep, ]
pancreas_endocrine = pancreas_endocrine[, 4:ncol(pancreas_endocrine)]

PACA_AU_donor_to_subtype = read.csv("../data/mutation_data/icgc/PACA_AU_donor_to_subtype.csv")
PACA_CA_donor_to_subtype = read.csv("../data/mutation_data/icgc/PACA_CA_donor_to_subtype.csv")
adeno_donor_subtypes = rbind(PACA_AU_donor_to_subtype, PACA_CA_donor_to_subtype)
df = data.frame(adeno_donor_subtypes, row.names = adeno_donor_subtypes[["donor_id"]])
adeno_donor_to_subtype = df[colnames(pancreas_adenoca), "X"]
adeno_donor_to_subtype[adeno_donor_to_subtype == ""] = "NA"
adeno_donor_to_subtype = c(adeno_donor_to_subtype, rep("not_adeno",
                                                       ncol(pancreas_endocrine)))
subtype = c(rep("Adeno", ncol(pancreas_adenoca)),
            rep("Endocrine", ncol(pancreas_endocrine)))

subtype_colors <- c(
  "Endocrine" = "red",
  "Adeno" = "blue"
)

adeno_subtype_colors <- c(
  "Pancreatic ductal carcinoma" = "red",
  "Invasive carcinoma arising in IPMN" = "green",
  "Carcinoma, adenosquamous" = "blue",
  "Adenocarcinoma, mucinous" = "black",
  "Acinar cell carcinoma" = "purple",
  "Adenocarcinoma" = "pink",
  "NA" = "orange",
  "not_adeno" = "cyan"
)
pancreas = cbind(pancreas_adenoca, pancreas_endocrine)
column_ha <- HeatmapAnnotation(subtype = subtype, 
                               adeno_subtype = adeno_donor_to_subtype,
                               counts = log10(colSums(pancreas)),
                               col = list(subtype = 
                                            subtype_colors,
                                          adeno_subtype = adeno_subtype_colors))
cor_matrix <- cor(pancreas, method = "pearson")
cor_distance <- as.dist(1 - cor_matrix)
hc <- hclust(cor_distance, method = "ward.D")
Heatmap(cor_matrix, 
        top_annotation = column_ha,
        name = "correlation",
        cluster_columns = hc,
        cluster_rows = hc,
        show_column_names = F,
        show_row_names = F,
        show_row_dend = F,
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")))



per_patient_density = log(colSums(pancreas_adenoca))
# lowest = per_patient_density == min(per_patient_density)
# pancreas_adenoca = pancreas_adenoca[, !lowest]
subtype_colors <- c(
  "Pancreatic ductal carcinoma" = "red",
  "Invasive carcinoma arising in IPMN" = "green",
  "Carcinoma, adenosquamous" = "blue",
  "Adenocarcinoma, mucinous" = "black",
  "Acinar cell carcinoma" = "purple",
  "Adenocarcinoma" = "pink"
)

dataset_colors = c("PACA-AU" = "brown",
                   "PACA-CA" = "cyan")

PACA_AU_donor_to_subtype = read.csv("../data/mutation_data/icgc/PACA_AU_donor_to_subtype.csv")
PACA_CA_donor_to_subtype = read.csv("../data/mutation_data/icgc/PACA_CA_donor_to_subtype.csv")
donor_subtypes = rbind(PACA_AU_donor_to_subtype, PACA_CA_donor_to_subtype)
df = data.frame(donor_subtypes, row.names = donor_subtypes[["donor_id"]])
donor_to_subtype = df[colnames(pancreas_adenoca), "X"]
dataset_annotation = rep("null", ncol(pancreas_adenoca))
dataset_annotation[colnames(pancreas_adenoca) %in% 
  PACA_AU_donor_to_subtype[["donor_id"]]] = "PACA-AU"
dataset_annotation[colnames(pancreas_adenoca) %in% 
                     PACA_CA_donor_to_subtype[["donor_id"]]] = "PACA-CA"

pancreas_metadata_br = get_filtered_cells(cell_num_filter = 50, 
                                          annotation = "Bingren_remove_same_celltype_indexing", 
                                          dataset ="Bingren", 
                                          tissue_col_name = "tissue", 
                                          tissue = "pancreas")
pancreas_metadata_sh = get_filtered_cells(cell_num_filter = 50, 
                                          annotation = "Shendure_remove_unknown_unsure", 
                                          dataset = "Shendure", 
                                          tissue_col_name = "tissue_name", 
                                          tissue = "pancreas")

pancreas_epithelial_cells = c("pancreas Acinar cells SH", 
                              "pancreas Ductal cells SH",
                              "pancreas Ductal Cell (Pancreatic) BR",
                              "pancreas Pancreatic Acinar Cell BR",
                              "pancreas Pancreatic Beta Cell BR",
                              "pancreas Pancreatic Delta,Gamma cell BR")
scatac_df_panc_epithelial = scatac_df[, pancreas_epithelial_cells]

stomach_metadata_br = get_filtered_cells(cell_num_filter = 50, 
                                          annotation = "Bingren_remove_same_celltype_indexing", 
                                          dataset ="Bingren", 
                                          tissue_col_name = "tissue", 
                                          tissue = "stomach")
stomach_metadata_sh = get_filtered_cells(cell_num_filter = 50, 
                                          annotation = "Shendure_remove_unknown_unsure", 
                                          dataset = "Shendure", 
                                          tissue_col_name = "tissue_name", 
                                          tissue = "stomach")
stomach_goblet = c("stomach Goblet cells SH", 
                   "stomach Foveolar Cell BR")

scatac_df_stomach_mucus = scatac_df[, stomach_goblet]
scatac_df_stomach_panc = cbind(scatac_df_stomach_mucus, scatac_df_panc_epithelial)

corrs_pearson = cor(scatac_df_stomach_panc, pancreas_adenoca)
pancreas_rna_annotation = read.table("../data/Annotation_Pancreatic_CPTAC.tsv",
                                     header=1)
rownames(pancreas_rna_annotation) = pancreas_rna_annotation[["Donor"]]
pancreas_rna_subtypes = pancreas_rna_annotation[colnames(pancreas_adenoca), "Subtype"]

rna_subtype_colors <- c(
  "Pancreatic_Progenitor" = "red",
  "Acinar" = "green",
  "Fibroblast" = "blue",
  "Squamous" = "black",
  "Immune" = "purple",
  "Endocrine" = "pink",
  "Neuroendocrine" = "yellow"
)

per_patient_density = log(colSums(pancreas_adenoca))
continuous_color_palette = colorRamp2(c(min(per_patient_density), 
                                      max(per_patient_density)),
                                      c("white", "red"))
# per_patient_density_colors = continuous_color_palette(per_patient_density)
# combined_annotation = cbind(subtype = donor_to_subtype, dataset = dataset_annotation,
#                             rna_subtype = pancreas_rna_subtypes, 
#                             per_patient_density = per_patient_density)
combined_colors = list(subtype = subtype_colors, dataset = dataset_colors,
                       rna_subtype = rna_subtype_colors, 
                       per_patient_density = continuous_color_palette)
# ha <- HeatmapAnnotation(subtype = donor_to_subtype, col = list(subtype = subtype_colors))

ha <- HeatmapAnnotation(subtype = donor_to_subtype, dataset = dataset_annotation,
                        rna_subtype = pancreas_rna_subtypes, 
                        per_patient_density = per_patient_density, col = combined_colors)
column_order <- order(scale(corrs_pearson)["stomach Goblet cells SH", ])
pdf("goblet.pdf", width=10, height=10)
Heatmap(scale(corrs_pearson), 
        col = RColorBrewer::brewer.pal(9, "RdBu"),
        column_names_gp = grid::gpar(fontsize = 2),
        row_names_gp = grid::gpar(fontsize = 8),
        # top_annotation = ha,
        # column_order = column_order,
        # cell_fun = create_cell_fun(corrs = corrs_pearson, fs=5),
        show_heatmap_legend = F)
dev.off()

mat = t(scale(corrs_pearson))
column_hclust <- hclust(dist(mat))
k <- 2  # number of clusters
col_clusters <- cutree(column_hclust, k)
delta_gamma_cluster = names(col_clusters)[col_clusters == 2]
delta_gamma_cluster = pancreas_adenoca[, delta_gamma_cluster]
delta_gamma_cluster = as.data.frame(rowSums(delta_gamma_cluster))
rownames(delta_gamma_cluster) = chr_keep
colnames(delta_gamma_cluster) = "pancreas_delta_gamma"
write.csv(delta_gamma_cluster, file = "../data/processed_data/hierarchically_clustered_mutations.csv")

column_order <- order(scale(corrs_pearson)["stomach Foveolar Cell BR", ])
pdf("foveolar.pdf", width=45, height=10)
Heatmap(scale(corrs_pearson), 
        col = RColorBrewer::brewer.pal(9, "RdBu"),
        column_names_gp = grid::gpar(fontsize = 2),
        row_names_gp = grid::gpar(fontsize = 8),
        # top_annotation = ha,
        column_order = column_order,
        cell_fun = create_cell_fun(corrs = corrs_pearson, fs=5),
        show_heatmap_legend = F)
dev.off()

Heatmap(scale(corrs_pearson)[grepl("goblet|acinar|ductal", rownames(corrs_pearson),
                                   ignore.case = T), ], 
        col = RColorBrewer::brewer.pal(9, "RdBu"),
        column_names_gp = grid::gpar(fontsize = 2),
        row_names_gp = grid::gpar(fontsize = 8),
        top_annotation = ha,
        # column_order = column_order,
        # cell_fun = create_cell_fun(corrs = corrs_pearson, fs=5),
        show_heatmap_legend = F)

plot_count_distribution(pancreas_adenoca, "panc_adenoca_counts.pdf")

#### ColoRect-AdenoCA ####
dataset = "Bingren"
annotation = "Bingren_remove_same_celltype_indexing"

colon_sigmoid_metadata_br = get_filtered_cells(cell_num_filter = 50, 
                                          annotation = "Bingren_remove_same_celltype_indexing", 
                                          dataset ="Bingren", 
                                          tissue_col_name = "tissue", 
                                          tissue = "colon_sigmoid")
colon_transverse_metadata_br = get_filtered_cells(cell_num_filter = 50, 
                                               annotation = "Bingren_remove_same_celltype_indexing", 
                                               dataset ="Bingren", 
                                               tissue_col_name = "tissue", 
                                               tissue = "colon_transverse")
intestine_metadata_sh = get_filtered_cells(cell_num_filter = 50, 
                                         annotation = "Shendure_remove_unknown_unsure", 
                                         dataset = "Shendure", 
                                         tissue_col_name = "tissue_name", 
                                         tissue = "intestine")
colon_transverse_cells = c("colon_transverse Colon Epithelial Cell BR",
                           "colon_transverse Colonic Goblet Cell BR",
                           "colon_transverse Enterochromaffin Cell BR",
                           "colon_transverse Small Intestinal Enterocyte BR")

colon_transverse_cells = c("colon_transverse Colon Epithelial Cell BR",
                           "colon_transverse Colonic Goblet Cell BR",
                           "colon_transverse Enterochromaffin Cell BR")

intestine_cells = c("intestine Chromaffin cells SH",
                    "intestine Intestinal epithelial cells SH",
                    "intestine Mesothelial cells SH")
stomach_goblet = c("stomach Goblet cells SH", 
                   "stomach Foveolar Cell BR")

colon_intestine_stomach_cells = c(colon_transverse_cells, intestine_cells, stomach_goblet)
colon_intestine_scatac_df = scatac_df[, colon_intestine_stomach_cells]
colon_adenoca = read.csv("../data/mutation_data/ColoRect-AdenoCA.txt", sep="\t")
colon_adenoca = colon_adenoca[chr_ranges %in% chr_keep, ]
colon_adenoca = colon_adenoca[, 4:ncol(colon_adenoca)]

corrs_pearson = cor(colon_intestine_scatac_df, colon_adenoca)
column_order <- order(scale(corrs_pearson)["stomach Foveolar Cell BR", ])
Heatmap(scale(corrs_pearson), 
        col = RColorBrewer::brewer.pal(9, "RdBu"),
        column_names_gp = grid::gpar(fontsize = 2),
        row_names_gp = grid::gpar(fontsize = 8),
        column_order = column_order,
        cell_fun = create_cell_fun(corrs = corrs_pearson, fs=3.7))

column_order <- order(scale(corrs_pearson)["stomach Goblet cells SH", ])
Heatmap(scale(corrs_pearson), 
        col = RColorBrewer::brewer.pal(9, "RdBu"),
        column_names_gp = grid::gpar(fontsize = 2),
        row_names_gp = grid::gpar(fontsize = 8),
        column_order = column_order,
        cell_fun = create_cell_fun(corrs = corrs_pearson, fs=3.7))

column_order <- order(scale(corrs_pearson)["intestine Intestinal epithelial cells SH", ])
Heatmap(scale(corrs_pearson), 
        col = RColorBrewer::brewer.pal(9, "RdBu"),
        column_names_gp = grid::gpar(fontsize = 2),
        row_names_gp = grid::gpar(fontsize = 8),
        column_order = column_order,
        cell_fun = create_cell_fun(corrs = corrs_pearson, fs=3.7))

column_order <- order(scale(corrs_pearson)["colon_transverse Colonic Goblet Cell BR", ])
Heatmap(scale(corrs_pearson), 
        col = RColorBrewer::brewer.pal(9, "RdBu"),
        column_names_gp = grid::gpar(fontsize = 2),
        row_names_gp = grid::gpar(fontsize = 8),
        column_order = column_order,
        cell_fun = create_cell_fun(corrs = corrs_pearson, fs=3.7))

column_order <- order(scale(corrs_pearson)["intestine Intestinal epithelial cells SH", ])
Heatmap(scale(corrs_pearson), 
        col = RColorBrewer::brewer.pal(9, "RdBu"),
        column_names_gp = grid::gpar(fontsize = 2),
        row_names_gp = grid::gpar(fontsize = 8),
        column_order = column_order,
        cell_fun = create_cell_fun(corrs = corrs_pearson, fs=3.7))

plot_count_distribution(colon_adenoca, "colorect_adenoca_counts.pdf")


cancer_samples_atac = read.csv("../data/processed_data/binned_atac/binned_atac_COAD-US.csv",
                                row.names = 1)
cancer_samples_atac = cancer_samples_atac[chr_keep, ]
colon_intestine_scatac_df = colon_intestine_scatac_df[, colnames(colon_intestine_scatac_df) !=
                                                        "colon_transverse Enterochromaffin Cell BR"]
colon_intestine_scatac_df = colon_intestine_scatac_df[, colnames(colon_intestine_scatac_df) !=
                                                        "colon_transverse Colonic Goblet Cell BR"]

corrs_pearson = cor(colon_intestine_scatac_df, cancer_samples_atac)
pdf("colon_atac_corrs.pdf", width=15, height=10)
Heatmap(corrs_pearson, 
        col = RColorBrewer::brewer.pal(9, "RdBu"),
        column_names_gp = grid::gpar(fontsize = 2),
        row_names_gp = grid::gpar(fontsize = 8),
        cell_fun = create_cell_fun(corrs = corrs_pearson, fs=3))
dev.off()
#### Brain ####
scatac_df_gl = t(readRDS("../data/processed_data/count_overlap_data/combined_count_overlaps/Greenleaf_brain_same+as+paper+but+Early+RG+Late+RG-RG_Unk-rm/Greenleaf_brain_combined_count_overlaps.rds"))
colnames(scatac_df_gl) = paste(colnames(scatac_df_gl), "GL_Br")
scatac_df_gl = scatac_df_gl[rownames(scatac_df_gl) %in% chr_keep, ]

gbm = read.csv("../data/mutation_data/CNS-GBM.txt", sep="\t")
gbm = gbm[, 4:ncol(gbm)]
gbm = gbm[chr_ranges %in% chr_keep, ]

astrocytoma = read.csv("../data/mutation_data/CNS-PiloAstro.txt", sep="\t")
astrocytoma = astrocytoma[, 4:ncol(astrocytoma)]
astrocytoma = astrocytoma[chr_ranges %in% chr_keep, ]

oligo = read.csv("../data/mutation_data/CNS-Oligo.txt", sep="\t")
oligo = oligo[, 4:ncol(oligo)]
oligo = oligo[chr_ranges %in% chr_keep, ]

medullo = read.csv("../data/mutation_data/CNS-Medullo.txt", sep="\t")
medullo = medullo[, 4:ncol(medullo)]
medullo = medullo[chr_ranges %in% chr_keep, ]

brain_mutations = cbind(gbm, astrocytoma, oligo, medullo)

subtype = c(rep("GBM", ncol(gbm)),
                     rep("Astrocytoma", ncol(astrocytoma)),
                     rep("Oligo", ncol(oligo)),
                     rep("Medullo", ncol(medullo)))

subtype_colors <- c(
  "GBM" = "red",
  "Astrocytoma" = "green",
  "Oligo" = "blue",
  "Medullo" = "purple"
)

column_ha <- HeatmapAnnotation(subtype = subtype, 
                               counts = log10(colSums(brain_mutations)),
                               col = list(subtype = 
                                            subtype_colors))
cor_matrix <- cor(brain_mutations, method = "pearson")
cor_distance <- as.dist(1 - cor_matrix)
hc <- hclust(cor_distance, method = "ward.D")
Heatmap(cor_matrix, 
        top_annotation = column_ha,
        name = "correlation",
        cluster_columns = hc,
        cluster_rows = hc,
        show_column_names = F,
        show_row_names = F,
        show_row_dend = F,
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")))


column_ha <- HeatmapAnnotation(counts = log10(colSums(medullo)))
cor_matrix <- cor(medullo, method = "pearson")
cor_distance <- as.dist(1 - cor_matrix)
hc <- hclust(cor_distance, method = "ward.D")
Heatmap(cor_matrix, 
        top_annotation = column_ha,
        name = "correlation",
        cluster_columns = hc,
        cluster_rows = hc,
        show_column_names = F,
        show_row_names = F,
        show_row_dend = F,
        col = colorRamp2(c(-0.2, 0, 0.2), c("blue", "white", "red")))



dataset = "Greenleaf_brain"
annotation = "Greenleaf_brain_same+as+paper+but+Early+RG+Late+RG-RG_Unk-rm"
brain_metadata_gl = get_filtered_cells(cell_num_filter = 1, 
                                       annotation = annotation, 
                                       dataset = dataset)
brain_cells = c("brain RG GL_Br",
                "brain mGPC GL_Br",
                "brain OPC/Oligo GL_Br",
                "brain MG GL_Br")

scatac_df_brain = scatac_df_gl[, brain_cells]
corrs_pearson = cor(scatac_df_brain, brain_mutations)

subtype_colors <- c(
  "GBM" = "red",
  "Astrocytoma" = "green",
  "Oligo" = "blue"
)

ha <- HeatmapAnnotation(subtype = donor_to_subtype, col = list(subtype = subtype_colors))
Heatmap(scale(corrs_pearson), 
        col = RColorBrewer::brewer.pal(9, "RdBu"),
        column_names_gp = grid::gpar(fontsize = 2),
        row_names_gp = grid::gpar(fontsize = 8),
        top_annotation = ha)

corrs_pearson = cor(scatac_df_brain, astrocytoma)
Heatmap(scale(corrs_pearson), 
        col = RColorBrewer::brewer.pal(9, "RdBu"),
        column_names_gp = grid::gpar(fontsize = 2),
        row_names_gp = grid::gpar(fontsize = 8))

corrs_pearson = cor(scatac_df_brain, gbm)
Heatmap(scale(corrs_pearson), 
        col = RColorBrewer::brewer.pal(9, "RdBu"),
        column_names_gp = grid::gpar(fontsize = 2),
        row_names_gp = grid::gpar(fontsize = 8))

corrs_pearson = cor(scatac_df_brain, oligo)
Heatmap(scale(corrs_pearson), 
        col = RColorBrewer::brewer.pal(9, "RdBu"),
        column_names_gp = grid::gpar(fontsize = 2),
        row_names_gp = grid::gpar(fontsize = 8))

medullo = read.csv("../data/mutation_data/CNS-Medullo.txt",
               sep="\t")
medullo = medullo[, 4:ncol(medullo)]
medullo = medullo[chr_ranges %in% chr_keep, ]

donor_to_subtype = read.csv("../data/processed_data/medullo_donors_subtyped.csv")
donor_to_subtype = data.frame(donor_to_subtype, row.names = "X")
donor_to_subtype = donor_to_subtype[colnames(medullo), ]

cells_for_medullo = c("cerebellum Astrocytes SH",
                      "brain GluN GL_Br",
                      "cerebellum Granule neurons SH",
                      "brain Limbic system neurons SH",
                      "cerebellum Purkinje neurons SH",
                      "brain IN GL_Br",
                      "brain SKOR2_NPSR1 positive cells SH")
scatac_df_sh = scatac_df_sh[chr_ranges %in% chr_keep, ]
scatac_df_brain = cbind(scatac_df_sh, scatac_df_gl)
scatac_df_brain = scatac_df_brain[, cells_for_medullo]
corrs_pearson = cor(scatac_df_brain, medullo)

subtype_colors <- c(
  "Large Cell" = "red",
  "Desmoplastic" = "green",
  "Medullo NOS" = "blue"
)
ha <- HeatmapAnnotation(subtype = donor_to_subtype, 
                        col = list(subtype = subtype_colors))

Heatmap(scale(corrs_pearson), 
        col = RColorBrewer::brewer.pal(9, "RdBu"),
        column_names_gp = grid::gpar(fontsize = 2),
        row_names_gp = grid::gpar(fontsize = 8),
        top_annotation = ha)

##### Kidney #####
# scatac_df_yang = t(readRDS("../data/processed_data/count_overlap_data/combined_count_overlaps/Yang_kidney_remove_cell_number_distinctions/Yang_kidney_combined_count_overlaps.rds"))
# scatac_df_yang = scatac_df_yang[chr_keep, ]
# cancer_samples_atac = read.csv("../data/processed_data/binned_atac/binned_atac_KIRP-US.csv",
#                                row.names = 1)
# cancer_samples_atac = cancer_samples_atac[chr_keep, ]

kidney = read.csv("../data/processed_data/mutations_with_subtypes/kidney_all.csv")

subtype = kidney[["subtype"]]

subtype_colors <- c(
  "Adenocarcinoma, clear cell type" = "red",
  "Adenocarcinoma, papillary type" = "green",
  "Adenocarcinoma, chromophobe type" = "blue"
)

kidney = kidney[, 3:ncol(kidney)]
kidney = t(kidney)
column_ha <- HeatmapAnnotation(subtype = subtype, 
                               counts = log10(colSums(kidney)),
                               col = list(subtype = 
                                            subtype_colors))
cor_matrix <- cor(kidney, method = "pearson")
cor_distance <- as.dist(1 - cor_matrix)
hc <- hclust(cor_distance, method = "average")
Heatmap(cor_matrix, 
        top_annotation = column_ha,
        name = "correlation",
        cluster_columns = hc,
        cluster_rows = hc,
        show_column_names = F,
        show_row_names = F,
        show_row_dend = F,
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")))

if (!file.exists("../data/processed_data/cancer_atac_50k_var_features.rds")) {
  cancer_samples_atac = readRDS("../data/normalized_atac_pan_peak_set.rds")
  features=rownames(cancer_samples_atac)
  cancer_samples_atac = cancer_samples_atac[, grepl("KIRP", colnames(cancer_samples_atac))]
  sample_names = unique(substr(colnames(cancer_samples_atac), 1, 41))
  cancer_samples_atac <- sapply(sample_names, function(i) {
    rowMeans(cancer_samples_atac[,
            grepl(i, colnames(cancer_samples_atac))])
  })
  cancer_samples_atac = as.data.frame(cancer_samples_atac)
  seurat_data = CreateSeuratObject(cancer_samples_atac)
  seurat_data = FindVariableFeatures(seurat_data, nfeatures = 50000,
                                     selection.method="vst")
  var_features <- VariableFeatures(seurat_data)
  var_features = gsub("-", "_", var_features)
  cancer_samples_atac = cancer_samples_atac[var_features, ]
  saveRDS(cancer_samples_atac, "../data/processed_data/cancer_atac_50k_var_features.rds")
} else {
  cancer_samples_atac = readRDS("../data/processed_data/cancer_atac_50k_var_features.rds")
}

if (!file.exists("../data/processed_data/cancer_atac_50k_var_features.rds")) {
  scatac_df_yang = t(readRDS("../data/processed_data/count_overlap_data/combined_count_overlaps/default_annotation/interval_ranges_yang_Yang_kidney_combined_count_overlaps.rds"))
  cell_types = colnames(scatac_df_yang)
  scatac_df_yang = cpm(scatac_df_yang, log=T, prior.count=5)
  scatac_df_yang = normalize.quantiles(scatac_df_yang)
  colnames(scatac_df_yang) = cell_types
  rownames(scatac_df_yang) = features
  seurat_data = CreateSeuratObject(counts = scatac_df_yang)
  seurat_data = FindVariableFeatures(seurat_data, nfeatures = 50000,
                                     selection.method="vst")
  var_features <- VariableFeatures(seurat_data)
  var_features = gsub("-", "_", var_features)
  scatac_df_yang = scatac_df_yang[var_features, ]
  saveRDS(scatac_df_yang, "../data/processed_data/scatac_50k_var_features.rds")
} else {
  scatac_df_yang = readRDS("../data/processed_data/scatac_50k_var_features.rds")
}

features_keep = intersect(rownames(scatac_df_yang), rownames(cancer_samples_atac))
scatac_df_yang = scatac_df_yang[features_keep, ]
cancer_samples_atac = cancer_samples_atac[features_keep, ]
corrs_pearson = cor(scatac_df_yang, cancer_samples_atac)

pdf("kidney_atac_corrs.pdf", width=15, height=10)
Heatmap(corrs_pearson, 
        col = RColorBrewer::brewer.pal(9, "RdBu"),
        column_names_gp = grid::gpar(fontsize = 2),
        row_names_gp = grid::gpar(fontsize = 8),
        cell_fun = create_cell_fun(corrs = corrs_pearson, fs=3))
dev.off()


#### Skin ####
skin = read.csv("../data/processed_data/mutations_with_subtypes/all_melanoma.csv",
                row.names = 1)

subtype = skin[["subtype"]]
subtype[subtype == ""] = "NA"
subtype_colors <- c(
  "NA" = "red",
  "Superficial Spreading" = "green",
  "Acral Lentiginuous" = "blue",
  "Desmoplastic" = "purple",
  "Nodular" = "black",
  "Mucosal Lentiginous" = "orange"
)

skin = skin[, grepl("chr[0-9]", colnames(skin))]
skin = t(skin)
column_ha <- HeatmapAnnotation(subtype = subtype, 
                               counts = log10(colSums(skin)),
                               col = list(subtype = 
                                            subtype_colors))
cor_matrix <- cor(skin, method = "pearson")
cor_distance <- as.dist(1 - cor_matrix)
hc <- hclust(cor_distance, method = "ward.D")
Heatmap(cor_matrix, 
        top_annotation = column_ha,
        name = "correlation",
        cluster_columns = hc,
        cluster_rows = hc,
        show_column_names = F,
        show_row_names = F,
        show_row_dend = F,
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")))

