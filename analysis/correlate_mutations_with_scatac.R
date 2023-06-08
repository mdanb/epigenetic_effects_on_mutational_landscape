library(ComplexHeatmap)
library(RColorBrewer)
library(tidyverse)
library(circlize)
library(edgeR)
library(preprocessCore)
library(Seurat)

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



scatac_df_br = t(readRDS("../data/processed_data/count_overlap_data/combined_count_overlaps/Bingren_remove_same_celltype_indexing/Bingren_combined_count_overlaps.rds"))
scatac_df_sh = t(readRDS("../data/processed_data/count_overlap_data/combined_count_overlaps/Shendure_remove_unknown_unsure/Shendure_combined_count_overlaps.rds"))
colnames(scatac_df_sh) = paste(colnames(scatac_df_sh), "SH")
colnames(scatac_df_br) = paste(colnames(scatac_df_br), "BR")
scatac_df = cbind(scatac_df_sh, scatac_df_br)
scatac_df = scatac_df[rownames(scatac_df) %in% chr_keep, ]

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

#### Meso ####
scatac_df_lung = t(readRDS("../data/processed_data/count_overlap_data/combined_count_overlaps/Tsankov_separate_fibroblasts/Tsankov_count_filter_1_combined_count_overlaps.rds"))
scatac_df_lung = scatac_df_lung[rownames(scatac_df_lung) %in% chr_keep, ]

mutations_df_meso = read.csv("../data/processed_data/mesothelioma.csv",
                             row.names = 1)
scatac_df_lung = scatac_df_lung[, !(colnames(scatac_df_lung) == "distal lung Immune")]

mutations_df_meso = mutations_df_meso[rownames(mutations_df_meso) %in% chr_keep, ]
corrs_pearson = cor(scatac_df_lung, mutations_df_meso)
Heatmap(scale(corrs_pearson), 
        col = RColorBrewer::brewer.pal(9, "RdBu"),
        column_names_gp = grid::gpar(fontsize = 8),
        row_names_gp = grid::gpar(fontsize = 8))

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

scatac_df_lung = t(readRDS("../data/processed_data/count_overlap_data/combined_count_overlaps/Tsankov_separate_fibroblasts/Tsankov_count_filter_1_combined_count_overlaps.rds"))
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

gbm = read.csv("../data/mutation_data/CNS-GBM.txt",
                            sep="\t")
gbm = gbm[, 4:ncol(gbm)]
gbm = gbm[chr_ranges %in% chr_keep, ]

astrocytoma = read.csv("../data/mutation_data/CNS-PiloAstro.txt",
               sep="\t")
astrocytoma = astrocytoma[, 4:ncol(astrocytoma)]
astrocytoma = astrocytoma[chr_ranges %in% chr_keep, ]

oligo = read.csv("../data/mutation_data/CNS-Oligo.txt",
                       sep="\t")
oligo = oligo[, 4:ncol(oligo)]
oligo = oligo[chr_ranges %in% chr_keep, ]

brain_mutations = cbind(gbm, astrocytoma, oligo)

donor_to_subtype = c(rep("GBM", ncol(gbm)),
                     rep("Astrocytoma", ncol(astrocytoma)),
                     rep("Oligo", ncol(oligo)))
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

if (!file.exists("../data/processed_data/cancer_atac_50k_var_features.rds")) {
  cancer_samples_atac = readRDS("../data/normalized_atac_pan_peak_set.rds")
  raw = readRDS("../data/raw_atac_pan_peak_set.rds")[grepl("KIRP",
                                                           colnames(raw))]
  features=rownames(cancer_samples_atac)
  cancer_samples_atac = cancer_samples_atac[, grepl("KIRP", colnames(cancer_samples_atac))]
  seurat_data = CreateSeuratObject(counts = raw)
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
  scatac_df_yang = t(readRDS("../data/processed_data/count_overlap_data/combined_count_overlaps/Yang_kidney_remove_cell_number_distinctions/interval_ranges_yang_Yang_kidney_combined_count_overlaps.rds"))
  rownames(scatac_df_yang) = features
  cell_types = colnames(scatac_df_yang)
  seurat_data = CreateSeuratObject(counts = scatac_df_yang)
  seurat_data = FindVariableFeatures(seurat_data, nfeatures = 50000,
                                     selection.method="vst")
  var_features <- VariableFeatures(seurat_data)
  var_features = gsub("-", "_", var_features)
  scatac_df_yang = cpm(scatac_df_yang, log=T, prior.count=5)
  scatac_df_yang = normalize.quantiles(scatac_df_yang)
  colnames(scatac_df_yang) = cell_types
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


