library(ComplexHeatmap)
library(RColorBrewer)

chr_keep = read.csv("../data/processed_data/chr_keep.csv")[["chr"]]
chr_ranges = unlist(read.csv("../data/processed_data/chr_ranges.csv"))

#### Meso ####
scatac_df = t(readRDS("../data/processed_data/count_overlap_data/combined_count_overlaps/Tsankov_separate_fibroblasts/Tsankov_count_filter_1_combined_count_overlaps.rds"))
scatac_df = scatac_df[rownames(scatac_df) %in% chr_keep, ]

mutations_df_meso = read.csv("../data/processed_data/mesothelioma.csv",
                             row.names = 1)
scatac_df = scatac_df[, !(colnames(scatac_df) == "distal lung Immune")]

mutations_df_meso = mutations_df_meso[rownames(mutations_df_meso) %in% chr_keep, ]
corrs_spearman = cor(scatac_df, mutations_df_meso, method = "spearman")
corrs_pearson = cor(scatac_df, mutations_df_meso)
Heatmap(scale(corrs_spearman), 
        col = RColorBrewer::brewer.pal(9, "RdBu"),
        column_names_gp = grid::gpar(fontsize = 8),
        row_names_gp = grid::gpar(fontsize = 8))
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
                                scatac = scatac_df[, "distal lung Mesothelium"])
df_bottom_mesothelium = data.frame(muts = mutations_df_meso["mesomics_bottom_r1_10th_perc"],
                                scatac = scatac_df[, "distal lung"])
df_top_fibroblast = data.frame(muts = mutations_df_meso["mesomics_top_r1_90th_perc"],
                                scatac = scatac_df)
df_bottom_mesothelium = data.frame(muts = mutations_df_meso["mesomics_bottom_r1_10th_perc"],
                                   scatac = scatac_df)

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
scatac_df = t(readRDS("../data/processed_data/count_overlap_data/combined_count_overlaps/Tsankov_separate_fibroblasts/Tsankov_count_filter_1_combined_count_overlaps.rds"))
scatac_df = scatac_df[rownames(scatac_df) %in% chr_keep, ]

subtype_colors <- c(
  "Biphasic" = "red",
  "Epithelioid" = "green",
  "Sarcomatoid" = "blue"
)

ha <- HeatmapAnnotation(subtype = subtypes, col = list(subtype = subtype_colors))

corrs_pearson = cor(scatac_df, meso_individuals)
Heatmap(scale(corrs_pearson), 
        col = RColorBrewer::brewer.pal(9, "RdBu"),
        column_names_gp = grid::gpar(fontsize = 4),
        row_names_gp = grid::gpar(fontsize = 8),
        top_annotation = ha)

scatac_df = scatac_df[, !(colnames(scatac_df) == "distal lung Immune")]
corrs_pearson = cor(scatac_df, meso_individuals)

Heatmap(scale(corrs_pearson), 
        col = RColorBrewer::brewer.pal(9, "RdBu"),
        column_names_gp = grid::gpar(fontsize = 4),
        row_names_gp = grid::gpar(fontsize = 8),
        top_annotation = ha)

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

#### Lung-AdenoCA ####
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

scatac_df_lung = scatac_df_lung[, !(colnames(scatac_df_lung) == "distal lung Immune")]
corrs_pearson = cor(scatac_df_lung, lung_adenoca)
ha <- HeatmapAnnotation(subtype = adeno_subtypes, col = list(subtype = subtype_colors))
Heatmap(scale(corrs_pearson), 
        col = RColorBrewer::brewer.pal(9, "RdBu"),
        column_names_gp = grid::gpar(fontsize = 2),
        row_names_gp = grid::gpar(fontsize = 8),
        top_annotation = ha)

subtype_colors <- c(
  "AT.NKX2-1" = "red",
  "mTOR.SOX2" = "green",
  "COL.SOX4" = "blue",
  "KRT13,4.FOXA1" = "black",
  "CellCycle" = "purple",
  "S100A.HES2" = "cyan"
)

corrs_pearson = cor(scatac_df_lung, lung_squamous)
ha <- HeatmapAnnotation(subtype = lung_subtypes, col = list(subtype = subtype_colors))
Heatmap(scale(corrs_pearson), 
        col = RColorBrewer::brewer.pal(9, "RdBu"),
        column_names_gp = grid::gpar(fontsize = 2),
        row_names_gp = grid::gpar(fontsize = 8),
        top_annotation = ha)

####################

scatac_df_br = t(readRDS("../data/processed_data/count_overlap_data/combined_count_overlaps/Bingren_remove_same_celltype_indexing/Bingren_combined_count_overlaps.rds"))
scatac_df_sh = t(readRDS("../data/processed_data/count_overlap_data/combined_count_overlaps/Shendure_remove_unknown_unsure/Shendure_combined_count_overlaps.rds"))
colnames(scatac_df_sh) = paste(colnames(scatac_df_sh), "SH")
colnames(scatac_df_br) = paste(colnames(scatac_df_br), "BR")
scatac_df = cbind(scatac_df_sh, scatac_df_br)
scatac_df = scatac_df[rownames(scatac_df) %in% chr_keep, ]

#### Pancreas-AdenoCA ####
pancreas_adenoca = read.csv("../data/mutation_data/Panc-AdenoCA.txt",
                            sep="\t")
pancreas_adenoca = pancreas_adenoca[chr_ranges %in% chr_keep, ]
pancreas_adenoca = pancreas_adenoca[, 4:ncol(pancreas_adenoca)]

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
combined_annotation = cbind(subtype = donor_to_subtype, dataset = dataset_annotation,
                            rna_subtype = pancreas_rna_subtypes)
combined_colors = list(subtype = subtype_colors, dataset = dataset_colors,
                       rna_subtype = rna_subtype_colors)
# ha <- HeatmapAnnotation(subtype = donor_to_subtype, col = list(subtype = subtype_colors))
ha <- HeatmapAnnotation(df = data.frame(combined_annotation), col = combined_colors)

column_order <- order(scale(corrs_pearson)["stomach Goblet cells SH", ])
Heatmap(scale(corrs_pearson), 
        col = RColorBrewer::brewer.pal(9, "RdBu"),
        column_names_gp = grid::gpar(fontsize = 2),
        row_names_gp = grid::gpar(fontsize = 8),
        top_annotation = ha,
        column_order = column_order)
column_order <- order(scale(corrs_pearson)["stomach Foveolar Cell BR", ])
Heatmap(scale(corrs_pearson), 
        col = RColorBrewer::brewer.pal(9, "RdBu"),
        column_names_gp = grid::gpar(fontsize = 2),
        row_names_gp = grid::gpar(fontsize = 8),
        top_annotation = ha,
        column_order = column_order)


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
column_order <- order(scale(corrs_pearson["stomach Foveolar Cell BR", ]))
Heatmap(scale(corrs_pearson), 
        col = RColorBrewer::brewer.pal(9, "RdBu"),
        column_names_gp = grid::gpar(fontsize = 2),
        row_names_gp = grid::gpar(fontsize = 8),
        column_order = column_order)

column_order <- order(corrs_pearson["stomach Goblet cells SH", ])
Heatmap(scale(corrs_pearson), 
        col = RColorBrewer::brewer.pal(9, "RdBu"),
        column_names_gp = grid::gpar(fontsize = 2),
        row_names_gp = grid::gpar(fontsize = 8),
        column_order = column_order)

column_order <- order(corrs_pearson["intestine Intestinal epithelial cells SH", ])
Heatmap(scale(corrs_pearson), 
        col = RColorBrewer::brewer.pal(9, "RdBu"),
        column_names_gp = grid::gpar(fontsize = 2),
        row_names_gp = grid::gpar(fontsize = 8),
        column_order = column_order)

column_order <- order(corrs_pearson["colon_transverse Colonic Goblet Cell BR", ])
Heatmap(scale(corrs_pearson), 
        col = RColorBrewer::brewer.pal(9, "RdBu"),
        column_names_gp = grid::gpar(fontsize = 2),
        row_names_gp = grid::gpar(fontsize = 8),
        column_order = column_order)
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
