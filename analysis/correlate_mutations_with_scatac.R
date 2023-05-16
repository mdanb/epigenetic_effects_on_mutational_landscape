library(ComplexHeatmap)
library(RColorBrewer)
chr_keep = read.csv("../processed_data/chr_keep.csv")[["chr"]]
scatac_df = t(readRDS("../processed_data/count_overlap_data/combined_count_overlaps/Tsankov_separate_fibroblasts/Tsankov_count_filter_1_combined_count_overlaps.rds"))
mutations_df = read.csv("../processed_data/mesothelioma.csv",
                        row.names = 1)
scatac_df = scatac_df[rownames(scatac_df) %in% chr_keep, ]
scatac_df = scatac_df[, !(colnames(scatac_df) == "distal lung Immune")]
mutations_df = mutations_df[rownames(mutations_df) %in% chr_keep, ]
corrs_spearman = cor(scatac_df, mutations_df, method = "spearman")
corrs_pearson = cor(scatac_df, mutations_df)
Heatmap(scale(corrs_spearman), 
        col = RColorBrewer::brewer.pal(9, "RdBu"),
        column_names_gp = grid::gpar(fontsize = 8),
        row_names_gp = grid::gpar(fontsize = 8))
Heatmap(scale(corrs_pearson), 
        col = RColorBrewer::brewer.pal(9, "RdBu"),
        column_names_gp = grid::gpar(fontsize = 8),
        row_names_gp = grid::gpar(fontsize = 8))

meso_paper = corrs_pearson[, c("mesomics_top_r1_90th_perc", "mesomics_bottom_r1_10th_perc")]
Heatmap(scale(meso_paper), 
        col = RColorBrewer::brewer.pal(9, "RdBu"),
        column_names_gp = grid::gpar(fontsize = 8),
        row_names_gp = grid::gpar(fontsize = 8),
        cluster_columns = F,
        name = "Pearson's r, \nscaled per column")

df_top_fibroblast = data.frame(muts = mutations_df["mesomics_top_r1_90th_perc"],
                                scatac = scatac_df[, "distal lung Mesothelium"])
df_bottom_mesothelium = data.frame(muts = mutations_df["mesomics_bottom_r1_10th_perc"],
                                scatac = scatac_df[, "distal lung"])
df_top_fibroblast = data.frame(muts = mutations_df["mesomics_top_r1_90th_perc"],
                                scatac = scatac_df)
df_bottom_mesothelium = data.frame(muts = mutations_df["mesomics_bottom_r1_10th_perc"],
                                   scatac = scatac_df)

