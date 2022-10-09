#Rscript create_frag_counts_and_tss_filtered_v_R_plots.R \
#--cancer_type "Skin.Melanoma" \
#--boxplot_cell_types "Skin Sun Exposed Melanocyte (BR)-Skin Melanocyte (BR)-Skin Sun Exposed Fibroblast (Epithelial) (BR)-Skin Fibroblast (Epithelial) (BR)-Skin Keratinocyte 1 (BR)-Skin Sun Exposed Keratinocyte 1 (BR)-Skin T Lymphocyte 1 (CD8+) (BR)-Skin Sun Exposed T Lymphocyte 1 (CD8+) (BR)-Skin T lymphocyte 2 (CD4+) (BR)-Skin Sun Exposed T lymphocyte 2 (CD4+) (BR)-Skin Macrophage (General,Alveolar) (BR)-Skin Sun Exposed Macrophage (General,Alveolar) (BR)" \
#--tissue_for_tsse_filtered_cell_types "Bing Ren-Skin,Bing Ren-Skin Sun Exposed" \
#--tsse_filtered_cell_types "Fibroblast (Epithelial)-Melanocyte,Fibroblast (Epithelial)-Melanocyte" \
#--plot_filename "melanoma_num_frags_vs_correlation.png" \
#--plot_x_ticks "1000,10000,50000,100000,150000,250000,300000,400000,500000,600000"

#Rscript create_frag_counts_and_tss_filtered_v_R_plots.R \
#--cancer_type "Eso.AdenoCA" \
#--boxplot_cell_types "Stomach Goblet cells (SH)-Stomach Foveolar Cell (BR)-Colon Transverse Colon Epithelial Cell 2 (BR)-Esophagus Muscularis Foveolar Cell (BR)-Small Intestine Small Intestinal Enterocyte (BR)" \
#--tsse_filtered_cell_types "Goblet cells,Foveolar Cell,Colon Epithelial Cell 2,Foveolar Cell,Small Intestinal Enterocyte" \
#--tissue_for_tsse_filtered_cell_types "Shendure-Stomach,Bing Ren-Stomach,Bing Ren-Colon Transverse,Bing Ren-Esophagus Muscularis,Bing Ren-Small Intestine" \
#--plot_filename "esophagus_adenoca_num_frags_vs_correlation.png" \
#--plot_x_ticks "1000,10000,50000,100000,150000,250000,300000,400000,500000,600000"

Rscript create_frag_counts_and_tss_filtered_v_R_plots.R \
--cancer_type "Lung.AdenoCA" \
--boxplot_cell_types "Distal Lung AT2 (TS)-Lung Alveolar Type 2 (AT2) Cell (BR)-Lung Bronchiolar and alveolar epithelial cells (SH)" \
--tissue_for_tsse_filtered_cell_types "Tsankov-Lung,Bing Ren-Lung,Shendure-Lung" \
--tsse_filtered_cell_types "AT2,Alveolar Type 2 (AT2) Cell,Bronchiolar and alveolar epithelial cells" \
--plot_filename "lung_adenoca_at2_only_frags_vs_correlation.png" \
--plot_x_ticks "1000,10000,50000,100000,150000,250000,300000,400000,500000,600000,700000,800000,900000,1000000,2000000,3000000,4000000,5000000,6000000,8000000,10000000,15000000,20000000,30000000,40000000,50000000"

#Rscript create_frag_counts_and_tss_filtered_v_R_plots.R \
#--cancer_type "ColoRect.AdenoCA" \
#--boxplot_cell_types "Intestine Intestinal epithelial cells (SH)-Stomach Goblet cells (SH)-Placenta PAEP_MECOM positive cells (SH)-Colon Transverse Colon Epithelial Cell 2 (BR)-Colon Transverse Colon Epithelial Cell 2 (BR)-Liver Erythroblasts (SH)" \
#--tsse_filtered_cell_types "shendure-Intestine,shendure-Stomach,shendure-Placenta,bing_ren-Colon Transverse,bing_ren-Liver" \
#--plot_filename "colorect_adenoca_num_frags_vs_correlation.png" \


