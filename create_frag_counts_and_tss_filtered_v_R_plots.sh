# Melanoma Skin cells only
#if ! [ -f "figures/melanoma_skin_cells_only_num_frags_vs_correlation.png" ]; then      #	Rscript create_frag_counts_and_tss_filtered_v_R_plots.R \
#	Rscript create_frag_counts_and_tss_filtered_v_R_plots.R \
#	--cancer_type "Skin.Melanoma" \
#	--boxplot_cell_types "Skin Sun Exposed Melanocyte (BR)-Skin Melanocyte (BR)-Skin Sun Exposed Fibroblast (Epithelial) (BR)-Skin Fibroblast (Epithelial) (BR)-Skin Keratinocyte 1 (BR)-Skin Sun Exposed Keratinocyte 1 (BR)-Skin T Lymphocyte 1 (CD8+) (BR)-Skin Sun Exposed T Lymphocyte 1 (CD8+) (BR)-Skin T lymphocyte 2 (CD4+) (BR)-Skin Sun Exposed T lymphocyte 2 (CD4+) (BR)-Skin Macrophage (General,Alveolar) (BR)-Skin Sun Exposed Macrophage (General,Alveolar) (BR)" \
#	--tissue_for_tsse_filtered_cell_types "Bing Ren-Skin,Bing Ren-Skin Sun Exposed" \
#	--tsse_filtered_cell_types "Fibroblast (Epithelial)-Melanocyte,Fibroblast (Epithelial)-Melanocyte" \
#	--plot_filename "melanoma_skin_cells_only_num_frags_vs_correlation.png" \
#	--plot_x_ticks "1000,10000,50000,100000,150000,250000,300000,400000,500000,1000000,1500000,2000000,5000000,10000000,15000000,20000000"
#fi

# Melanoma
if ! [ -f "figures/melanoma_num_frags_vs_correlation.png" ]; then
	Rscript create_frag_counts_and_tss_filtered_v_R_plots.R \
	--cancer_type "Skin.Melanoma" \
	--boxplot_cell_types "Skin Sun Exposed Melanocyte (BR),Heart Schwann cells (SH),Muscle Type II Skeletal Myocyte (BR),Stomach Stromal cells (SH),Lung Ciliated epithelial cells (SH)" \
	--tissue_for_tsse_filtered_cell_types "Bing Ren-Skin Sun Exposed,Shendure-Heart,Bing Ren-Muscle,Shendure-Stomach,Shendure-Lung" \
	--tsse_filtered_cell_types "Melanocyte,Schwann cells,Type II Skeletal Myocyte,Stromal cells,Ciliated epithelial cells" \
	--plot_filename "melanoma_num_frags_vs_correlation.png" \
	--plot_x_ticks "1000,10000,50000,100000,150000,250000,300000,400000,500000,1000000,1500000,2000000,5000000,10000000,15000000,20000000"
fi

if ! [ -f "figures/melanoma_BR_only_num_frags_vs_correlation.png" ]; then               
    	Rscript create_frag_counts_and_tss_filtered_v_R_plots.R \
    	--cancer_type "Skin.Melanoma" \
    	--boxplot_cell_types "Skin Sun Exposed Melanocyte (BR),Muscle Type II Skeletal Myocyte (BR),Artery Aorta Smooth Muscle (General) (BR),Small Intestine Smooth Muscle (General) (BR),Adrenal Gland Cortical Epithelial-like (BR),Esophagus Ge Junction Macrophage (General) (BR),Thyroid Fibroblast (General) (BR),Artery Aorta Fibroblast (General) (BR),Adipose Omentum Fibroblast (General) (BR),Artery Aorta Vascular Smooth Muscle 1 (BR)" \
    	--tissue_for_tsse_filtered_cell_types "Bing Ren-Skin Sun Exposed,Bing Ren-Muscle,Bing Ren-Artery Aorta,Bing Ren-Small Intestine,Bing Ren-Adrenal Gland,Bing Ren-Esophagus Ge Junction,Bing Ren-Thyroid,Bing Ren-Artery Aorta,Bing Ren-Adipose Omentum,Bing Ren-Artery Aorta" \
    	--tsse_filtered_cell_types "Melanocyte,Type II Skeletal Myocyte,Smooth Muscle (General),Smooth Muscle (General),Cortical Epithelial-like,Macrophage (General),Fibroblast (General),Fibroblast (General),Fibroblast (General),Vascular Smooth Muscle 1" \
    	--plot_filename "melanoma_BR_only_num_frags_vs_correlation.png" \
    	--plot_x_ticks "1000,10000,50000,100000,150000,250000,300000,400000,500000,1000000,1500000,2000000,5000000,10000000,15000000,20000000"
fi                                                                              

# Eso.AdenoCA
if ! [ -f "figures/esophagus_adenoca_num_frags_vs_correlation.png" ]; then
	Rscript create_frag_counts_and_tss_filtered_v_R_plots.R \
	--cancer_type "Eso.AdenoCA" \
	--boxplot_cell_types "Stomach Goblet cells (SH),Stomach Foveolar Cell (BR),Colon Transverse Colon Epithelial Cell 2 (BR),Esophagus Muscularis Foveolar Cell (BR),Small Intestine Small Intestinal Enterocyte (BR)" \
	--tsse_filtered_cell_types "Goblet cells,Foveolar Cell,Colon Epithelial Cell 2,Foveolar Cell,Small Intestinal Enterocyte" \
	--tissue_for_tsse_filtered_cell_types "Shendure-Stomach,Bing Ren-Stomach,Bing Ren-Colon Transverse,Bing Ren-Esophagus Muscularis,Bing Ren-Small Intestine" \
	--plot_filename "esophagus_adenoca_num_frags_vs_correlation.png" \
	--plot_x_ticks "1000,10000,50000,100000,150000,250000,300000,400000,500000,1000000,1500000,2000000,5000000,10000000,15000000,20000000"
fi

# Eso.AdenoCA                                                                   
if ! [ -f "figures/esophagus_adenoca_BR_only_num_frags_vs_correlation.png" ]; then     
   Rscript create_frag_counts_and_tss_filtered_v_R_plots.R \
   --cancer_type "Eso.AdenoCA" \
   --boxplot_cell_types "Lung Club Cell (BR),Stomach Foveolar Cell (BR),Colon Transverse Colon Epithelial Cell 2 (BR),Esophagus Muscularis Foveolar Cell (BR),Small Intestine Small Intestinal Enterocyte (BR)" \
   --tsse_filtered_cell_types "Club Cell,Foveolar Cell,Colon Epithelial Cell 2,Foveolar Cell,Small Intestinal Enterocyte" \
   --tissue_for_tsse_filtered_cell_types "Bing Ren-Lung,Bing Ren-Stomach,Bing Ren-Colon Transverse,Bing Ren-Esophagus Muscularis,Bing Ren-Small Intestine" \
   --plot_filename "esophagus_adenoca_BR_only_num_frags_vs_correlation.png" \
   --plot_x_ticks "1000,10000,50000,100000,150000,250000,300000,400000,500000,1000000,1500000,2000000,5000000,10000000,15000000,20000000"
fi                                                                             


# Lung AdenoCA AT2 only
if ! [ -f "figures/lung_adenoca_at2_only_frags_vs_correlation.png" ]; then
	Rscript create_frag_counts_and_tss_filtered_v_R_plots.R \
	--cancer_type "Lung.AdenoCA" \
	--boxplot_cell_types "Distal Lung AT2 (TS),Lung Alveolar Type 2 (AT2) Cell (BR),Lung Bronchiolar and alveolar epithelial cells (SH)" \
	--tissue_for_tsse_filtered_cell_types "Tsankov-Lung,Bing Ren-Lung,Shendure-Lung" \
	--tsse_filtered_cell_types "AT2,Alveolar Type 2 (AT2) Cell,Bronchiolar and alveolar epithelial cells" \
	--plot_filename "lung_adenoca_at2_only_frags_vs_correlation.png" \
	--plot_x_ticks "1000,10000,50000,100000,150000,250000,300000,400000,500000,1000000,1500000,2000000,5000000,10000000,15000000,20000000"
fi

# Lung AdenoCA
if ! [ -f "figures/lung_adenoca_frags_vs_correlation.png" ]; then
	Rscript create_frag_counts_and_tss_filtered_v_R_plots.R \
	--cancer_type "Lung.AdenoCA" \
	--boxplot_cell_types "Distal Lung AT2 (TS)-Lung Ciliated epithelial cells (SH)-Lung Club Cell (BR)-Muscle Type II Skeletal Myocyte (BR)-Lung Bronchiolar and alveolar epithelial cells (SH)" \
	--tissue_for_tsse_filtered_cell_types "Tsankov-Lung,Shendure-Lung,Bing Ren-Lung,Bing Ren-Muscle,Shendure-Lung" \
	--tsse_filtered_cell_types "AT2,Ciliated epithelial cells,Club Cell,Type II Skeletal Myocyte,Bronchiolar and alveolar epithelial cells" \
	--plot_filename "lung_adenoca_frags_vs_correlation.png" \
	--plot_x_ticks "1000,10000,50000,100000,150000,250000,300000,400000,500000,1000000,1500000,2000000,5000000,10000000,15000000,20000000"
fi

if ! [ -f "figures/lung_adenoca_BR_only_num_frags_vs_correlation.png" ]; then               
    Rscript create_frag_counts_and_tss_filtered_v_R_plots.R \
    --cancer_type "Lung.AdenoCA" \
    --boxplot_cell_types "Lung Club Cell (BR),Muscle Type II Skeletal Myocyte (BR),Heart Lv Endothelial Cell (Myocardial) (BR),Small Intestine Smooth Muscle (General) (BR),Artery Tibial Vascular Smooth Muscle 1 (BR)" \
    --tissue_for_tsse_filtered_cell_types "Bing Ren-Lung,Bing Ren-Muscle,Bing Ren-Heart Lv,Bing Ren-Small Intestine,Bing Ren-Artery Tibial" \
    --tsse_filtered_cell_types "Club Cell,Type II Skeletal Myocyte,Endothelial Cell (Myocardial),Smooth Muscle (General),Vascular Smooth Muscle 1" \
    --plot_filename "lung_adenoca_BR_only_frags_vs_correlation.png" \
    --plot_x_ticks "1000,10000,50000,100000,150000,250000,300000,400000,500000,1000000,1500000,2000000,5000000,10000000,15000000,20000000"
fi

# Breast AdenoCA
if ! [ -f "figures/breast_adenoca_num_frags_vs_correlation.png" ]; then
	Rscript create_frag_counts_and_tss_filtered_v_R_plots.R \
	--cancer_type "Breast.AdenoCA" \
	--boxplot_cell_types "Mammary Tissue Basal Epithelial (Mammary) (BR),Stomach Goblet cells (SH),Mammary Tissue Mammary Luminal Epithelial Cell 1 (BR),Mammary Tissue Mammary Luminal Epithelial Cell 2 (BR),Esophagus Mucosa Airway Goblet Cell (BR)" \
	--tissue_for_tsse_filtered_cell_types "Bing Ren-Mammary Tissue,Shendure-Stomach,Bing Ren-Mammary Tissue,Bing Ren-Mammary Tissue,Bing Ren-Esophagus Mucosa" \
	--tsse_filtered_cell_types "Basal Epithelial (Mammary),Goblet cells,Mammary Luminal Epithelial Cell 1,Mammary Luminal Epithelial Cell 2,Airway Goblet Cell" \
	--plot_filename "breast_adenoca_num_frags_vs_correlation.png" \
	--plot_x_ticks "1000,10000,50000,100000,150000,250000,300000,400000,500000,1000000,1500000,2000000,5000000,10000000,15000000,20000000"
fi

if ! [ -f "figures/breast_adenoca_BR_only_num_frags_vs_correlation.png" ]; then         
    Rscript create_frag_counts_and_tss_filtered_v_R_plots.R \
    --cancer_type "Breast.AdenoCA" \
    --boxplot_cell_types "Mammary Tissue Basal Epithelial (Mammary) (BR),Mammary Tissue Mammary Luminal Epithelial Cell 1 (BR),Mammary Tissue Mammary Luminal Epithelial Cell 2 (BR),Esophagus Mucosa Airway Goblet Cell (BR),Heart Atrial Appendage Cardiac Pericyte 1 (BR)" \
    --tissue_for_tsse_filtered_cell_types "Bing Ren-Mammary Tissue,Bing Ren-Mammary Tissue,Bing Ren-Mammary Tissue,Bing Ren-Esophagus Mucosa,Bing Ren-Heart Atrial Appendage" \
    --tsse_filtered_cell_types "Basal Epithelial (Mammary),Mammary Luminal Epithelial Cell 1,Mammary Luminal Epithelial Cell 2,Airway Goblet Cell,Cardiac Pericyte 1" \
    --plot_filename "breast_adenoca_BR_only_num_frags_vs_correlation.png" \
    --plot_x_ticks "1000,10000,50000,100000,150000,250000,300000,400000,500000,1000000,1500000,2000000,5000000,10000000,15000000,20000000"
fi

# CNS GBM
if ! [ -f "figures/cns_gbm_frags_vs_correlation.png" ]; then
	Rscript create_frag_counts_and_tss_filtered_v_R_plots.R \
	--cancer_type "CNS.GBM" \
	--boxplot_cell_types "Brain Astrocytes/Oligodendrocytes (SH)-Brain Astrocytes (SH)-Nerve Tibial Schwann Cell (General) (BR)-Heart Lv Cardiac Pericyte 1 (BR)-Muscle Type I Skeletal Myocyte (BR)" \
	--tsse_filtered_cell_types "Astrocytes/Oligodendrocytes,Astrocytes,Schwann Cell (General),Cardiac Pericyte 1,Type I Skeletal Myocyte" \
	--tissue_for_tsse_filtered_cell_types "Shendure-Brain,Shendure-Brain,Bing Ren-Nerve Tibial,Bing Ren-Heart Lv,Bing Ren-Muscle" \
	--plot_filename "cns_gbm_frags_vs_correlation.png" \
	--plot_x_ticks "1000,10000,50000,100000,150000,250000,300000,400000,500000,1000000,1500000,2000000,5000000,10000000,15000000"
fi

if ! [ -f "figures/cns_gbm_BR_only_num_frags_vs_correlation.png" ]; then                    
    Rscript create_frag_counts_and_tss_filtered_v_R_plots.R \
    --cancer_type "CNS.GBM" \
    --boxplot_cell_types "Nerve Tibial Schwann Cell (General) (BR),Heart Lv Cardiac Pericyte 1 (BR),Muscle Type I Skeletal Myocyte (BR),Lung Club Cell (BR),Heart Lv Cardiac Pericyte 2 (BR)" \
    --tsse_filtered_cell_types "Schwann Cell (General),Cardiac Pericyte 1,Type I Skeletal Myocyte,Club Cell,Cardiac Pericyte 2" \
    --tissue_for_tsse_filtered_cell_types "Bing Ren-Nerve Tibial,Bing Ren-Heart Lv,Bing Ren-Muscle,Bing Ren-Lung,Bing Ren-Heart Lv" \
    --plot_filename "cns_gbm_BR_only_num_frags_vs_correlation.png" \
    --plot_x_ticks "1000,10000,50000,100000,150000,250000,300000,400000,500000,1000000,1500000,2000000,5000000,10000000,15000000,20000000"
fi

# ColoRect AdenoCA
if ! [ -f "figures/colorect_adenoca_num_frags_vs_correlation.png" ]; then
	Rscript create_frag_counts_and_tss_filtered_v_R_plots.R \
	--cancer_type "ColoRect.AdenoCA" \
	--boxplot_cell_types "Intestine Intestinal epithelial cells (SH),Stomach Goblet cells (SH),Placenta PAEP_MECOM positive cells (SH),Colon Transverse Colon Epithelial Cell 2 (BR),Liver Erythroblasts (SH)" \
	--tsse_filtered_cell_types "Intestinal epithelial cells,Goblet cells,PAEP_MECOM positive cells,Colon Epithelial Cell 2,Erythroblasts" \
	--tissue_for_tsse_filtered_cell_types "Shendure-Intestine,Shendure-Stomach,Shendure-Placenta,Bing Ren-Colon Transverse,Shendure-Liver" \
	--plot_filename "colorect_adenoca_num_frags_vs_correlation.png" \
	--plot_x_ticks "1000,10000,50000,100000,150000,250000,300000,400000,500000,1000000,1500000,2000000,5000000"
fi

if ! [ -f "figures/colorect_adenoca_BR_only_num_frags_vs_correlation.png" ]; then       
    Rscript create_frag_counts_and_tss_filtered_v_R_plots.R \
    --cancer_type "ColoRect.AdenoCA" \
    --boxplot_cell_types "Colon Transverse Colon Epithelial Cell 2 (BR),Colon Transverse Colonic Goblet Cell (BR),Colon Transverse Colon Epithelial Cell 1 (BR)" \
    --tsse_filtered_cell_types "Colon Epithelial Cell 2;Colonic Goblet Cell;Colon Epithelial Cell 1" \
    --tissue_for_tsse_filtered_cell_types "Bing Ren-Colon Transverse" \
    --plot_filename "colorect_adenoca_BR_only_num_frags_vs_correlation.png" \
    --plot_x_ticks "1000,10000,50000,100000,150000,250000,300000,400000,500000,1000000,1500000,2000000,5000000,10000000,15000000,20000000"
fi

# Lung SCC
if ! [ -f "figures/lung_scc_num_frags_vs_correlation.png" ]; then
	Rscript create_frag_counts_and_tss_filtered_v_R_plots.R \
	--cancer_type "Lung.SCC" \
	--boxplot_cell_types "Proximal Lung Basal (TS),Esophagus Mucosa Esophageal Epithelial Cell (BR),Artery Aorta Smooth Muscle (General) (BR),Proximal Lung Secretory (TS),Stomach Goblet cells (SH)" \
	--tsse_filtered_cell_types "Basal,Esophageal Epithelial Cell,Smooth Muscle (General),Secretory,Goblet cells" \
	--tissue_for_tsse_filtered_cell_types "Tsankov-Lung,Bing Ren-Esophagus Mucosa,Bing Ren-Artery Aorta,Tsankov-Lung,Shendure-Stomach" \
	--plot_filename "lung_scc_num_frags_vs_correlation.png" \
	--plot_x_ticks "1000,10000,50000,100000,150000,250000,300000,400000,500000,1000000,1500000,2000000,5000000,10000000,15000000,20000000"
fi

if ! [ -f "figures/lung_scc_BR_only_num_frags_vs_correlation.png" ]; then               
    Rscript create_frag_counts_and_tss_filtered_v_R_plots.R \
    --cancer_type "Lung.SCC" \
    --boxplot_cell_types "Esophagus Mucosa Esophageal Epithelial Cell (BR),Artery Aorta Smooth Muscle (General) (BR),Lung Club Cell (BR),Small Intestine T Lymphocyte 1 (CD8+) (BR),Vagina Smooth Muscle (General) (BR)" \
    --tsse_filtered_cell_types "Esophageal Epithelial Cell,Smooth Muscle (General),Club Cell,T Lymphocyte 1 (CD8+),Smooth Muscle (General)" \
    --tissue_for_tsse_filtered_cell_types "Bing Ren-Esophagus Mucosa,Bing Ren-Artery Aorta,Bing Ren-Lung,Bing Ren-Small Intestine,Bing Ren-Vagina" \
    --plot_filename "lung_scc_BR_only_num_frags_vs_correlation.png" \
    --plot_x_ticks "1000,10000,50000,100000,150000,250000,300000,400000,500000,1000000,1500000,2000000,5000000,10000000,15000000,20000000"
fi

# Liver HCC
if ! [ -f "figures/liver_hcc_num_frags_vs_correlation.png" ]; then
	Rscript create_frag_counts_and_tss_filtered_v_R_plots.R \
	--cancer_type "Liver.HCC" \
	--boxplot_cell_types "Liver Hepatoblasts (SH),Liver Hepatocyte (BR),Artery Aorta Smooth Muscle (General) (BR),Eye Horizontal cells/Amacrine cells? (SH),Nerve Tibial Endothelial Cell (General) 1 (BR)" \
	--tsse_filtered_cell_types "Hepatoblasts,Hepatocyte,Smooth Muscle (General),Horizontal cells/Amacrine cells?,Endothelial Cell (General) 1" \
	--tissue_for_tsse_filtered_cell_types "Shendure-Liver,Bing Ren-Liver,Bing Ren-Artery Aorta,Shendure-Eye,Bing Ren-Nerve Tibial" \
	--plot_filename "liver_hcc_num_frags_vs_correlation.png" \
	--plot_x_ticks "1000,10000,50000,100000,150000,250000,300000,400000,500000,1000000,1500000,2000000,5000000"
fi

if ! [ -f "figures/liver_hcc_BR_only_num_frags_vs_correlation.png" ]; then              
    Rscript create_frag_counts_and_tss_filtered_v_R_plots.R \
    --cancer_type "Liver.HCC" \
    --boxplot_cell_types "Liver Hepatocyte (BR),Artery Aorta Smooth Muscle (General) (BR),Nerve Tibial Endothelial Cell (General) 1 (BR),Small Intestine Small Intestinal Goblet Cell (BR),Muscle Type II Skeletal Myocyte (BR)" \
    --tsse_filtered_cell_types "Hepatocyte,Smooth Muscle (General),Endothelial Cell (General) 1,Small Intestinal Goblet Cell,Type II Skeletal Myocyte" \
    --tissue_for_tsse_filtered_cell_types "Bing Ren-Liver,Bing Ren-Artery Aorta,Bing Ren-Nerve Tibial,Bing Ren-Small Intestine,Bing Ren-Muscle" \
    --plot_filename "liver_hcc_BR_only_num_frags_vs_correlation.png" \
    --plot_x_ticks "1000,10000,50000,100000,150000,250000,300000,400000,500000,1000000,1500000,2000000,5000000,10000000,15000000,20000000"
fi 
