#!/bin/bash

#Rscript plot_top_features_iteratively.R --all_cells --bing_ren --cell_number_filter=100 --cancer_types=Skin-Melanoma,Liver-HCC,ColoRect-AdenoCA,CNS-GBM,Eso-AdenoCA,Lung-AdenoCA,Lung-SCC 
#Rscript plot_top_features.R --all_cells --combined_datasets --cell_number_filter=1 --cancer_types=Skin-Melanoma,Liver-HCC,ColoRect-AdenoCA,CNS-GBM,Eso-AdenoCA,Lung-AdenoCA,Lung-SCC --bar_plot_num_features=20

#Rscript plot_top_features.R --all_cells --bing_ren --cell_number_filter=1 --cancer_types=Skin-Melanoma,Liver-HCC,ColoRect-AdenoCA,CNS-GBM,Eso-AdenoCA,Lung-AdenoCA,Lung-SCC --bar_plot_num_features=20 --tss_filtered --tss_filtered_num_fragment_filter=100000

#Rscript plot_top_features.R --all_cells --tsankov --cell_number_filter=1 --cancer_types=epithelioid,sarcomatoid --bar_plot_num_features=20

# Lung AdenoCA, Lung SCC 100k filter, Broad only
#Rscript plot_top_features.R --all_cells --tsankov --cell_number_filter=1 --cancer_types=Lung-AdenoCA,Lung-SCC --bar_plot_num_features=20 --tss_fragment_filter=100000

# Lung AdenoCA, Lung SCC, Broad only
# Rscript plot_top_features.R --all_cells --tsankov --cell_number_filter=1 --cancer_types=Lung-AdenoCA,Lung-SCC --bar_plot_num_features=20 

# 7 cancer types, all data, 100k filter, 1 cell number filter
# Rscript plot_top_features.R --cancer_types=Skin-Melanoma,Liver-HCC,ColoRect-AdenoCA,CNS-GBM,Eso-AdenoCA,Lung-AdenoCA,Lung-SCC --all_cells --bing_ren --tsankov --shendure --cell_number_filter=1 --tss_fragment_filter=100000 --bar_plot_num_features=20

# 7 cancer types, all data, no fragment filter, 1 cell number filter
# Rscript plot_top_features.R --cancer_types=Skin-Melanoma,Liver-HCC,ColoRect-AdenoCA,CNS-GBM,Eso-AdenoCA,Lung-AdenoCA,Lung-SCC --all_cells --bing_ren --tsankov --shendure --cell_number_filter=1

# 7 cancer types, all data, no fragment filter, 100 cell number filter
# Rscript plot_top_features.R --cancer_types=Liver-HCC,ColoRect-AdenoCA,CNS-GBM,Eso-AdenoCA,Skin-Melanoma,Lung-AdenoCA,Lung-SCC --all_cells --bing_ren --tsankov --shendure --cell_number_filter=100 --bar_plot_num_features=20

# 7 cancer types, all data, 100k fragment filter, 100 cell number filter
# Rscript plot_top_features.R --cancer_types=Skin-Melanoma,Liver-HCC,ColoRect-AdenoCA,CNS-GBM,Eso-AdenoCA,Lung-AdenoCA,Lung-SCC --all_cells --bing_ren --tsankov --shendure --cell_number_filter=100 --tss_fragment_filter=100000 --bar_plot_num_features=20


# 7 cancer types, Bing Ren only, 1 cell filter, no fragment filter
# Rscript plot_top_features.R --cancer_types Skin-Melanoma,Liver-HCC,ColoRect-AdenoCA,CNS-GBM,Eso-AdenoCA,Lung-AdenoCA,Lung-SCC --all_cells --bing_ren --cell_number_filter=1

# 7 cancer types, Bing Ren only, 1 cell filter, 100k fragment filter
# Rscript plot_top_features.R --cancer_types Skin-Melanoma,Liver-HCC,ColoRect-AdenoCA,CNS-GBM,Eso-AdenoCA,Lung-AdenoCA,Lung-SCC --all_cells --bing_ren --cell_number_filter=1 --tss_fragment_filter=100000

# 7 cancer types, Bing Ren only, 1 cell filter, 200k fragment filter
# Rscript plot_top_features.R --cancer_types Skin-Melanoma,Liver-HCC,ColoRect-AdenoCA,CNS-GBM,Eso-AdenoCA,Lung-AdenoCA,Lung-SCC --all_cells --bing_ren --cell_number_filter=1 --tss_fragment_filter=200000

# 7 cancer types, Bing Ren only, 1 cell filter, 300k fragment filter
# Rscript plot_top_features.R --cancer_types Skin-Melanoma,Liver-HCC,ColoRect-AdenoCA,CNS-GBM,Eso-AdenoCA,Lung-AdenoCA,Lung-SCC --all_cells --bing_ren --cell_number_filter=1 --tss_fragment_filter=300000

# 7 cancer types, Bing Ren only, 1 cell filter, 400k fragment filter
# Rscript plot_top_features.R --cancer_types Skin-Melanoma,Liver-HCC,ColoRect-AdenoCA,CNS-GBM,Eso-AdenoCA,Lung-AdenoCA,Lung-SCC --all_cells --bing_ren --cell_number_filter=1 --tss_fragment_filter=400000

# 7 cancer types, Bing Ren only, 1 cell filter, 500k fragment filter
Rscript plot_top_features.R --cancer_types Skin-Melanoma,Liver-HCC,ColoRect-AdenoCA,CNS-GBM,Eso-AdenoCA,Lung-AdenoCA,Lung-SCC --all_cells --bing_ren --cell_number_filter=1 --tss_fragment_filter=500000

# 7 cancer types, Bing Ren only, 1 cell filter, 1m fragment filter
Rscript plot_top_features.R --cancer_types Skin-Melanoma,Liver-HCC,ColoRect-AdenoCA,CNS-GBM,Eso-AdenoCA,Lung-AdenoCA,Lung-SCC --all_cells --bing_ren --cell_number_filter=1 --tss_fragment_filter=1000000


# Biliary-AdenoCA, Bladder-TCC, Bone-Benign, Bone-Epith, Bone-Osteosarc, Breast-AdenoCA, Breast-DCIS, all data, 1 cell filter, no fragment filter
# Rscript plot_top_features.R --cancer_types Biliary-AdenoCA,Bladder-TCC,Bone-Benign,Bone-Epith,Bone-Osteosarc,Breast-AdenoCA,Breast-DCIS --all_cells --bing_ren --tsankov --shendure --cell_number_filter=1

# Breast-LobularCA, CNS-Medullo, CNS-Oligo, CNS-PiloAstro, Cervix-AdenoCA, Cervix-SCC, Head-SCC, all data, 1 cell filter, no fragment filter
# Rscript plot_top_features.R --cancer_types Breast-LobularCA,CNS-Medullo,CNS-Oligo,CNS-PiloAstro,Cervix-AdenoCA,Cervix-SCC,Head-SCC --all_cells --bing_ren --tsankov --shendure --cell_number_filter=1

# Kidney-ChRCC, Kidney-RCC, Lymph-BNHL, Lymph-CLL, Myeloid-AML, Myeloid-MDS, Myeloid-MPN, all data, 1 cell filter, no fragment filter
# Rscript plot_top_features.R --cancer_types Kidney-ChRCC,Kidney-RCC,Lymph-BNHL,Lymph-CLL,Myeloid-AML,Myeloid-MDS,Myeloid-MPN --all_cells --bing_ren --tsankov --shendure --cell_number_filter=1

# Ovary-AdenoCA, Panc-AdenoCA, Panc-Endocrine, Prost-AdenoCA, SoftTissue-Leiomyo, SoftTissue-Liposarc, Stomach-AdenoCA, all data, 1 cell filter, no fragment filter
# Rscript plot_top_features.R --cancer_types Ovary-AdenoCA,Panc-AdenoCA,Panc-Endocrine,Prost-AdenoCA,SoftTissue-Leiomyo,SoftTissue-Liposarc,Stomach-AdenoCA --all_cells --bing_ren --tsankov --shendure --cell_number_filter=1

# Thy-AdenoCA, Uterus-AdenoCA, all data, 1 cell filter, no fragment filter
# Rscript plot_top_features.R --cancer_types Thy-AdenoCA,Uterus-AdenoCA --all_cells --bing_ren --tsankov --shendure --cell_number_filter=1


