#!/bin/bash

#python prep_dfs_for_feature_importance_plots.py --cancer_types Skin-Melanoma Liver-HCC ColoRect-AdenoCA CNS-GBM Eso-AdenoCA Lung-AdenoCA Lung-SCC --all_cells --bing_ren --cell_number_filter=1 --tss_filtered --tss_filtered_num_fragment_filter=100000
#python prep_dfs_for_feature_importance_plots.py --cancer_types Skin-Melanoma Liver-HCC ColoRect-AdenoCA CNS-GBM Eso-AdenoCA Lung-AdenoCA Lung-SCC --all_cells --combined_datasets --cell_number_filter=1

#python prep_dfs_for_feature_importance_plots.py --cancer_types sarcomatoid epithelioid --all_cells --tsankov --cell_number_filter=1
#python prep_dfs_for_feature_importance_plots.py --cancer_types Lung-AdenoCA Lung-SCC --all_cells --tsankov --cell_number_filter=1
#python prep_dfs_for_feature_importance_plots.py --cancer_types Skin-Melanoma Liver-HCC ColoRect-AdenoCA CNS-GBM Eso-AdenoCA Lung-AdenoCA Lung-SCC --all_cells --combined_datasets --cell_number_filter=1

# Lung AdenoCA, Lung SCC 100k filter, Broad only
# python3 prep_dfs_for_feature_importance_plots.py --cancer_types Lung-AdenoCA Lung-SCC --all_cells --tsankov --cell_number_filter=1 --tss_fragment_filter=100000

# Lung AdenoCA, Lung SCC, Broad only
# python3 prep_dfs_for_feature_importance_plots.py --cancer_types Lung-AdenoCA Lung-SCC --all_cells --tsankov --cell_number_filter=1

# 7 cancer types, all data, 1 cell filter, 100k filter
# python3 prep_dfs_for_feature_importance_plots.py --cancer_types Skin-Melanoma Liver-HCC ColoRect-AdenoCA CNS-GBM Eso-AdenoCA Lung-AdenoCA Lung-SCC --all_cells --bing_ren --tsankov --shendure --cell_number_filter=1 --tss_fragment_filter=100000

# 7 cancer types, all data, 1 cell filter, no fragment filter
# python3 prep_dfs_for_feature_importance_plots.py --cancer_types Skin-Melanoma Liver-HCC ColoRect-AdenoCA CNS-GBM Eso-AdenoCA Lung-AdenoCA Lung-SCC --all_cells --bing_ren --tsankov --shendure --cell_number_filter=1

# 7 cancer types, all data, 100 cell filter, no fragment filter
# python3 prep_dfs_for_feature_importance_plots.py --cancer_types Liver-HCC ColoRect-AdenoCA CNS-GBM Eso-AdenoCA Skin-Melanoma Lung-AdenoCA Lung-SCC --all_cells --bing_ren --tsankov --shendure --cell_number_filter=100

# 7 cancer types, all data, 100 cell filter, 100k fragment filter
# python3 prep_dfs_for_feature_importance_plots.py --cancer_types Skin-Melanoma Liver-HCC ColoRect-AdenoCA CNS-GBM Eso-AdenoCA Lung-AdenoCA Lung-SCC --all_cells --bing_ren --tsankov --shendure --cell_number_filter=100 --tss_fragment_filter=100000

# 7 cancer types, Bing Ren only, 1 cell filter, no fragment filter
# python3 prep_dfs_for_feature_importance_plots.py --cancer_types Skin-Melanoma Liver-HCC ColoRect-AdenoCA CNS-GBM Eso-AdenoCA Lung-AdenoCA Lung-SCC --all_cells --bing_ren --cell_number_filter=1

# 7 cancer types, Bing Ren only, 1 cell filter, 100k fragment filter
# python3 prep_dfs_for_feature_importance_plots.py --cancer_types Skin-Melanoma Liver-HCC ColoRect-AdenoCA CNS-GBM Eso-AdenoCA Lung-AdenoCA Lung-SCC --all_cells --bing_ren --cell_number_filter=1 --tss_fragment_filter=100000

# 7 cancer types, Bing Ren only, 1 cell filter, 200k fragment filter
# python3 prep_dfs_for_feature_importance_plots.py --cancer_types Skin-Melanoma Liver-HCC ColoRect-AdenoCA CNS-GBM Eso-AdenoCA Lung-AdenoCA Lung-SCC --all_cells --bing_ren --cell_number_filter=1 --tss_fragment_filter=200000

# 7 cancer types, Bing Ren only, 1 cell filter, 300k fragment filter
# python3 prep_dfs_for_feature_importance_plots.py --cancer_types Skin-Melanoma Liver-HCC ColoRect-AdenoCA CNS-GBM Eso-AdenoCA Lung-AdenoCA Lung-SCC --all_cells --bing_ren --cell_number_filter=1 --tss_fragment_filter=300000

# 7 cancer types, Bing Ren only, 1 cell filter, 400k fragment filter
# python3 prep_dfs_for_feature_importance_plots.py --cancer_types Skin-Melanoma Liver-HCC ColoRect-AdenoCA CNS-GBM Eso-AdenoCA Lung-AdenoCA Lung-SCC --all_cells --bing_ren --cell_number_filter=1 --tss_fragment_filter=400000

# 7 cancer types, Bing Ren only, 1 cell filter, 500k fragment filter
# python3 prep_dfs_for_feature_importance_plots.py --cancer_types Skin-Melanoma Liver-HCC ColoRect-AdenoCA CNS-GBM Eso-AdenoCA Lung-AdenoCA Lung-SCC --all_cells --bing_ren --cell_number_filter=1 --tss_fragment_filter=500000

# 7 cancer types, Bing Ren only, 1 cell filter, 1M fragment filter
# python3 prep_dfs_for_feature_importance_plots.py --cancer_types Skin-Melanoma Liver-HCC ColoRect-AdenoCA CNS-GBM Eso-AdenoCA Lung-AdenoCA Lung-SCC --all_cells --bing_ren --cell_number_filter=1 --tss_fragment_filter=1000000

# Biliary-AdenoCA, Bladder-TCC, Bone-Benign, Bone-Epith, Bone-Osteosarc, Breast-AdenoCA, Breast-DCIS, all data, 1 cell filter, no fragment filter
# python3 prep_dfs_for_feature_importance_plots.py --cancer_types Biliary-AdenoCA Bladder-TCC Bone-Benign Bone-Epith Bone-Osteosarc Breast-AdenoCA Breast-DCIS --all_cells --bing_ren --tsankov --shendure --cell_number_filter=1

# Breast-LobularCA, CNS-Medullo, CNS-Oligo, CNS-PiloAstro, Cervix-AdenoCA, Cervix-SCC, Head-SCC, all data, 1 cell filter, no fragment filter
# python3 prep_dfs_for_feature_importance_plots.py --cancer_types Breast-LobularCA CNS-Medullo CNS-Oligo CNS-PiloAstro Cervix-AdenoCA Cervix-SCC Head-SCC --all_cells --bing_ren --tsankov --shendure --cell_number_filter=1

# Kidney-ChRCC, Kidney-RCC, Lymph-BNHL, Lymph-CLL, Myeloid-AML, Myeloid-MDS, Myeloid-MPN, all data, 1 cell filter, no fragment filter
# python3 prep_dfs_for_feature_importance_plots.py --cancer_types Kidney-ChRCC Kidney-RCC Lymph-BNHL Lymph-CLL Myeloid-AML Myeloid-MDS Myeloid-MPN --all_cells --bing_ren --tsankov --shendure --cell_number_filter=1

# Ovary-AdenoCA, Panc-AdenoCA, Panc-Endocrine, Prost-AdenoCA, SoftTissue-Leiomyo, SoftTissue-Liposarc, Stomach-AdenoCA, all data, 1 cell filter, no fragment filter
# python3 prep_dfs_for_feature_importance_plots.py --cancer_types Ovary-AdenoCA Panc-AdenoCA Panc-Endocrine Prost-AdenoCA SoftTissue-Leiomyo SoftTissue-Liposarc Stomach-AdenoCA --all_cells --bing_ren --tsankov --shendure --cell_number_filter=1

# Thy-AdenoCA, Uterus-AdenoCA, all data, 1 cell filter, no fragment filter
# python3 prep_dfs_for_feature_importance_plots.py --cancer_types Thy-AdenoCA Uterus-AdenoCA --all_cells --bing_ren --tsankov --shendure --cell_number_filter=1


#############################################
#############################################
#############################################

################ Mesothelioma ###############

#python3 prep_dfs_for_feature_importance_plots.py --cancer_types sarcomatoid epithelioid --all_cells --datasets Tsankov --cell_number_filter=1 --waddell_sarc_biph
#python3 prep_dfs_for_feature_importance_plots.py --cancer_types sarcomatoid --all_cells --datasets Tsankov --cell_number_filter=1 --waddell_sarc
#python3 prep_dfs_for_feature_importance_plots.py --cancer_types sarcomatoid --all_cells --datasets Tsankov --cell_number_filter=1 --waddell_sarc_tsankov_sarc
#python3 prep_dfs_for_feature_importance_plots.py --cancer_types sarcomatoid --all_cells --datasets Tsankov --cell_number_filter=1 --waddell_sarc_biph_tsankov_sarc_biph


############################################
############################################
############################################

### BR, SH, TS, GL_Br, GL_BlBm, Y_k ###

python3 prep_dfs_for_feature_importance_plots.py --cancer_types CNS-GBM CNS-Medullo CNS-Oligo CNS-PiloAstro Head-SCC Kidney-ChRCC --all_cells --datasets Bingren Greenleaf_brain Greenleaf_pbmc_bm Shendure Tsankov Yang_kidney --cell_number_filter=1

python3 prep_dfs_for_feature_importance_plots.py --cancer_types Kidney-RCC Lymph-BNHL Lymph-CLL Myeloid-AML Myeloid-MDS Myeloid-MPN --all_cells --datasets Bingren Greenleaf_brain Greenleaf_pbmc_bm Shendure Tsankov Yang_kidney --cell_number_filter=1

python3 prep_dfs_for_feature_importance_plots.py --cancer_types Myeloid-AML Myeloid-MDS Myeloid-MPN --all_cells --datasets Bingren Greenleaf_brain Greenleaf_pbmc_bm Shendure Tsankov Yang_kidney --cell_number_filter=1

