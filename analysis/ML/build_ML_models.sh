#!/bin/bash

#python3 build_ML_model.py --cancer_types Skin-Melanoma Liver-HCC ColoRect-AdenoCA CNS-GBM Eso-AdenoCA Lung-AdenoCA Lung-SCC --all_cells --bing_ren --scATAC_cell_number_filter=1

#python3 build_ML_model.py --cancer_types Skin-Melanoma Liver-HCC ColoRect-AdenoCA CNS-GBM Eso-AdenoCA Lung-AdenoCA Lung-SCC --all_cells --bing_ren --scATAC_cell_number_filter=100


#python3 build_ML_model.py --cancer_types epithelioid sarcomatoid --all_cells --tsankov --scATAC_cell_number_filter=100 --meso

#ython3 build_ML_model.py --cancer_types sarcomatoid --all_cells --tsankov --scATAC_cell_number_filter=100 --meso_waddell_only

#python3 build_ML_model.py --cancer_types sarcomatoid_w786 --all_cells --tsankov --scATAC_cell_number_filter=100 --meso_waddell_and_broad_only

#python3 build_ML_model.py --cancer_types sarcomatoid --all_cells --tsankov --scATAC_cell_number_filter=1 --meso --tss_filtered --tss_filtered_num_fragment_filter=50000
#python3 build_ML_model.py --cancer_types epithelioid sarcomatoid --all_cells --tsankov --scATAC_cell_number_filter=1 --meso --tss_filtered --tss_filtered_num_fragment_filter=100000
#python3 build_ML_model.py --cancer_types sarcomatoid --all_cells --tsankov --scATAC_cell_number_filter=1 --meso --tss_filtered --tss_filtered_num_fragment_filter=300000

#python3 build_ML_model.py --cancer_types sarcomatoid --all_cells --tsankov --scATAC_cell_number_filter=1 --meso_waddell_only --tss_filtered --tss_filtered_num_fragment_filter=50000
#python3 build_ML_model.py --cancer_types epithelioid sarcomatoid --all_cells --tsankov --scATAC_cell_number_filter=1 --meso_waddell_only --tss_filtered --tss_filtered_num_fragment_filter=100000
#python3 build_ML_model.py --cancer_types sarcomatoid --all_cells --tsankov --scATAC_cell_number_filter=1 --meso_waddell_only --tss_filtered --tss_filtered_num_fragment_filter=300000

#python3 build_ML_model.py --cancer_types sarcomatoid_w786 --all_cells --tsankov --scATAC_cell_number_filter=1 --meso_waddell_and_broad_only --tss_filtered --tss_filtered_num_fragment_filter=50000
#python3 build_ML_model.py --cancer_types sarcomatoid_w786 --all_cells --tsankov --scATAC_cell_number_filter=1 --meso_waddell_and_broad_only --tss_filtered --tss_filtered_num_fragment_filter=100000
#python3 build_ML_model.py --cancer_types sarcomatoid_w786 --all_cells --tsankov --scATAC_cell_number_filter=1 --meso_waddell_and_broad_only --tss_filtered --tss_filtered_num_fragment_filter=300000

#python3 build_ML_model.py --cancer_types Lung-AdenoCA Lung-SCC --all_cells --tsankov --scATAC_cell_number_filter=1 

#python3 build_ML_model.py --cancer_types Ovary-AdenoCA CNS-PiloAstro CNS-Oligo Panc-Endocrine Kidney-RCC Prost-AdenoCA Thy-AdenoCA Lymph-BNHL Uterus-AdenoCA Breast-AdenoCA Panc-AdenoCA Head-SCC CNS-Medullo SoftTissue-Leiomyo Cervix-SCC Lymph-CLL SoftTissue-Liposarc Kidney-ChRCC Stomach-AdenoCA Bladder-TCC Myeloid-AML Biliary-AdenoCA Breast-LobularCA Cervix-AdenoCA Bone-Osteosarc Breast-DCIS Myeloid-MPN Myeloid-MDS Bone-Benign Bone-Epith --all_cells --combined_datasets --scATAC_cell_number_filter=100


#python3 build_ML_model.py --cancer_types Ovary-AdenoCA CNS-PiloAstro CNS-Oligo Panc-Endocrine Kidney-RCC Prost-AdenoCA Thy-AdenoCA Lymph-BNHL Uterus-AdenoCA Breast-AdenoCA Panc-AdenoCA Head-SCC CNS-Medullo SoftTissue-Leiomyo Cervix-SCC Lymph-CLL SoftTissue-Liposarc Kidney-ChRCC Stomach-AdenoCA Bladder-TCC Myeloid-AML Biliary-AdenoCA Breast-LobularCA Cervix-AdenoCA Bone-Osteosarc Breast-DCIS Myeloid-MPN Myeloid-MDS Bone-Benign Bone-Epith --all_cells --combined_datasets --scATAC_cell_number_filter=1

#python3 build_ML_model.py --cancer_types sarcomatoid --all_cells --tsankov --scATAC_cell_number_filter=1 --meso_waddell_biph_786_846 --tss_fragment_filter=300000
#python3 build_ML_model.py --cancer_types sarcomatoid --all_cells --tsankov --scATAC_cell_number_filter=1 --meso_waddell_biph_786_846 --tss_filtered --tss_filtered_num_fragment_filter=100000
#python3 build_ML_model.py --cancer_types sarcomatoid --all_cells --tsankov --scATAC_cell_number_filter=1 --meso_waddell_biph_786_846 --tss_filtered --tss_filtered_num_fragment_filter=50000


##############################################
##############################################
##############################################


# Lung AdenoCA, Lung SCC 100k filter, Broad only
#python3 build_ML_model.py --cancer_types Lung-AdenoCA Lung-SCC --all_cells --tsankov --scATAC_cell_number_filter=1 --tss_fragment_filter=100000

# Lung AdenoCA, Lung SCC, Broad only, no cell or fragment filter
# python3 build_ML_model.py --cancer_types Lung-AdenoCA Lung-SCC --all_cells --tsankov --scATAC_cell_number_filter=1

# 7 cancer types, all data, 100k filter, no cell filter
# python3 build_ML_model.py --cancer_types Skin-Melanoma Liver-HCC ColoRect-AdenoCA CNS-GBM Eso-AdenoCA Lung-AdenoCA Lung-SCC --all_cells --bing_ren --tsankov --shendure --scATAC_cell_number_filter=1 --tss_fragment_filter=100000

# 7 cancer types, all data, 1 cell filter, no fragment filter
#python3 build_ML_model.py --cancer_types Skin-Melanoma Liver-HCC ColoRect-AdenoCA CNS-GBM Eso-AdenoCA Lung-AdenoCA Lung-SCC --all_cells --bing_ren --tsankov --shendure --scATAC_cell_number_filter=1

# 7 cancer types, all data, 100 cell filter, no fragment filter
#python3 build_ML_model.py --cancer_types Skin-Melanoma Liver-HCC ColoRect-AdenoCA CNS-GBM Eso-AdenoCA Lung-AdenoCA Lung-SCC --all_cells --bing_ren --tsankov --shendure --scATAC_cell_number_filter=100

# 7 cancer types, all data, 100 cell filter, 100k fragment filter
# python3 build_ML_model.py --cancer_types Skin-Melanoma Liver-HCC ColoRect-AdenoCA CNS-GBM Eso-AdenoCA Lung-AdenoCA Lung-SCC --all_cells --bing_ren --tsankov --shendure --scATAC_cell_number_filter=100 --tss_fragment_filter=100000

# 7 cancer types, Bing Ren only, 1 cell filter, no fragment filter
# python3 build_ML_model.py --cancer_types Skin-Melanoma Liver-HCC ColoRect-AdenoCA CNS-GBM Eso-AdenoCA Lung-AdenoCA Lung-SCC --all_cells --scATAC_cell_number_filter=1 --bing_ren 

# 7 cancer types, Bing Ren only, 1 cell filter, 100k fragment filter
#python3 build_ML_model.py --cancer_types Skin-Melanoma Liver-HCC ColoRect-AdenoCA CNS-GBM Eso-AdenoCA Lung-AdenoCA Lung-SCC --all_cells --scATAC_cell_number_filter=1 --bing_ren --tss_fragment_filter=100000

# 7 cancer types, Bing Ren only, 1 cell filter, 200k fragment filter
# python3 build_ML_model.py --cancer_types Skin-Melanoma Liver-HCC ColoRect-AdenoCA CNS-GBM Eso-AdenoCA Lung-AdenoCA Lung-SCC --all_cells --scATAC_cell_number_filter=1 --bing_ren --tss_fragment_filter=200000

# 7 cancer types, Bing Ren only, 1 cell filter, 300k fragment filter
#python3 build_ML_model.py --cancer_types Skin-Melanoma Liver-HCC ColoRect-AdenoCA CNS-GBM Eso-AdenoCA Lung-AdenoCA Lung-SCC --all_cells --scATAC_cell_number_filter=1 --bing_ren --tss_fragment_filter=300000

# 7 cancer types, Bing Ren only, 1 cell filter, 400k fragment filter
# python3 build_ML_model.py --cancer_types Skin-Melanoma Liver-HCC ColoRect-AdenoCA CNS-GBM Eso-AdenoCA Lung-AdenoCA Lung-SCC --all_cells --scATAC_cell_number_filter=1 --bing_ren --tss_fragment_filter=400000

# 7 cancer types, Bing Ren only, 1 cell filter, 500k fragment filter
# python3 build_ML_model.py --cancer_types Skin-Melanoma Liver-HCC ColoRect-AdenoCA CNS-GBM Eso-AdenoCA Lung-AdenoCA Lung-SCC --all_cells --scATAC_cell_number_filter=1 --bing_ren --tss_fragment_filter=500000

# 7 cancer types, Bing Ren only, 1 cell filter, 1M fragment filter
# python3 build_ML_model.py --cancer_types Skin-Melanoma Liver-HCC ColoRect-AdenoCA CNS-GBM Eso-AdenoCA Lung-AdenoCA Lung-SCC --all_cells --scATAC_cell_number_filter=1 --bing_ren --tss_fragment_filter=1000000

# Biliary-AdenoCA, Bladder-TCC, Bone-Benign, Bone-Epith, Bone-Osteosarc, Breast-AdenoCA, Breast-DCIS, all data, 1 cell filter, no fragment filter
# python3 build_ML_model.py --cancer_types Biliary-AdenoCA Bladder-TCC Bone-Benign Bone-Epith Bone-Osteosarc Breast-AdenoCA Breast-DCIS --all_cells --bing_ren --tsankov --shendure --scATAC_cell_number_filter=1

# Breast-LobularCA, CNS-Medullo, CNS-Oligo, CNS-PiloAstro, Cervix-AdenoCA, Cervix-SCC, Head-SCC, all data, 1 cell filter, no fragment filter
# python3 build_ML_model.py --cancer_types Breast-LobularCA CNS-Medullo CNS-Oligo CNS-PiloAstro Cervix-AdenoCA Cervix-SCC Head-SCC --all_cells --bing_ren --tsankov --shendure --scATAC_cell_number_filter=1

# Kidney-ChRCC, Kidney-RCC, Lymph-BNHL, Lymph-CLL, Myeloid-AML, Myeloid-MDS, Myeloid-MPN, all data, 1 cell filter, no fragment filter
# python3 build_ML_model.py --cancer_types Kidney-ChRCC Kidney-RCC Lymph-BNHL Lymph-CLL Myeloid-AML Myeloid-MDS Myeloid-MPN --all_cells --bing_ren --tsankov --shendure --scATAC_cell_number_filter=1

# Ovary-AdenoCA, Panc-AdenoCA, Panc-Endocrine, Prost-AdenoCA, SoftTissue-Leiomyo, SoftTissue-Liposarc, Stomach-AdenoCA, all data, 1 cell filter, no fragment filter
# python3 build_ML_model.py --cancer_types Ovary-AdenoCA Panc-AdenoCA Panc-Endocrine Prost-AdenoCA SoftTissue-Leiomyo SoftTissue-Liposarc Stomach-AdenoCA --all_cells --bing_ren --tsankov --shendure --scATAC_cell_number_filter=1

# Thy-AdenoCA, Uterus-AdenoCA, all data, 1 cell filter, no fragment filter
# python3 build_ML_model.py --cancer_types Thy-AdenoCA Uterus-AdenoCA --all_cells --bing_ren --tsankov --shendure --scATAC_cell_number_filter=1

# 7 cancer types, with Greenleaf, 1 cell filter, no frag filter
# python3 build_ML_model.py --cancer_types Skin-Melanoma Liver-HCC ColoRect-AdenoCA CNS-GBM Eso-AdenoCA Lung-AdenoCA Lung-SCC --all_cells --datasets Bingren Greenleaf_brain Greenleaf_pbmc_bm Shendure Tsankov --scATAC_cell_number_filter=1

# Biliary-AdenoCA, Bladder-TCC, Bone-Benign, Bone-Epith, Bone-Osteosarc, Breast-AdenoCA, Breast-DCIS, all data with Greenleaf, 1 cell filter, no fragment filter
# python3 build_ML_model.py --cancer_types Biliary-AdenoCA Bladder-TCC Bone-Benign Bone-Epith Bone-Osteosarc Breast-AdenoCA Breast-DCIS --all_cells --datasets Bingren Greenleaf_brain Greenleaf_pbmc_bm Shendure Tsankov --scATAC_cell_number_filter=1

# Breast-LobularCA, CNS-Medullo, CNS-Oligo, CNS-PiloAstro, Cervix-AdenoCA, Cervix-SCC, Head-SCC, all data with Greenleaf, 1 cell filter, no fragment filter
# python3 build_ML_model.py --cancer_types Breast-LobularCA CNS-Medullo CNS-Oligo CNS-PiloAstro Cervix-AdenoCA Cervix-SCC Head-SCC --all_cells --datasets Bingren Greenleaf_brain Greenleaf_pbmc_bm Shendure Tsankov --scATAC_cell_number_filter=1

# Kidney-ChRCC, Kidney-RCC, Lymph-BNHL, Lymph-CLL, Myeloid-AML, Myeloid-MDS, Myeloid-MPN, all data with Greenleaf, 1 cell filter, no fragment filter
# python3 build_ML_model.py --cancer_types Kidney-ChRCC Kidney-RCC Lymph-BNHL Lymph-CLL Myeloid-AML Myeloid-MDS Myeloid-MPN --all_cells --datasets Bingren Greenleaf_brain Greenleaf_pbmc_bm Shendure Tsankov --scATAC_cell_number_filter=1

# Ovary-AdenoCA, Panc-AdenoCA, Panc-Endocrine, Prost-AdenoCA, SoftTissue-Leiomyo, SoftTissue-Liposarc, Stomach-AdenoCA, all data with Greenleaf, 1 cell filter, no fragment filter
python3 build_ML_model.py --cancer_types Ovary-AdenoCA Panc-AdenoCA Panc-Endocrine Prost-AdenoCA SoftTissue-Leiomyo SoftTissue-Liposarc Stomach-AdenoCA --all_cells --datasets Bingren Greenleaf_brain Greenleaf_pbmc_bm Shendure Tsankov --scATAC_cell_number_filter=1

# Thy-AdenoCA, Uterus-AdenoCA, all data, 1 cell filter, no fragment filter
# python3 build_ML_model.py --cancer_types Thy-AdenoCA Uterus-AdenoCA --all_cells --datasets Bingren Greenleaf_brain Greenleaf_pbmc_bm Shendure Tsankov --scATAC_cell_number_filter=1

