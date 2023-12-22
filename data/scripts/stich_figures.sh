#!/bin/bash

root="../../figures/models/XGB/Skin-Melanoma/scATAC_source_Bingren_Greenleaf_colon_Greenleaf_pbmc_bm_Shendure_Tsankov_Yang_kidney_cell_number_filter_100_annotation_finalized_annotation_seed_all_seeds/backwards_elimination_results/"

cancer_types=("Skin-Melanoma" "Liver-HCC" "ColoRect-AdenoCA" "multiple_myeloma" "Eso-AdenoCA" "CNS-GBM" "Lung-AdenoCA" "Lung-SCC")

for i in "${!cancer_types[@]}"; do
    cancer_types[$i]="${root}${cancer_types[$i]}_feature_importance_with_10_5_2_features_top_5_features.pdf
done

for i in "${cancer_types[@]}"; do
    echo "$i"
done
