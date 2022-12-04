#!/bin/bash

python prep_dfs_for_feature_importance_plots.R --cancer_types Skin-Melanoma Liver-HCC ColoRect-AdenoCA CNS-GBM Eso-AdenoCA Lung-AdenoCA Lung-SCC --all_cells --bing_ren --cell_number_filter=1
python prep_dfs_for_feature_importance_plots.R --cancer_types Skin-Melanoma Liver-HCC ColoRect-AdenoCA CNS-GBM Eso-AdenoCA Lung-AdenoCA Lung-SCC --all_cells --combined_datasets --cell_number_filter=1
