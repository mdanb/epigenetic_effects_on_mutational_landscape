#!/bin/bash

python prep_df_for_pie_charts.py --cancer_types Skin-Melanoma Liver-HCC ColoRect-AdenoCA CNS-GBM Eso-AdenoCA Lung-AdenoCA Lung-SCC --all_cells --bing_ren --cell_number_filter=100
python prep_df_for_pie_charts.py --cancer_types Skin-Melanoma Liver-HCC ColoRect-AdenoCA CNS-GBM Eso-AdenoCA Lung-AdenoCA Lung-SCC --all_cells --combined_datasets --cell_number_filter=100
