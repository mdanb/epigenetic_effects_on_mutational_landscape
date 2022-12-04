#!/bin/bash

#Rscript plot_top_features_iteratively.R --all_cells --bing_ren --cell_number_filter=100 --cancer_types=Skin-Melanoma,Liver-HCC,ColoRect-AdenoCA,CNS-GBM,Eso-AdenoCA,Lung-AdenoCA,Lung-SCC 
Rscript plot_top_features.R --all_cells --combined_datasets --cell_number_filter=100 --cancer_types=Skin-Melanoma,Breast-AdenoCA,Liver-HCC,ColoRect-AdenoCA,CNS-GBM,Eso-AdenoCA,Lung-AdenoCA,Lung-SCC

