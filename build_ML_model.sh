#!/bin/bash

python3 build_ML_model.py --cancer_types Skin-Melanoma Liver-HCC ColoRect-AdenoCA CNS-GBM Eso-AdenoCA Lung-AdenoCA Lung-SCC --all_cells --combined_datasets
