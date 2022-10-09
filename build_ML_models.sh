#!/bin/bash

#python3 build_ML_model.py --cancer_types Skin-Melanoma Liver-HCC ColoRect-AdenoCA CNS-GBM Eso-AdenoCA Lung-AdenoCA Lung-SCC --all_cells --bing_ren --scATAC_cell_number_filter=1

python3 build_ML_model.py --cancer_types Skin-Melanoma Liver-HCC ColoRect-AdenoCA CNS-GBM Eso-AdenoCA Lung-AdenoCA Lung-SCC --all_cells --bing_ren --scATAC_cell_number_filter=100

#python3 build_ML_model.py --cancer_types Ovary-AdenoCA CNS-PiloAstro CNS-Oligo Panc-Endocrine Kidney-RCC Prost-AdenoCA Thy-AdenoCA Lymph-BNHL Uterus-AdenoCA Breast-AdenoCA Panc-AdenoCA Head-SCC CNS-Medullo SoftTissue-Leiomyo Cervix-SCC Lymph-CLL SoftTissue-Liposarc Kidney-ChRCC Stomach-AdenoCA Bladder-TCC Myeloid-AML Biliary-AdenoCA Breast-LobularCA Cervix-AdenoCA Bone-Osteosarc Breast-DCIS Myeloid-MPN Myeloid-MDS Bone-Benign Bone-Epith --all_cells --combined_datasets --scATAC_cell_number_filter=100


#python3 build_ML_model.py --cancer_types Ovary-AdenoCA CNS-PiloAstro CNS-Oligo Panc-Endocrine Kidney-RCC Prost-AdenoCA Thy-AdenoCA Lymph-BNHL Uterus-AdenoCA Breast-AdenoCA Panc-AdenoCA Head-SCC CNS-Medullo SoftTissue-Leiomyo Cervix-SCC Lymph-CLL SoftTissue-Liposarc Kidney-ChRCC Stomach-AdenoCA Bladder-TCC Myeloid-AML Biliary-AdenoCA Breast-LobularCA Cervix-AdenoCA Bone-Osteosarc Breast-DCIS Myeloid-MPN Myeloid-MDS Bone-Benign Bone-Epith --all_cells --combined_datasets --scATAC_cell_number_filter=1

