#!/bin/bash


#cancer_types=("Skin-Melanoma" "Liver-HCC" "ColoRect-AdenoCA" "multiple_myeloma" "Eso-AdenoCa" "CNS-GBM" "Lung-AdenoCA" "Lung-SCC" "waddell_combined" "SCLC")
#
#pdfs=()
#
#for i in "${!cancer_types[@]}"; do
#    root="../../figures/models/XGB/${cancer_types[$i]}/scATAC_source_Bingren_Greenleaf_colon_Greenleaf_pbmc_bm_Shendure_Tsankov_Yang_kidney_cell_number_filter_100_annotation_finalized_annotation_seed_all_seeds/backwards_elimination_results/"
#    figure="${root}${cancer_types[$i]}_feature_importance_with_10_5_2_features_top_5_features.pdf"
#    pdfs+=("$figure")
#done
#
#pdftk "${pdfs[@]}" cat output combined_feat_imp.pdf
#pdfjam --nup 1x10 --landscape combined_feat_imp.pdf --outfile stitched_feat_imp.pdf --fitpaper true
#
#
#
#cancer_types=("kidney_papillary" "kidney_clear_cell" "Kidney-ChRCC" "Panc-AdenoCA" "Panc-Endocrine")
#
#pdfs=()
#
#for i in "${!cancer_types[@]}"; do
#    root="../../figures/models/XGB/${cancer_types[$i]}/scATAC_source_Bingren_Greenleaf_colon_Greenleaf_pbmc_bm_Shendure_Tsankov_Yang_kidney_cell_number_filter_100_annotation_finalized_annotation_seed_all_seeds/backwards_elimination_results/"
#    figure="${root}${cancer_types[$i]}_feature_importance_with_10_5_2_features_top_5_features.pdf"
#    pdfs+=("$figure")
#done
#
#root="../../figures/models/XGB/msi_high/scATAC_source_Bingren_Greenleaf_colon_Greenleaf_pbmc_bm_Shendure_Tsankov_Yang_kidney_cell_number_filter_100_annotation_Greenleaf_colon_remove_cancer_merge_normal_unaffected_seed_all_seeds/backwards_elimination_results/"
#figure="${root}msi_high_feature_importance_with_10_5_2_features_top_5_features.pdf"
#pdfs+=("$figure")
#
#root="../../figures/models/XGB/Panc-Endocrine/scATAC_source_Bingren_Shendure_cell_number_filter_100_tissues_to_consider_pancreas_stomach_annotation_finalized_annotation_seed_all_seeds/backwards_elimination_results/"
#figure="${root}Panc-Endocrine_feature_importance_with_10_5_2_features_top_5_features.pdf"
#pdfs+=("$figure")
#
#pdftk "${pdfs[@]}" cat output combined_feat_imp_fig3.pdf
#pdfjam --nup 1x7 --landscape combined_feat_imp_fig3.pdf --outfile stitched_feat_imp_fig3.pdf --fitpaper true
#

#pdfs=("Lymph-CLL" "Myeloid-MPN")
#
#for i in "${!cancer_types[@]}"; do
#    root="../../figures/models/XGB/${cancer_types[$i]}/scATAC_source_Bingren_Greenleaf_colon_Greenleaf_pbmc_bm_Shendure_Tsankov_Yang_kidney_cell_number_filter_100_annotation_intermediate_blood_bm_annotation_seed_all_seeds/backwards_elimination_results/"
#    figure="${root}${cancer_types[$i]}_feature_importance_with_10_5_2_features_top_5_features.pdf"
#    pdfs+=("$figure")
#done
#
#root="../../figures/models/XGB/mss/scATAC_source_Bingren_Greenleaf_colon_Greenleaf_pbmc_bm_Shendure_Tsankov_Yang_kidney_cell_number_filter_100_annotation_finalized_annotation_seed_all_seeds/backwards_elimination_results/"
#figure="${root}mss_high_feature_importance_with_10_5_2_features_top_5_features.pdf"
#pdfs+=("$figure")
#
#pdftk "${pdfs[@]}" cat output combined_feat_imp_fig2.pdf
#pdfjam --nup 1x7 --landscape combined_feat_imp_fig2.pdf --outfile stitched_feat_imp_fig2.pdf --fitpaper true
#





#cancer_types=("CNS-GBM" "CNS-Medullo" "CNS-PiloAstro" "CNS-Oligo")
#
#for i in "${!cancer_types[@]}"; do
#	root="../../figures/models/XGB/${cancer_types[$i]}/scATAC_source_Bingren_Greenleaf_brain_Greenleaf_colon_Greenleaf_pbmc_bm_Shendure_Tsankov_Yang_kidney_cell_number_filter_100_annotation_finalized_annotation_seed_all_seeds/backwards_elimination_results/"
#    	#echo $i
#    	figure="${root}${cancer_types[$i]}_feature_importance_with_10_5_2_features_top_5_features.pdf"
#    	pdfs+=("$figure")
#done


#pdftk "${pdfs[@]}" cat output combined_feat_imp_fig4.pdf
#pdfjam --nup 1x7 --landscape combined_feat_imp_fig4.pdf --outfile stitched_feat_imp_fig4.pdf --fitpaper true



#cancer_types=("Lymph-CLL" "Myeloid-MPN")
#
#for i in "${!cancer_types[@]}"; do
#        root="../../figures/models/XGB/${cancer_types[$i]}/scATAC_source_Bingren_Greenleaf_colon_Greenleaf_pbmc_bm_Shendure_Tsankov_Yang_kidney_cell_number_filter_100_annotation_intermediate_blood_bm_annotation_seed_all_seeds/backwards_elimination_results/"
#        #echo $i
#        figure="${root}${cancer_types[$i]}_feature_importance_with_10_5_2_features_top_5_features.pdf"
#        pdfs+=("$figure")
#done
#
#root="../../figures/models/XGB/mss/scATAC_source_Bingren_Greenleaf_colon_Greenleaf_pbmc_bm_Shendure_Tsankov_Yang_kidney_cell_number_filter_100_annotation_finalized_annotation_seed_all_seeds/backwards_elimination_results/"
##echo $i
#figure="${root}mss_feature_importance_with_10_5_2_features_top_5_features.pdf"
#pdfs+=("$figure")
#
#pdftk "${pdfs[@]}" cat output combined_feat_imp_blood.pdf
#pdfjam --nup 1x7 --landscape combined_feat_imp_blood.pdf --outfile stitched_feat_imp_blood.pdf --fitpaper true
#


#cancer_types=("Bone-Osteosarc" "Breast-LobularCa" "Cervix-SCC" "Head-SCC" "Ovary-AdenoCA" "Prost-AdenoCA" "Uterus-AdenoCA")
#
#for i in "${!cancer_types[@]}"; do
#        root="../../figures/models/XGB/${cancer_types[$i]}/scATAC_source_Bingren_Greenleaf_colon_Greenleaf_pbmc_bm_Shendure_Tsankov_Yang_kidney_cell_number_filter_100_annotation_finalized_annotation_seed_all_seeds/backwards_elimination_results/"
#        figure="${root}${cancer_types[$i]}_feature_importance_with_10_5_2_features_top_5_features.pdf"
#        pdfs+=("$figure")
#done
#
#pdftk "${pdfs[@]}" cat output combined_feat_imp_remaining.pdf
#pdfjam --nup 1x7 --landscape combined_feat_imp_remaining.pdf --outfile stitched_feat_imp_remaining.pdf --fitpaper true
#

cancer_types=("Breast-AdenoCa" "Bone-Leiomyo" "Thy-AdenoCA" "Biliary-AdenoCA" "Bladder-TCC" "Stomach-AdenoCA")
for i in "${!cancer_types[@]}"; do
        root="../../figures/models/XGB/${cancer_types[$i]}/scATAC_source_Bingren_Greenleaf_colon_Greenleaf_pbmc_bm_Shendure_Tsankov_Yang_kidney_cell_number_filter_100_annotation_finalized_annotation_seed_all_seeds/backwards_elimination_results/"
        figure="${root}${cancer_types[$i]}_feature_importance_with_10_5_2_features_top_5_features.pdf"
        pdfs+=("$figure")
done

cancer_types=("Lymph-BNHL" "Myeloid-AML")
for i in "${!cancer_types[@]}"; do
        root="../../figures/models/XGB/${cancer_types[$i]}/scATAC_source_Bingren_Greenleaf_colon_Greenleaf_pbmc_bm_Shendure_Tsankov_Yang_kidney_cell_number_filter_100_annotation_intermediate_blood_bm_annotation_seed_all_seeds/backwards_elimination_results/"
        figure="${root}${cancer_types[$i]}_feature_importance_with_10_5_2_features_top_5_features.pdf"
        pdfs+=("$figure")
done


pdftk "${pdfs[@]}" cat output combined_feat_imp_fig4.pdf
pdfjam --nup 1x8 --landscape combined_feat_imp_fig4.pdf --outfile stitched_feat_imp_fig4.pdf --fitpaper true

