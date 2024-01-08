#!/bin/bash


cancer_types=("Skin-Melanoma" "Liver-HCC" "ColoRect-AdenoCA" "multiple_myeloma" "Eso-AdenoCa" "CNS-GBM" "Lung-AdenoCA" "Lung-SCC")
#cancer_types=("Skin-Melanoma" "ColoRect-AdenoCA" "CNS-GBM" "Lung-AdenoCA" "Lung-SCC")

side_by_side_pdfs=()

for i in "${!cancer_types[@]}"; do
    root="../../figures/models/XGB/${cancer_types[$i]}/scATAC_source_Bingren_Greenleaf_colon_Greenleaf_pbmc_bm_Shendure_Tsankov_Yang_kidney_cell_number_filter_100_annotation_finalized_annotation_seed_all_seeds/backwards_elimination_results/"
    figure1="${root}${cancer_types[$i]}_top_feature_appearances.pdf"
    figure2="${root}${cancer_types[$i]}_top_feature_test_set_perf_with_10_5_2_1_features.pdf"
    output_pdf="${root}${cancer_types[$i]}_side_by_side.pdf"
    pdfjam --fitpaper true --trim "0mm 40mm 0mm 40mm" --nup 2x1 "$figure1" "$figure2" --outfile "$output_pdf"

    side_by_side_pdfs+=("$output_pdf")
done

pdftk "${side_by_side_pdfs[@]}" cat output combined.pdf
pdfjam --nup 1x8 --landscape combined.pdf --outfile stitched.pdf --fitpaper true
