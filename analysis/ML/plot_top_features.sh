#Rscript plot_top_features.R  --datasets="Tsankov" --cancer_types=Lung-SCC --cell_number_filter=30 --ML_model=XGB --annotation=finalized_annotation --top_features_to_plot=2,5,10,15,19 --seed=3 --feature_importance_method=permutation_importance

#Rscript plot_top_features.R  --datasets="Tsankov" --cancer_types=Lung-SCC --top_features_to_plot=1,2,5,10,15 --cell_number_filter=1 --ML_model=XGB --annotation=finalized_annotation --feature_importance_method=permutation_importance --robustness_analysis --robustness_seed_range=1-5
#Rscript plot_top_features.R --robustness_analysis --datasets="Tsankov" --cancer_types=Lung-AdenoCA --cell_number_filter=30 --ML_model=XGB --annotation=finalized_annotation --top_features_to_plot=2,5,10,15,19 --robustness_test_perf_boxplot

#Rscript plot_top_features.R --datasets="Rawlins_fetal_lung,Wang_lung" --cancer_types=SCLC --cell_number_filter=100 --ML_model=XGB --annotation=triangle_cells_combine_meso_and_neuroendocrine_no_schwann --top_features_to_plot="1,2,5,10,15" --feature_importance_method=permutation_importance  --seed_range=1-5

#Rscript plot_top_features.R --datasets="Rawlins_fetal_lung,Tsankov" --cancer_types=SCLC --cell_number_filter=1 --ML_model=XGB --annotation=remove_basal --top_features_to_plot="1,2,5,10,15" --feature_importance_method=permutation_importance  --seed_range=1-5

#Rscript plot_top_features.R --datasets="Rawlins_fetal_lung,Tsankov" --cancer_types="sclcxtop_500-scale_TRUE-norm_by_mut_counts_TRUE-pca_n_components_30-dims_for_dim_reduction_1to15-regress_by_counts_FALSE-min_mutations_per_sample0_clustering_res_1xcluster_0,sclcxtop_500-scale_TRUE-norm_by_mut_counts_TRUE-pca_n_components_30-dims_for_dim_reduction_1to15-regress_by_counts_FALSE-min_mutations_per_sample0_clustering_res_1xcluster_1,sclcxtop_500-scale_TRUE-norm_by_mut_counts_TRUE-pca_n_components_30-dims_for_dim_reduction_1to15-regress_by_counts_FALSE-min_mutations_per_sample0_clustering_res_1xcluster_2" --cell_number_filter=1 --ML_model=XGB --annotation=finalized_annotation --top_features_to_plot="1,2,5,10,15" --feature_importance_method=permutation_importance  --seed_range=1-5


#Rscript plot_top_features.R --datasets="Tsankov,Wang_lung" --cancer_types=Lung-AdenoCA,Lung-SCC --cell_number_filter=100 --ML_model=XGB --annotation=finalized_annotation --top_features_to_plot="1,2,5,10,15" --feature_importance_method=permutation_importance  --seed_range=1-5

#Rscript plot_top_features.R --datasets="Rawlins_fetal_lung,Wang_lung" --cancer_types=combined_meso --cell_number_filter=100 --ML_model=XGB --annotation="triangle_cells_combine_meso_and_neuroendocrine_no_schwann" --top_features_to_plot="1,2,5,10,15" --feature_importance_method=permutation_importance  --seed_range=1-5

#Rscript plot_top_features.R --datasets="Rawlins_fetal_lung,Tsankov" --cancer_types="sclcxtop_500-scale_TRUE-norm_by_mut_counts_TRUE-pca_n_components_30-dims_for_dim_reduction_1to15-regress_by_counts_FALSE-min_mutations_per_sample0_clustering_res_0.7xcluster_0,sclcxtop_500-scale_TRUE-norm_by_mut_counts_TRUE-pca_n_components_30-dims_for_dim_reduction_1to15-regress_by_counts_FALSE-min_mutations_per_sample0_clustering_res_0.7xcluster_1" --cell_number_filter=1 --ML_model=XGB --annotation=Tsankov_basal_refined --top_features_to_plot="1,2,5,10,15" --feature_importance_method=permutation_importance  --seed_range=1-5


#Rscript plot_top_features.R --datasets="Rawlins_fetal_lung,Tsankov" --cancer_types=SCLC.A,SCLC.N,SCLC.P,SCLC.Y --cell_number_filter=1 --ML_model=XGB --annotation=Tsankov_basal_refined --top_features_to_plot="1,2,5,10,15" --feature_importance_method=permutation_importance --seed_range=1-5 --fold_for_test_set=1

#Rscript plot_top_features.R --datasets="Rawlins_fetal_lung,Tsankov" --cancer_types=SCLC.A,SCLC.N,SCLC.P,SCLC.Y --cell_number_filter=1 --ML_model=XGB --annotation=triangle_cells_combine_meso_and_neuroendocrine_no_schwann --top_features_to_plot="1,2,5,10,15" --feature_importance_method=permutation_importance --seed_range=1-5 --fold_for_test_set=1

#Rscript plot_top_features.R --datasets="Bingren,Shendure,Tsankov,Rawlins_fetal_lung,Yang_kidney,Greenleaf_brain,Greenleaf_pbmc_bm,Greenleaf_colon" --cancer_types=ColoRect-AdenoCA --cell_number_filter=100 --ML_model=XGB --annotation=remove_cancer_polyp_merge_normal_unaffected --top_features_to_plot="1,2,5,10" --feature_importance_method=permutation_importance --seed_range=1-10 --folds_for_test_set=1-10 --top_features_to_plot_feat_imp="1,2,5,10" --robustness_analysis --feat_imp_min_n_robustness=50
#Rscript plot_top_features.R --cancer_types=SCLC.Y --datasets=Tsankov,Rawlins_fetal_lung --cell_number_filter 1 --annotation finalized_annotation --seed_range 1-5 --top_features_to_plot "1,2,5,10,15" --feature_importance_method permutation_importance --folds_for_test_set 1-5 --robustness_analysis --ML_model=XGB --top_features_to_plot_feat_imp "1,2,5,10" --feat_imp_min_n_robustness=5

#Rscript plot_top_features.R --datasets="Rawlins_fetal_lung,Tsankov" --cancer_types='blum_50_top_10_perc' --cell_number_filter=100 --ML_model=XGB --annotation=finalized_annotation --top_features_to_plot="1,2,5,10,15" --feature_importance_method=permutation_importance  --robustness_seed_range=1-100 --robustness_analysis --top_features_to_plot_feat_imp="1,2,5" --skip_seeds_robustness=10,70


#Rscript plot_top_features.R --cancer_types="Skin-Melanoma,Liver-HCC,ColoRect-AdenoCA,multipe_myeloma,Eso-AdenoCA,CNS-GBM,Lung-AdenoCA,Lung-SCC" --ML_model=XGB --annotation=finalized_annotation --robustness_analysis --seed_range=1-10 --feature_importance_method=permutation_importance --folds_for_test_set=1-10 --grid_analysis --top_features_to_plot=1 --grid_cell_types="skin Melanocyte BR,liver Hepatoblasts SH,normal_colon Stem GL_Co,BR,bonemarrow B GL_BlBm,stomach Foveolar Cell BR,cerebrum Astrocytes-Oligodendrocytes SH,lung AT2 TS,lung Basal TS,mammary_tissue Basal Epithelial (Mammary) BR"

Rscript plot_top_features.R --cancer_types="Skin-Melanoma,Liver-HCC,ColoRect-AdenoCA,multiple_myeloma,Eso-AdenoCA,CNS-GBM,Lung-AdenoCA,Lung-SCC" --ML_model=XGB --annotation=finalized_annotation --robustness_analysis --seed_range=1-10 --feature_importance_method=permutation_importance --folds_for_test_set=1-10 --grid_analysis --top_features_to_plot=1 --grid_cell_types="skin_sun_exposed Melanocyte BR,liver Hepatoblasts SH,normal_colon Stem GL_Co,bonemarrow B GL_BlBm,stomach Goblet cells SH,cerebrum Astrocytes-Oligodendrocytes SH,lung AT2 TS,lung Basal TS"
#Rscript plot_top_features.R --cancer_types="Breast-AdenoCa,Lymph-BNHL,Myeloid-AML,Bone-Leiomyo,Thy-AdenoCA,Uterus-AdenoCA" --ML_model=XGB --annotation=finalized_annotation --robustness_analysis --seed_range=1-10 --feature_importance_method=permutation_importance --folds_for_test_set=1-10 --grid_analysis --top_features_to_plot=1 --grid_cell_types="mammary_tissue Basal Epithelial (Mammary) BR,bonemarrow B GL_BlBm,bonemarrow Early.Baso GL_BlBm,stomach Stromal cells SH,thyroid Thyroid Follicular Cell BR,placenta PAEP_MECOM positive cells SH"


#Rscript plot_top_features.R --cancer_types="Lymph-BNHL" --ML_model=XGB --annotation=finalized_annotation --robustness_analysis --seed_range=1-10 --feature_importance_method=permutation_importance --folds_for_test_set=1-10 --grid_analysis --top_features_to_plot=1 --grid_cell_types="mammary_tissue Basal Epithelial (Mammary) BR,bonemarrow B GL_BlBm,bonemarrow,stomach Stromal cells SH,thyroid Thyroid Follicular Cell BR"

#Rscript plot_top_features.R --cancer_types="Myeloid-AML" --ML_model=XGB --annotation=finalized_annotation --robustness_analysis --seed_range=1-10 --feature_importance_method=permutation_importance --folds_for_test_set=1-10 --grid_analysis --top_features_to_plot=1 --grid_cell_types="mammary_tissue Basal Epithelial (Mammary) BR,bonemarrow B GL_BlBm,bonemarrow,stomach Stromal cells SH,thyroid Thyroid Follicular Cell BR"

#Rscript plot_top_features.R --cancer_types="SoftTissue-Leiomyo" --ML_model=XGB --annotation=finalized_annotation --robustness_analysis --seed_range=1-10 --feature_importance_method=permutation_importance --folds_for_test_set=1-10 --grid_analysis --top_features_to_plot=1 --grid_cell_types="mammary_tissue Basal Epithelial (Mammary) BR,bonemarrow B GL_BlBm,bonemarrow,stomach Stromal cells SH,thyroid Thyroid Follicular Cell BR"

#Rscript plot_top_features.R --cancer_types="Thy-AdenoCA" --ML_model=XGB --annotation=finalized_annotation --robustness_analysis --seed_range=1-10 --feature_importance_method=permutation_importance --folds_for_test_set=1-10 --grid_analysis --top_features_to_plot=1 --grid_cell_types="mammary_tissue Basal Epithelial (Mammary) BR,bonemarrow B GL_BlBm,bonemarrow,stomach Stromal cells SH,thyroid Thyroid Follicular Cell BR"







#Rscript plot_top_features.R --top_features_to_plot_feat_imp="1,2,5" --robustness_analysis --datasets="Bingren,Greenleaf_brain,Greenleaf_pbmc_bm,Rawlins_fetal_lung,Shendure,Tsankov,Yang_kidney" --cancer_types=Lung-SCC,Lung-AdenoCA --cell_number_filter=100 --ML_model=XGB --annotation=finalized_annotation --top_features_to_plot="1,2,5,10,15" --robustness_seed_range=1-100 --feature_importance_method=permutation_importance
