#python3 build_ML_model.py --cancer_types sarcomatoid_waddell_mesomics --datasets Bingren Greenleaf_brain Shendure Yang_kidney --scATAC_cell_number_filter=30 --annotation_dir=finalized_annotation --meso --ML_model=XGB --seed=42:42 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50

#python3 build_ML_model.py --test_backward_selection_iter 17 --cancer_types sarcomatoid_waddell_mesomics --datasets Tsankov --scATAC_cell_number_filter=30 --annotation_dir=finalized_annotation --meso --ML_model=XGB --seed_range=1-100 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --iters_dont_skip=17

#python3 build_ML_model.py --cancer_types sarcomatoid_waddell_mesomics --datasets Tsankov --scATAC_cell_number_filter=30 --annotation_dir=finalized_annotation --meso --ML_model=XGB --seed_range=101-101 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --iters_dont_skip=17

#python3 build_ML_model.py --test_backward_selection_iter 15 --cancer_types sarcomatoid_waddell_mesomics --datasets Tsankov --scATAC_cell_number_filter=30 --annotation_dir=finalized_annotation --meso --ML_model=XGB --seed_range=1-100 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --iters_dont_skip=17

#python3 build_ML_model.py --cancer_types SCLC --datasets Greenleaf_brain Bingren Shendure Tsankov --scATAC_cell_number_filter=30 --annotation_dir=finalized_annotation --ML_model=XGB --seed=42:42 --SCLC --n_optuna_trials_prebackward_selection=0 --n_optuna_trials_backward_selection=0

#python3 build_ML_model.py --cancer_types Lung-AdenoCA Lung-SCC --datasets Tsankov --scATAC_cell_number_filter=30 --annotation_dir=finalized_annotation --ML_model=XGB --seed_range=1-1 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --iters_dont_skip=17

#python3 /broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML/build_ML_model.py --cancer_types Lung-AdenoCA --datasets Tsankov --scATAC_cell_number_filter 30 --annotation_dir finalized_annotation --seed_range 6-10 --iters_dont_skip 17 15 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50

#python3 /broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML/build_ML_model.py --cancer_types Lung-SCC --datasets Tsankov --scATAC_cell_number_filter 30 --annotation_dir finalized_annotation --seed_range 6-10 --iters_dont_skip 17 15 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50

#python3 /broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML/build_ML_model.py --cancer_types Lung-AdenoCA --datasets Tsankov --scATAC_cell_number_filter 30 --annotation_dir finalized_annotation --seed_range 6-10 --iters_dont_skip 17 15 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50

#python3 /broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML/build_ML_model.py --cancer_types Lung-AdenoCA --datasets Tsankov --scATAC_cell_number_filter 30 --annotation_dir finalized_annotation --seed_range 1-5 --iters_dont_skip 17 15 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50

#python3 /broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML/build_ML_model.py --cancer_types Lung-SCC --datasets Tsankov --scATAC_cell_number_filter 30 --annotation_dir finalized_annotation --seed_range 1-100 --top_features_to_plot 2 5 10 15 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50


##########

#python3 build_ML_model.py --cancer_types CNS-GBM CNS-Medullo CNS-PiloAstro CNS-Oligo --datasets Bingren Shendure Greenleaf_brain Tsankov Yang_kidney --scATAC_cell_number_filter=1 --annotation_dir=finalized_annotation --ML_model=XGB --seed_range=42-42 --feature_importance_method=permutation_importance --make_plots --top_features_to_plot 2 5 10 15 20 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50

#python3 build_ML_model.py --histologically_subtyped_mutations --cancer_types kidney_papillary kidney_clear_cell --datasets Bingren Shendure Greenleaf_brain Tsankov Yang_kidney --scATAC_cell_number_filter=1 --annotation_dir=finalized_annotation --ML_model=XGB --seed_range=42-42 --feature_importance_method=permutation_importance --make_plots --top_features_to_plot 2 5 10 15 20 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50

#python3 build_ML_model.py --cancer_types Kidney-ChRCC --datasets Bingren Shendure Greenleaf_brain Tsankov Yang_kidney --scATAC_cell_number_filter=1 --annotation_dir=finalized_annotation --ML_model=XGB --seed_range=42-42 --feature_importance_method=permutation_importance --make_plots --top_features_to_plot 2 5 10 15 20 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50

#python3 build_ML_model.py --cancer_types medullo_desmoplastic medullo_large_cell pancreatic_acinar_cell pancreatic_mucinuous --histologically_subtyped_mutations --datasets Bingren Shendure Greenleaf_brain Tsankov Yang_kidney --scATAC_cell_number_filter=1 --annotation_dir=finalized_annotation --ML_model=XGB --seed_range=42-42 --feature_importance_method=permutation_importance --make_plots --top_features_to_plot 2 5 10 15 20 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50

#python3 build_ML_model.py --make_plots --cancer_types pancreatic_adenosquamous pancreatic_invasive_IPMN pancreatic_ductal --histologically_subtyped_mutations --datasets Bingren Shendure Greenleaf_brain Tsankov Yang_kidney --scATAC_cell_number_filter=1 --annotation_dir=finalized_annotation --ML_model=XGB --seed_range=42-42 --feature_importance_method=permutation_importance --make_plots --top_features_to_plot 2 5 10 15 20 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50

#python3 build_ML_model.py --meso --cancer_types blum_top_10_perc blum_bottom_10_perc nmf_top_10_perc nmf_bottom_10_perc blum_50_top_10_perc blum_50_bottom_10_perc --datasets Tsankov --scATAC_cell_number_filter=1 --annotation_dir=finalized_annotation --ML_model=XGB --seed_range=42-42 --feature_importance_method=permutation_importance --make_plots --top_features_to_plot 2 5 10 15 20 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50

#python3 build_ML_model.py --make_plots --datasets Bingren Shendure Greenleaf_brain Tsankov Yang_kidney --scATAC_cell_number_filter=1 --annotation_dir=finalized_annotation --ML_model=XGB --seed_range=42-42 --feature_importance_method=permutation_importance --top_features_to_plot 2 5 10 15 20 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --cancer_types Pancreas_RNA_Endocrine Pancreas_RNA_Fibroblast Pancreas_RNA_Immune Pancreas_RNA_Neuroendocrine Pancreas_RNA_Pancreatic_Progenitor Pancreas_RNA_Squamous --RNA_subtyped

#python3 build_ML_model.py --make_plots --de_novo_seurat_clustering --datasets Bingren Shendure Greenleaf_brain Tsankov Yang_kidney --scATAC_cell_number_filter=1 --annotation_dir=finalized_annotation --ML_model=XGB --seed_range=42-42 --feature_importance_method=permutation_importance --top_features_to_plot 2 5 10 15 20 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --cancer_types ICGC_pancreas_no_missing_subtypesxtop_2000-mut_count_thresh_0-scale_TRUE-norm_by_mut_counts_TRUE-pca_n_components_30-dims_for_UMAP_1to15-regress_by_counts_FALSE-harmony_correction_FALSE_clustering_res_0.5xcluster_0 ICGC_pancreas_no_missing_subtypesxtop_2000-mut_count_thresh_0-scale_TRUE-norm_by_mut_counts_TRUE-pca_n_components_30-dims_for_UMAP_1to15-regress_by_counts_FALSE-harmony_correction_FALSE_clustering_res_0.5xcluster_1

#python3 build_ML_model.py --save_test_set_perf --cancer_types Lung-SCC --datasets Tsankov --scATAC_cell_number_filter=1 --annotation_dir=finalized_annotation --ML_model=XGB --seed_range=1-5 --feature_importance_method=permutation_importance --test_set_perf_num_features 1 2 5 10 15 20 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50
#python3 build_ML_model.py --meso --save_test_set_perf --cancer_types blum_bottom_10_perc --datasets Tsankov --scATAC_cell_number_filter=1 --annotation_dir=finalized_annotation --ML_model=XGB --seed_range=42-42 --feature_importance_method=permutation_importance --test_set_perf_num_features 1 2 5 10 15 20 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50

#python3 build_ML_model.py --make_plots --save_test_set_perf --cancer_types Biliary-AdenoCA Bladder-TCC Breast-AdenoCA Breast-LobularCA --datasets Bingren Shendure Greenleaf_brain Tsankov Yang_kidney --scATAC_cell_number_filter=1 --annotation_dir=finalized_annotation --ML_model=XGB --seed_range=42-42 --feature_importance_method=permutation_importance --test_set_perf_num_features 1 2 5 10 15 20 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --top_features_to_plot 2 5 10 15

#python3 build_ML_model.py --make_plots --save_test_set_perf --cancer_types Ovary-AdenoCA Prost-AdenoCA Skin-Melanoma Thy-AdenoCA Uterus-AdenoCA --datasets Bingren Shendure Greenleaf_brain Tsankov Yang_kidney --scATAC_cell_number_filter=1 --annotation_dir=finalized_annotation --ML_model=XGB --seed_range=42-42 --feature_importance_method=permutation_importance --test_set_perf_num_features 1 2 5 10 15 20 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --top_features_to_plot 2 5 10 15

#python3 build_ML_model.py --save_test_set_perf --cancer_types Liver-HCC Lymph-BNHL Lymph-CLL Myeloid-MDS Myeloid-MPN Myeloid-AML --datasets Bingren Shendure Greenleaf_pbmc_bm Greenleaf_brain Tsankov Yang_kidney --scATAC_cell_number_filter=1 --annotation_dir=finalized_annotation --ML_model=XGB --seed_range=42-42 --feature_importance_method=permutation_importance --test_set_perf_num_features 1 2 5 10 15 20 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50

#python3 build_ML_model.py --save_test_set_perf --cancer_types Bone-Osteosarc Bone-Epith --datasets Bingren Shendure Greenleaf_pbmc_bm Greenleaf_brain Tsankov Yang_kidney --scATAC_cell_number_filter=1 --annotation_dir=finalized_annotation --ML_model=XGB --seed_range=42-42 --feature_importance_method=permutation_importance --test_set_perf_num_features 1 2 5 10 15 20 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50

#python3 build_ML_model.py --make_plots --SCLC --save_test_set_perf --cancer_types SCLC --datasets Bingren Shendure Greenleaf_pbmc_bm Greenleaf_brain Rawlins_fetal_lung Tsankov Yang_kidney --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --ML_model=XGB --seed_range=42-42 --feature_importance_method=permutation_importance --test_set_perf_num_features 1 2 5 10 15 20 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --top_features_to_plot 2 5 10 15 20

#python3 build_ML_model.py --make_plots --SCLC --save_test_set_perf --cancer_types SCLC --datasets Rawlins_fetal_lung Tsankov --scATAC_cell_number_filter=100 --annotation_dir=triangle_cells_combine_meso_no_schwann --ML_model=XGB --seed_range=1-5 --feature_importance_method=permutation_importance --test_set_perf_num_features 1 2 5 10 15 20 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --top_features_to_plot 2 5 10 15 20

#python3 /broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML/build_ML_model.py --cancer_types Lung-AdenoCA --datasets Bingren Greenleaf_brain Greenleaf_pbmc_bm Rawlins_fetal_lung Shendure Tsankov Yang_kidney --scATAC_cell_number_filter 100 --annotation_dir finalized_annotation --seed_range 1-100 --top_features_to_plot 18 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method permutation_importance --test_set_perf_num_features 20 15 10 5 2 1 --save_test_set_perf --debug_bfs

#python3 /broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML/build_ML_model.py --cancer_types Lung-SCC --datasets Bingren Greenleaf_brain Greenleaf_pbmc_bm Rawlins_fetal_lung Shendure Tsankov Yang_kidney --scATAC_cell_number_filter 100 --annotation_dir finalized_annotation --seed_range 1-100 --top_features_to_plot 18 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method permutation_importance --test_set_perf_num_features 20 15 10 5 2 1 --save_test_set_perf --debug_bfs

#python3 build_ML_model.py --make_plots --meso --save_test_set_perf --cancer_types blum_top_10_perc nmf_top_10_perc blum_50_top_10_perc --datasets Rawlins_fetal_lung Tsankov --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --ML_model=XGB --seed_range=1-5 --feature_importance_method=permutation_importance --test_set_perf_num_features 1 2 5 10 15 20 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --top_features_to_plot 2 5 10 15 20

#python3 build_ML_model.py --make_plots --meso --save_test_set_perf --cancer_types blum_bottom_10_perc nmf_bottom_10_perc blum_50_bottom_10_perc --datasets Rawlins_fetal_lung Tsankov --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --ML_model=XGB --seed_range=40-45 --feature_importance_method=permutation_importance --test_set_perf_num_features 1 2 5 10 15 20 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --top_features_to_plot 2 5 10 15 20

#python3 build_ML_model.py --make_plots --SCLC --save_test_set_perf --cancer_types SCLC --datasets Rawlins_fetal_lung Wang_lung --scATAC_cell_number_filter=100 --annotation_dir=triangle_cells_combine_meso_and_neuroendocrine_no_schwann --ML_model=XGB --seed_range=1-5 --feature_importance_method=permutation_importance --test_set_perf_num_features 1 2 5 10 15 20 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --top_features_to_plot 2 5 10 15 20

#python3 build_ML_model.py --SCLC --save_test_set_perf --cancer_types SCLC --datasets Rawlins_fetal_lung --scATAC_cell_number_filter=100 --annotation_dir=triangle_cells_combine_meso_and_neuroendocrine_no_schwann --ML_model=XGB --seed_range=1-5 --feature_importance_method=permutation_importance --test_set_perf_num_features 1 2 5 10 13 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --top_features_to_plot 2 5 10 13

#python3 build_ML_model.py --SCLC --save_test_set_perf --cancer_types SCLC --datasets Tsankov --scATAC_cell_number_filter=1 --annotation_dir=remove_basal --ML_model=XGB --seed_range=1-5 --feature_importance_method=permutation_importance --test_set_perf_num_features 1 2 5 10 15 18 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --top_features_to_plot 1 2 5 10 15

#python3 build_ML_model.py --SCLC --save_test_set_perf --cancer_types SCLC --datasets Wang_lung --scATAC_cell_number_filter=100 --annotation_dir=remove_basal --ML_model=XGB --seed_range=1-5 --feature_importance_method=permutation_importance --test_set_perf_num_features 1 2 5 10 15 18 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --top_features_to_plot 1 2 5 10 15

#python3 build_ML_model.py --SCLC --save_test_set_perf --cancer_types SCLC --datasets Wang_lung Rawlins_fetal_lung --scATAC_cell_number_filter=100 --annotation_dir=remove_basal --ML_model=XGB --seed_range=1-5 --feature_importance_method=permutation_importance --test_set_perf_num_features 1 2 5 10 15 20 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --top_features_to_plot 1 2 5 10 15

#python3 build_ML_model.py --make_plots --de_novo_seurat_clustering --save_test_set_perf --cancer_types sclcxtop_500-scale_TRUE-norm_by_mut_counts_TRUE-pca_n_components_30-dims_for_dim_reduction_1to15-regress_by_counts_FALSE-min_mutations_per_sample0_clustering_res_1xcluster_0 sclcxtop_500-scale_TRUE-norm_by_mut_counts_TRUE-pca_n_components_30-dims_for_dim_reduction_1to15-regress_by_counts_FALSE-min_mutations_per_sample0_clustering_res_1xcluster_1 sclcxtop_500-scale_TRUE-norm_by_mut_counts_TRUE-pca_n_components_30-dims_for_dim_reduction_1to15-regress_by_counts_FALSE-min_mutations_per_sample0_clustering_res_1xcluster_2 --datasets Tsankov Rawlins_fetal_lung --scATAC_cell_number_filter=1 --annotation_dir=finalized_annotation --ML_model=XGB --seed_range=1-5 --feature_importance_method=permutation_importance --test_set_perf_num_features 1 2 5 10 15 20 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --top_features_to_plot 1 2 5 10 15

#python3 build_ML_model.py --make_plots --de_novo_seurat_clustering --save_test_set_perf --cancer_types sclcxtop_500-scale_TRUE-norm_by_mut_counts_TRUE-pca_n_components_30-dims_for_dim_reduction_1to15-regress_by_counts_FALSE-min_mutations_per_sample0_clustering_res_0.7xcluster_0 sclcxtop_500-scale_TRUE-norm_by_mut_counts_TRUE-pca_n_components_30-dims_for_dim_reduction_1to15-regress_by_counts_FALSE-min_mutations_per_sample0_clustering_res_0.7xcluster_1 --datasets Tsankov Rawlins_fetal_lung --scATAC_cell_number_filter=1 --annotation_dir=finalized_annotation --ML_model=XGB --seed_range=1-5 --feature_importance_method=permutation_importance --test_set_perf_num_features 1 2 5 10 15 20 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --top_features_to_plot 1 2 5 10 15

#python3 build_ML_model.py --meso --save_test_set_perf --cancer_types combined_meso --datasets Wang_lung Rawlins_fetal_lung --scATAC_cell_number_filter=100 --annotation_dir=triangle_cells_combine_meso_and_neuroendocrine_no_schwann --ML_model=XGB --seed_range=1-5 --feature_importance_method=permutation_importance --test_set_perf_num_features 1 2 5 10 15 20 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50

#python3 build_ML_model.py --save_test_set_perf --cancer_types Lung-SCC Lung-AdenoCA --datasets Wang_lung Tsankov --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --ML_model=XGB --seed_range=1-5 --feature_importance_method=permutation_importance --test_set_perf_num_features 1 2 5 10 15 20 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --fold_for_test_set=1

#python3 build_ML_model.py --meso --save_test_set_perf --cancer_types blum_50_top_10_perc --datasets Wang_lung Rawlins_fetal_lung --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --ML_model=XGB --seed_range=10-10 --feature_importance_method=permutation_importance --test_set_perf_num_features 1 2 5 10 15 20 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50

#python3 build_ML_model.py --meso --save_test_set_perf --cancer_types nmf_top_10_perc --datasets Tsankov Rawlins_fetal_lung --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --ML_model=XGB --seed_range=70-70 --feature_importance_method=permutation_importance --test_set_perf_num_features 1 2 5 10 15 20 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50

#python3 build_ML_model.py --meso --save_test_set_perf --cancer_types nmf_top_10_perc blum_top_10_perc blum_50_top_10_perc --datasets Wang_lung --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --ML_model=XGB --seed_range=1-1 --feature_importance_method=permutation_importance --test_set_perf_num_features 1 2 5 10 15 20 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50


#python3 build_ML_model.py --de_novo_seurat_clustering --save_test_set_perf --cancer_types sclcxtop_500-scale_TRUE-norm_by_mut_counts_TRUE-pca_n_components_30-dims_for_dim_reduction_1to15-regress_by_counts_FALSE-min_mutations_per_sample0_clustering_res_1xcluster_0 sclcxtop_500-scale_TRUE-norm_by_mut_counts_TRUE-pca_n_components_30-dims_for_dim_reduction_1to15-regress_by_counts_FALSE-min_mutations_per_sample0_clustering_res_1xcluster_1 sclcxtop_500-scale_TRUE-norm_by_mut_counts_TRUE-pca_n_components_30-dims_for_dim_reduction_1to15-regress_by_counts_FALSE-min_mutations_per_sample0_clustering_res_1xcluster_2 --datasets Tsankov Rawlins_fetal_lung --scATAC_cell_number_filter=1 --annotation_dir=Tsankov_basal_refined --ML_model=XGB --seed_range=1-5 --feature_importance_method=permutation_importance --test_set_perf_num_features 1 2 5 10 15 20 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --top_features_to_plot 1 2 5 10 15

#python3 build_ML_model.py --de_novo_seurat_clustering --save_test_set_perf --cancer_types sclcxtop_500-scale_TRUE-norm_by_mut_counts_TRUE-pca_n_components_30-dims_for_dim_reduction_1to15-regress_by_counts_FALSE-min_mutations_per_sample0_clustering_res_0.7xcluster_0 sclcxtop_500-scale_TRUE-norm_by_mut_counts_TRUE-pca_n_components_30-dims_for_dim_reduction_1to15-regress_by_counts_FALSE-min_mutations_per_sample0_clustering_res_0.7xcluster_1 --datasets Tsankov Rawlins_fetal_lung --scATAC_cell_number_filter=1 --annotation_dir=Tsankov_basal_refined --ML_model=XGB --seed_range=1-5 --feature_importance_method=permutation_importance --test_set_perf_num_features 1 2 5 10 15 20 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --top_features_to_plot 1 2 5 10 15


#python3 build_ML_model.py --SCLC --save_test_set_perf --cancer_types SCLC.P --datasets Rawlins_fetal_lung Tsankov --scATAC_cell_number_filter=1 --annotation_dir=triangle_cells_combine_meso_and_neuroendocrine_no_schwann --ML_model=XGB --seed_range=1-5 --feature_importance_method=permutation_importance --test_set_perf_num_features 1 2 5 10 15 20 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --top_features_to_plot 2 5 10 15 20 --fold_for_test_set 1

#python3 build_ML_model.py --SCLC --save_test_set_perf --cancer_types SCLC.A --datasets Rawlins_fetal_lung Tsankov --scATAC_cell_number_filter=1 --annotation_dir=triangle_cells_combine_meso_and_neuroendocrine_no_schwann --ML_model=XGB --seed_range=1-5 --feature_importance_method=permutation_importance --test_set_perf_num_features 1 2 5 10 15 20 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --top_features_to_plot 2 5 10 15 20 --fold_for_test_set 1

#python3 build_ML_model.py --SCLC --save_test_set_perf --cancer_types SCLC.N --datasets Rawlins_fetal_lung Tsankov --scATAC_cell_number_filter=1 --annotation_dir=triangle_cells_combine_meso_and_neuroendocrine_no_schwann --ML_model=XGB --seed_range=1-5 --feature_importance_method=permutation_importance --test_set_perf_num_features 1 2 5 10 15 20 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --top_features_to_plot 2 5 10 15 20 --fold_for_test_set 1

#python3 build_ML_model.py --SCLC --save_test_set_perf --cancer_types SCLC.Y --datasets Rawlins_fetal_lung Tsankov --scATAC_cell_number_filter=1 --annotation_dir=triangle_cells_combine_meso_and_neuroendocrine_no_schwann --ML_model=XGB --seed_range=1-5 --feature_importance_method=permutation_importance --test_set_perf_num_features 1 2 5 10 15 20 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --top_features_to_plot 2 5 10 15 20 --fold_for_test_set 1


#python3 build_ML_model.py --SCLC --save_test_set_perf --cancer_types SCLC.P --datasets Rawlins_fetal_lung Tsankov --scATAC_cell_number_filter=1 --annotation_dir=Tsankov_basal_refined --ML_model=XGB --seed_range=1-5 --feature_importance_method=permutation_importance --test_set_perf_num_features 1 2 5 10 15 20 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --top_features_to_plot 2 5 10 15 20 --fold_for_test_set 1

#python3 build_ML_model.py --SCLC --save_test_set_perf --cancer_types SCLC.A --datasets Rawlins_fetal_lung Tsankov --scATAC_cell_number_filter=1 --annotation_dir=Tsankov_basal_refined --ML_model=XGB --seed_range=1-5 --feature_importance_method=permutation_importance --test_set_perf_num_features 1 2 5 10 15 20 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --top_features_to_plot 2 5 10 15 20 --fold_for_test_set 1

#python3 build_ML_model.py --SCLC --save_test_set_perf --cancer_types SCLC.N --datasets Rawlins_fetal_lung Tsankov --scATAC_cell_number_filter=1 --annotation_dir=Tsankov_basal_refined --ML_model=XGB --seed_range=1-5 --feature_importance_method=permutation_importance --test_set_perf_num_features 1 2 5 10 15 20 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --top_features_to_plot 2 5 10 15 20 --fold_for_test_set 1

#python3 build_ML_model.py --SCLC --save_test_set_perf --cancer_types SCLC.Y --datasets Rawlins_fetal_lung Tsankov --scATAC_cell_number_filter=1 --annotation_dir=Tsankov_basal_refined --ML_model=XGB --seed_range=1-5 --feature_importance_method=permutation_importance --test_set_perf_num_features 1 2 5 10 15 20 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --top_features_to_plot 2 5 10 15 20 --fold_for_test_set 1

python3 /broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML/build_ML_model.py --save_test_set_perf --cancer_types SCLC --datasets Bingren Greenleaf_brain Greenleaf_pbmc_bm Rawlins_fetal_lung Shendure Tsankov Yang_kidney --scATAC_cell_number_filter 100 --annotation_dir finalized_annotation --seed_range 1-1 --top_features_to_plot 20 15 10 5 2 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method permutation_importance --fold_for_test_set 10 --SCLC --test_set_perf_num_features 1 2 5 10 15 20





#python3 /broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML/build_ML_model.py --save_test_set_perf --cancer_types epithelioid_mesomics sarcomatoid_mesomics biphasic_mesomics --datasets Bingren Greenleaf_brain Greenleaf_pbmc_bm Rawlins_fetal_lung Shendure Tsankov Yang_kidney --scATAC_cell_number_filter 100 --annotation_dir finalized_annotation --seed_range 1-1 --top_features_to_plot 19 15 10 5 2 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method permutation_importance --fold_for_test_set 1 --test_set_perf_num_features 1 2 5 10 15 20 --meso


#python3 /broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML/build_ML_model.py --save_test_set_perf --cancer_types combined_meso --datasets Bingren Shendure Greenleaf_brain Greenleaf_pbmc_bm Rawlins_fetal_lung Yang_kidney --scATAC_cell_number_filter 100 --annotation_dir finalized_annotation --seed_range 10-10 --top_features_to_plot 15 10 5 2 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method permutation_importance --fold_for_test_set 10 --test_set_perf_num_features 1 2 5 10 15 --meso
#python3 /broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML/build_ML_model.py --save_test_set_perf --cancer_types combined_meso --datasets Tsankov --scATAC_cell_number_filter 100 --annotation_dir finalized_annotation --seed_range 10-10 --top_features_to_plot 19 15 10 5 2 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method permutation_importance --fold_for_test_set 1 --test_set_perf_num_features 1 2 5 10 15 20 --meso

#python3 /broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML/build_ML_model.py --cancer_types woo_hg38 --datasets Tsankov --scATAC_cell_number_filter 100 --annotation_dir finalized_annotation --seed_range 1-5 --top_features_to_plot 15 10 5 2 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method permutation_importance --fold_for_test_set 1 --meso --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15

#python3 build_ML_model.py --save_test_set_perf --cancer_types Head-SCC Cervix-AdenoCA --datasets Bingren Shendure Greenleaf_brain Tsankov Yang_kidney --scATAC_cell_number_filter=1 --annotation_dir=finalized_annotation --ML_model=XGB --seed_range=42-42 --feature_importance_method=permutation_importance --test_set_perf_num_features 1 2 5 10 15 20 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50
