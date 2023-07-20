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


#python3 build_ML_model.py --datasets Bingren Shendure Greenleaf_brain Tsankov Yang_kidney --scATAC_cell_number_filter=1 --annotation_dir=finalized_annotation --ML_model=XGB --seed_range=42-42 --feature_importance_method=permutation_importance --top_features_to_plot 2 5 10 15 20 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --cancer_types Pancreas_RNA_Endocrine Pancreas_RNA_Fibroblast Pancreas_RNA_Immune Pancreas_RNA_Neuroendocrine Pancreas_RNA_Pancreatic_Progenitor Pancreas_RNA_Squamous --RNA_subtyped

#python3 build_ML_model.py --make_plots --de_novo_seurat_clustering --datasets Bingren Shendure Greenleaf_brain Tsankov Yang_kidney --scATAC_cell_number_filter=1 --annotation_dir=finalized_annotation --ML_model=XGB --seed_range=42-42 --feature_importance_method=permutation_importance --top_features_to_plot 2 5 10 15 20 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --cancer_types ICGC_pancreas_no_missing_subtypesxtop_2000-mut_count_thresh_0-scale_TRUE-norm_by_mut_counts_TRUE-pca_n_components_30-dims_for_UMAP_1to15-regress_by_counts_FALSE-harmony_correction_FALSE_clustering_res_0.5xcluster_0 ICGC_pancreas_no_missing_subtypesxtop_2000-mut_count_thresh_0-scale_TRUE-norm_by_mut_counts_TRUE-pca_n_components_30-dims_for_UMAP_1to15-regress_by_counts_FALSE-harmony_correction_FALSE_clustering_res_0.5xcluster_1
