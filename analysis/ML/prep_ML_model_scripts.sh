#python3 prep_ML_model_scripts.py --cancer_types mesomics_top_r1_90th_perc --datasets Tsankov --scATAC_cell_number_filter=30 --annotation_dir=finalized_annotation --seed_interval=1-100 --seed_interval_step=5 --iters_dont_skip 17 15 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --meso --cores=8 
#python3 prep_ML_model_scripts.py --cancer_types mesomics_top_r2_90th_perc --datasets Tsankov --scATAC_cell_number_filter=30 --annotation_dir=finalized_annotation --seed_interval=1-100 --seed_interval_step=5 --iters_dont_skip 17 15 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --meso --cores=8
#python3 prep_ML_model_scripts.py --cancer_types mesomics_top_r3_90th_perc --datasets Tsankov --scATAC_cell_number_filter=30 --annotation_dir=finalized_annotation --seed_interval=1-100 --seed_interval_step=5 --iters_dont_skip 17 15 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --meso --cores=8

#python3 prep_ML_model_scripts.py --cancer_types blum_top_10_perc --datasets Tsankov --meso --scATAC_cell_number_filter=30 --annotation_dir=finalized_annotation --seed_interval=1-100 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=2

#python3 prep_ML_model_scripts.py --cancer_types nmf_top_10_perc --datasets Tsankov --meso --scATAC_cell_number_filter=30 --annotation_dir=finalized_annotation --seed_interval=1-100 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=2

#python3 prep_ML_model_scripts.py --cancer_types Lung-SCC --datasets Tsankov --scATAC_cell_number_filter=30 --annotation_dir=finalized_annotation --seed_interval=1-100 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=2 --top_features_to_plot "19" "15" "10" "5" "2" --feature_importance_method=permutation_importance

#python3 prep_ML_model_scripts.py --cancer_types Lung-SCC --datasets Bingren Shendure Greenleaf_pbmc_bm Greenleaf_brain Rawlins_fetal_lung Tsankov Yang_kidney --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=8 --top_features_to_plot "20" "15" "10" "5" "2" --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 20 --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 20 --submit_jobs
#python3 prep_ML_model_scripts.py --cancer_types Lung-AdenoCA --datasets Bingren Shendure Greenleaf_pbmc_bm Greenleaf_brain Rawlins_fetal_lung Tsankov Yang_kidney --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=8 --top_features_to_plot "20" "15" "10" "5" "2" --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 20 --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 20 --submit_jobs

#python3 prep_ML_model_scripts.py --cancer_types Lung-AdenoCA --datasets Bingren Shendure Greenleaf_pbmc_bm Greenleaf_brain Rawlins_fetal_lung Tsankov Yang_kidney --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=8 --top_features_to_plot "20" "15" "10" "5" "2" --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 20 --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 20 --submit_jobs


#python3 prep_ML_model_scripts.py --de_novo_seurat_clustering --cancer_types sclcxtop_500-scale_TRUE-norm_by_mut_counts_TRUE-pca_n_components_30-dims_for_dim_reduction_1to15-regress_by_counts_FALSE-min_mutations_per_sample0_clustering_res_1xcluster_0 --datasets Bingren Shendure Greenleaf_pbmc_bm Greenleaf_brain Rawlins_fetal_lung Tsankov Yang_kidney --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=8 --top_features_to_plot "20" "15" "10" "5" "2" --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 20 --submit_jobs --cleanup

#python3 prep_ML_model_scripts.py --de_novo_seurat_clustering --cancer_types sclcxtop_500-scale_TRUE-norm_by_mut_counts_TRUE-pca_n_components_30-dims_for_dim_reduction_1to15-regress_by_counts_FALSE-min_mutations_per_sample0_clustering_res_1xcluster_1 --datasets Bingren Shendure Greenleaf_pbmc_bm Greenleaf_brain Rawlins_fetal_lung Tsankov Yang_kidney --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=8 --top_features_to_plot "20" "15" "10" "5" "2" --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 20 --submit_jobs --cleanup

#python3 prep_ML_model_scripts.py --de_novo_seurat_clustering --cancer_types sclcxtop_500-scale_TRUE-norm_by_mut_counts_TRUE-pca_n_components_30-dims_for_dim_reduction_1to15-regress_by_counts_FALSE-min_mutations_per_sample0_clustering_res_1xcluster_2 --datasets Bingren Shendure Greenleaf_pbmc_bm Greenleaf_brain Rawlins_fetal_lung Tsankov Yang_kidney --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=4 --top_features_to_plot "20" "15" "10" "5" "2" --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 20 --submit_jobs --cleanup



#python3 prep_ML_model_scripts.py --SCLC --cancer_types SCLC.P --datasets Rawlins_fetal_lung Tsankov --scATAC_cell_number_filter=1 --annotation_dir=finalized_annotation --seed_interval=1-5 --fold_for_test_set_range=1-5 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=4 --top_features_to_plot "20" "15" "10" "5" "2" --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 20 --submit_jobs
#python3 prep_ML_model_scripts.py --SCLC --cancer_types SCLC.A --datasets Rawlins_fetal_lung Tsankov --scATAC_cell_number_filter=1 --annotation_dir=finalized_annotation --seed_interval=1-5 --fold_for_test_set_range=1-5 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=4 --top_features_to_plot "20" "15" "10" "5" "2" --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 20 --submit_jobs

#python3 prep_ML_model_scripts.py --SCLC --cancer_types SCLC.N --datasets Rawlins_fetal_lung Tsankov --scATAC_cell_number_filter=1 --annotation_dir=finalized_annotation --seed_interval=1-5 --fold_for_test_set_range=1-5 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=4 --top_features_to_plot "20" "15" "10" "5" "2" --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 20 --submit_jobs
#python3 prep_ML_model_scripts.py --SCLC --cancer_types SCLC.Y --datasets Rawlins_fetal_lung Tsankov --scATAC_cell_number_filter=1 --annotation_dir=finalized_annotation --seed_interval=1-5 --fold_for_test_set_range=1-5 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=4 --top_features_to_plot "20" "15" "10" "5" "2" --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 20 --submit_jobs


#python3 prep_ML_model_scripts.py --histologically_subtyped_mutations --cancer_types colorectal_adeno_nos --datasets Greenleaf_colon --scATAC_cell_number_filter=100 --annotation_dir=remove_cancer_polyp_merge_normal_unaffected --seed_interval=1-5 --fold_for_test_set_range=1-5 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=4 --top_features_to_plot "10" "5" "2" --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 1 2 5 10 --submit_jobs

#python3 prep_ML_model_scripts.py --histologically_subtyped_mutations --cancer_types colorectal_adeno_mucinous --datasets Greenleaf_colon --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-5 --fold_for_test_set_range=1-5 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=4 --top_features_to_plot "10" "5" "2" --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 1 2 5 10 --submit_jobs


#python3 prep_ML_model_scripts.py --de_novo_seurat_clustering --cancer_types colorectalxtop_1000-scale_TRUE-norm_by_mut_counts_TRUE-pca_n_components_20-dims_for_dim_reduction_1to15-regress_by_counts_FALSE-min_mutations_per_sample0_clustering_res_0.85xcluster_0 --datasets Greenleaf_colon --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-5 --fold_for_test_set_range=1-5 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=4 --top_features_to_plot "10" "5" "2" --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 1 2 5 10 --submit_jobs

#python3 prep_ML_model_scripts.py --de_novo_seurat_clustering --cancer_types colorectalxtop_1000-scale_TRUE-norm_by_mut_counts_TRUE-pca_n_components_20-dims_for_dim_reduction_1to15-regress_by_counts_FALSE-min_mutations_per_sample0_clustering_res_0.85xcluster_1 --datasets Greenleaf_colon --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-5 --fold_for_test_set_range=1-5 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=4 --top_features_to_plot "10" "5" "2" --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 1 2 5 10 --submit_jobs


#python3 prep_ML_model_scripts.py --cancer_types blum_top_10_perc nmf_top_10_perc blum_50_top_10_perc --datasets Tsankov --meso --scATAC_cell_number_filter=1 --annotation_dir=finalized_annotation --seed_interval=1-100 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=4 --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 20 15 10 5 2 1

#python3 prep_ML_model_scripts.py --cancer_types blum_bottom_10_perc nmf_bottom_10_perc blum_50_bottom_10_perc --datasets Tsankov --meso --scATAC_cell_number_filter=1 --annotation_dir=finalized_annotation --seed_interval=1-100 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=4 --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 20 15 10 5 2 1

#python3 prep_ML_model_scripts.py --meso --cancer_types blum_top_10_perc nmf_top_10_perc blum_50_top_10_perc --datasets Rawlins_fetal_lung Tsankov --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=66-100 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=8 --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 20 15 10 5 2 1

#python3 prep_ML_model_scripts.py --meso --cancer_types combined_meso --datasets Bingren Shendure Greenleaf_pbmc_bm Greenleaf_brain Rawlins_fetal_lung Tsankov Yang_kidney --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=4 --top_features_to_plot "20" "15" "10" "5" "2" --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 20 --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 20 --submit_jobs

#python3 prep_ML_model_scripts.py --histologically_subtyped_mutations --cancer_types kidney_clear_cell --datasets Bingren Shendure Greenleaf_pbmc_bm Greenleaf_brain Rawlins_fetal_lung Tsankov Yang_kidney --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=4 --top_features_to_plot "20" "15" "10" "5" "2" --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 20 --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 20 --submit_jobs

#python3 prep_ML_model_scripts.py --meso --cancer_types epithelioid_waddell --datasets Tsankov --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=4 --top_features_to_plot "15" "10" "5" "2" --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 --submit_jobs




#python3 prep_ML_model_scripts.py --meso --cancer_types combined_mesomics_no_biphasic --datasets Tsankov --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-5 --fold_for_test_set_range=1-5 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=4 --top_features_to_plot "15" "10" "5" "2" --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 --submit_jobs

#python3 prep_ML_model_scripts.py --meso --cancer_types epithelioid_waddell --datasets Tsankov --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-5 --fold_for_test_set_range=1-5 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=4 --top_features_to_plot "15" "10" "5" "2" --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 --submit_jobs

#python3 prep_ML_model_scripts.py --meso --cancer_types sarcomatoid_waddell --datasets Tsankov --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-5 --fold_for_test_set_range=1-5 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=4 --top_features_to_plot "15" "10" "5" "2" --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 --submit_jobs

#python3 prep_ML_model_scripts.py --meso --cancer_types s786_s848_combined --datasets Tsankov --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-5 --fold_for_test_set_range=1-5 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=4 --top_features_to_plot "15" "10" "5" "2" --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 --submit_jobs

#python3 prep_ML_model_scripts.py --meso --cancer_types all_meso_datasets_combined --datasets Tsankov --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-5 --fold_for_test_set_range=1-5 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=4 --top_features_to_plot "15" "10" "5" "2" --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 --submit_jobs

#python3 prep_ML_model_scripts.py --meso --cancer_types all_meso_datasets_combined_no_biphasic --datasets Tsankov --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-5 --fold_for_test_set_range=1-5 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=4 --top_features_to_plot "15" "10" "5" "2" --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 --submit_jobs

#python3 prep_ML_model_scripts.py --meso --cancer_types all_sarcomatoid_combined --datasets Tsankov --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-5 --fold_for_test_set_range=1-5 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=4 --top_features_to_plot "15" "10" "5" "2" --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 --submit_jobs

#python3 prep_ML_model_scripts.py --meso --cancer_types all_epithelioid_combined --datasets Tsankov --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-5 --fold_for_test_set_range=1-5 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=4 --top_features_to_plot "15" "10" "5" "2" --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 --submit_jobs

#python3 prep_ML_model_scripts.py --meso --cancer_types s786 --datasets Tsankov --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-5 --fold_for_test_set_range=1-5 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=4 --top_features_to_plot "15" "10" "5" "2" --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 --submit_jobs

#python3 prep_ML_model_scripts.py --meso --cancer_types s848 --datasets Tsankov --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-5 --fold_for_test_set_range=1-5 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=4 --top_features_to_plot "15" "10" "5" "2" --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 --submit_jobs



#python3 prep_ML_model_scripts.py --histologically_subtyped_mutations --cancer_types melanoma_superificial_spreading --datasets Bingren Shendure Greenleaf_pbmc_bm Greenleaf_brain Rawlins_fetal_lung Tsankov Yang_kidney --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=8 --top_features_to_plot "20" "15" "10" "5" "2" --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 20 --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 20 --submit_jobs

#python3 prep_ML_model_scripts.py --histologically_subtyped_mutations --cancer_types melanoma_acral_lentiginuous --datasets Bingren Shendure Greenleaf_pbmc_bm Greenleaf_brain Rawlins_fetal_lung Tsankov Yang_kidney --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=8 --top_features_to_plot "20" "15" "10" "5" "2" --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 20 --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 20 --submit_jobs

#python3 prep_ML_model_scripts.py --histologically_subtyped_mutations --cancer_types melanoma_desmoplastic --datasets Bingren Shendure Greenleaf_pbmc_bm Greenleaf_brain Rawlins_fetal_lung Tsankov Yang_kidney --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=8 --top_features_to_plot "20" "15" "10" "5" "2" --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 20 --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 20 --submit_jobs

#python3 prep_ML_model_scripts.py --histologically_subtyped_mutations --cancer_types melanoma_nodular --datasets Bingren Shendure Greenleaf_pbmc_bm Greenleaf_brain Rawlins_fetal_lung Tsankov Yang_kidney --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=8 --top_features_to_plot "20" "15" "10" "5" "2" --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 20 --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 20 --submit_jobs

#python3 prep_ML_model_scripts.py --histologically_subtyped_mutations --cancer_types melanoma_mucosal_lentiginous --datasets Bingren Shendure Greenleaf_pbmc_bm Greenleaf_brain Rawlins_fetal_lung Tsankov Yang_kidney --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=8 --top_features_to_plot "20" "15" "10" "5" "2" --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 20 --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 20 --submit_jobs


#python3 prep_ML_model_scripts.py --save_test_set_perf --cancer_types ColoRect-AdenoCA --datasets Bingren Greenleaf_brain Greenleaf_colon Greenleaf_pbmc_bm Rawlins_fetal_lung Shendure Tsankov Yang_kidney --scATAC_cell_number_filter 100 --annotation_dir remove_cancer_polyp_merge_normal_unaffected --seed_interval=1-10 --top_features_to_plot 15 10 5 2 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method=permutation_importance --fold_for_test_set_range 1-10 --test_set_perf_num_features 1 2 5 10 15 20 --cores=8 --seed_interval_step=5 --submit_jobs 


#python3 prep_ML_model_scripts.py --SCLC --cancer_types SCLC --datasets Bingren Shendure Greenleaf_pbmc_bm Greenleaf_brain Greenleaf_colon Rawlins_fetal_lung Tsankov Yang_kidney --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=8 --top_features_to_plot "20" "15" "10" "5" "2" --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 20 --cleanup --submit_jobs

#python3 prep_ML_model_scripts.py --meso --cancer_types sarcomatoid_mesomics --datasets Tsankov --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=4 --top_features_to_plot "15" "10" "5" "2" --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 --submit_jobs --cleanup

#python3 prep_ML_model_scripts.py --meso --cancer_types sarcomatoid_mesomics --datasets Tsankov --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=4 --top_features_to_plot ""15" "10" "5" "2" --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 --submit_jobs

#python3 prep_ML_model_scripts.py --meso --cancer_types biphasic_mesomics --datasets Tsankov --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=4 --top_features_to_plot "15" "10" "5" "2" --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 --submit_jobs

#python3 prep_ML_model_scripts.py --save_test_set_perf --cancer_types Skin-Melanoma --datasets Bingren --scATAC_cell_number_filter 100 --annotation_dir Bingren_remove_same_celltype_indexing --seed_interval=1-10 --top_features_to_plot 15 10 5 2 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method=permutation_importance --fold_for_test_set_range 1-10 --test_set_perf_num_features 1 2 5 10 15 20 --cores=8 --seed_interval_step=5 --hundred_kb --submit_jobs


#for i in 1 $(seq 5 5 35); do
#	python3 prep_ML_model_scripts.py --tissues_to_consider adult_brain frontal_cortex cerebrum brain cerebellum --cancer_types CNS-GBM_n_$i --scATAC_cell_number_filter 100 --annotation_dir finalized_annotation --datasets Bingren Bingren_adult_brain Greenleaf_brain Shendure --seed_interval=1-10 --top_features_to_plot 1 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method=permutation_importance --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --subsampled_mutations --mem_per_core 1 --test_set_perf_num_features all --top_features_to_plot_feat_imp 1 2 5 10 --create_bash_scripts --submit_jobs
#done

#for i in 1 $(seq 5 5 40); do
#        python3 prep_ML_model_scripts.py --tissues_to_consider all --cancer_types Kidney-ChRCC_n_$i --scATAC_cell_number_filter 100 --annotation_dir finalized_annotation --datasets Bingren Greenleaf_pbmc_bm Greenleaf_colon Shendure Tsankov Yang_kidney --seed_interval=1-10 --top_features_to_plot 1 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method=permutation_importance --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --subsampled_mutations --mem_per_core 1 --test_set_perf_num_features all --top_features_to_plot_feat_imp 1 2 5 10 --create_bash_script --submit_jobs
#done

#for i in 1 $(seq 5 5 35); do
#        python3 prep_ML_model_scripts.py --tissues_to_consider lung fetal_lung --cancer_types Lung-SCC_n_$i --scATAC_cell_number_filter 100 --annotation_dir finalized_annotation --datasets Bingren Shendure Rawlins_fetal_lung Tsankov --seed_interval=1-10 --top_features_to_plot 1 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method=permutation_importance --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --subsampled_mutations --mem_per_core 1 --test_set_perf_num_features all --top_features_to_plot_feat_imp 1 2 5 10 --create_bash_scripts --cell_types_keep="lung Neuroendocrine-Tsankov" --submit_jobs
#done

python3 prep_ML_model_scripts.py --tissues_to_consider colon_sigmoid colon_transverse intestine normal_colon polyp_colon --cancer_types mss --scATAC_cell_number_filter 100 --annotation_dir Greenleaf_colon_remove_cancer_merge_normal_unaffected --datasets Bingren Greenleaf_colon Shendure --seed_interval=1-10 --top_features_to_plot 1 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method=permutation_importance --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --msi_high --mem_per_core 1 --test_set_perf_num_features all --top_features_to_plot_feat_imp 1 2 5 10 --create_bash_scripts --submit_jobs
#python3 prep_ML_model_scripts.py --cancer_types Lymph-BNHL --scATAC_cell_number_filter 100 --annotation_dir finalized_annotation --datasets Greenleaf_pbmc_bm Shendure Bingren  --seed_interval=1-10 --top_features_to_plot 1 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method=permutation_importance --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --woo_pcawg --mem_per_core 1 --top_features_to_plot_feat_imp 1 2 5 10 --create_bash_scripts --grid_analysis --grid_cell_types "bonemarrow B GL_BlBm,stomach Stromal cells SH,placenta PAEP_MECOM positive cells SH,thyroid Thyroid Follicular Cell BR"

#python3 prep_ML_model_scripts.py --cancer_types Bone-Leiomyo --scATAC_cell_number_filter 100 --annotation_dir finalized_annotation --datasets Greenleaf_pbmc_bm Shendure Bingren  --seed_interval=1-10 --top_features_to_plot 1 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method=permutation_importance --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --woo_pcawg --mem_per_core 1 --top_features_to_plot_feat_imp 1 2 5 10 --create_bash_scripts --grid_analysis --grid_cell_types "bonemarrow B GL_BlBm,stomach Stromal cells SH,placenta PAEP_MECOM positive cells SH,thyroid Thyroid Follicular Cell BR"
#python3 prep_ML_model_scripts.py --cancer_types Thy-AdenoCA --scATAC_cell_number_filter 100 --annotation_dir finalized_annotation --datasets Greenleaf_pbmc_bm Shendure Bingren  --seed_interval=1-10 --top_features_to_plot 1 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method=permutation_importance --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --woo_pcawg --mem_per_core 1 --top_features_to_plot_feat_imp 1 2 5 10 --create_bash_scripts --grid_analysis --grid_cell_types "bonemarrow B GL_BlBm,stomach Stromal cells SH,placenta PAEP_MECOM positive cells SH,thyroid Thyroid Follicular Cell BR"
#python3 prep_ML_model_scripts.py --cancer_types Uterus-AdenoCA --scATAC_cell_number_filter 100 --annotation_dir finalized_annotation --datasets Greenleaf_pbmc_bm Shendure Bingren  --seed_interval=1-10 --top_features_to_plot 1 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method=permutation_importance --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --woo_pcawg --mem_per_core 1 --top_features_to_plot_feat_imp 1 2 5 10 --create_bash_scripts --grid_analysis --grid_cell_types "bonemarrow B GL_BlBm,stomach Stromal cells SH,placenta PAEP_MECOM positive cells SH,thyroid Thyroid Follicular Cell BR"




#python3 prep_ML_model_scripts.py --cancer_types Melanomaxtop_2000-scale_TRUE-norm_by_mut_counts_TRUE-pca_n_components_30-dims_for_dim_reduction_1to15-regress_by_counts_FALSE-min_mutations_per_sample0_clustering_res_0.6xcluster_0 Melanomaxtop_2000-scale_TRUE-norm_by_mut_counts_TRUE-pca_n_components_30-dims_for_dim_reduction_1to15-regress_by_counts_FALSE-min_mutations_per_sample0_clustering_res_0.6xcluster_1 --datasets Bingren Shendure Greenleaf_pbmc_bm Greenleaf_brain Rawlins_fetal_lung Tsankov Yang_kidney --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-5 --fold_for_test_set_range=1-5 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=8 --top_features_to_plot "20" "15" "10" "5" "2" --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 20 --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 20 --submit_jobs
