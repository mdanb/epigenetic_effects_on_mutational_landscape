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



#python3 prep_ML_model_scripts.py --SCLC --cancer_types SCLC.P --datasets Bingren Shendure Greenleaf_pbmc_bm Greenleaf_brain Rawlins_fetal_lung Tsankov Yang_kidney --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=8 --top_features_to_plot "20" "15" "10" "5" "2" --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 20 

#python3 prep_ML_model_scripts.py --SCLC --cancer_types SCLC.A --datasets Bingren Shendure Greenleaf_pbmc_bm Greenleaf_brain Rawlins_fetal_lung Tsankov Yang_kidney --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=8 --top_features_to_plot "20" "15" "10" "5" "2" --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 20

#python3 prep_ML_model_scripts.py --SCLC --cancer_types SCLC.N --datasets Bingren Shendure Greenleaf_pbmc_bm Greenleaf_brain Rawlins_fetal_lung Tsankov Yang_kidney --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=8 --top_features_to_plot "20" "15" "10" "5" "2" --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 20

#python3 prep_ML_model_scripts.py --SCLC --cancer_types SCLC.Y --datasets Bingren Shendure Greenleaf_pbmc_bm Greenleaf_brain Rawlins_fetal_lung Tsankov Yang_kidney --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=8 --top_features_to_plot "20" "15" "10" "5" "2" --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 20


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

python3 prep_ML_model_scripts.py --histologically_subtyped_mutations --cancer_types melanoma_mucosal_lentiginous --datasets Bingren Shendure Greenleaf_pbmc_bm Greenleaf_brain Rawlins_fetal_lung Tsankov Yang_kidney --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=8 --top_features_to_plot "20" "15" "10" "5" "2" --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 20 --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 20 --submit_jobs






#python3 prep_ML_model_scripts.py --cancer_types Melanomaxtop_2000-scale_TRUE-norm_by_mut_counts_TRUE-pca_n_components_30-dims_for_dim_reduction_1to15-regress_by_counts_FALSE-min_mutations_per_sample0_clustering_res_0.6xcluster_0 Melanomaxtop_2000-scale_TRUE-norm_by_mut_counts_TRUE-pca_n_components_30-dims_for_dim_reduction_1to15-regress_by_counts_FALSE-min_mutations_per_sample0_clustering_res_0.6xcluster_1 --datasets Bingren Shendure Greenleaf_pbmc_bm Greenleaf_brain Rawlins_fetal_lung Tsankov Yang_kidney --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-5 --fold_for_test_set_range=1-5 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=8 --top_features_to_plot "20" "15" "10" "5" "2" --feature_importance_method=permutation_importance --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 20 --save_test_set_perf --test_set_perf_num_features 1 2 5 10 15 20 --submit_jobs
