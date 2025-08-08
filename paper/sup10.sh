### Sup fig 10c ###
#python3 ../analysis/ML/prep_ML_model_scripts.py --cancer_types lung_adeno_cptac --scATAC_cell_number_filter 100 --annotation_dir finalized_annotation --datasets Bingren Shendure Greenleaf_colon Greenleaf_pbmc_bm Tsankov Rawlins_fetal_lung Yang_kidney --seed_interval=1-10 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method=permutation_importance --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --CPTAC --mem_per_core 1000 --cell_types_keep="lung Neuroendocrine-Tsankov" --top_features_to_plot_feat_imp 10 5 2  --add_p_to_file --add_perf_to_file 
#python3 ../analysis/ML/prep_ML_model_scripts.py --cancer_types lung_squamous_cptac --scATAC_cell_number_filter 100 --annotation_dir finalized_annotation --datasets Bingren Shendure Greenleaf_colon Greenleaf_pbmc_bm Tsankov Rawlins_fetal_lung Yang_kidney --seed_interval=1-10 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method=permutation_importance --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --CPTAC --mem_per_core 1000 --cell_types_keep="lung Neuroendocrine-Tsankov" --top_features_to_plot_feat_imp 10 5 2  --add_p_to_file --add_perf_to_file
#python3 ../analysis/ML/prep_ML_model_scripts.py --cancer_types kidney_rcc_cptac --datasets Bingren Shendure Greenleaf_colon Greenleaf_pbmc_bm Tsankov Yang_kidney --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=8 --feature_importance_method=permutation_importance --test_set_perf_num_features all --top_features_to_plot_feat_imp 10 5 2 --CPTAC --add_p_to_file --add_perf_to_file  --mem_per_core=1000 
#python3 ../analysis/ML/prep_ML_model_scripts.py --cancer_types pancreas_infiltrating_duct_cptac --datasets Bingren Shendure Greenleaf_colon Greenleaf_pbmc_bm Tsankov Yang_kidney --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=8 --feature_importance_method=permutation_importance --test_set_perf_num_features all --top_features_to_plot_feat_imp 10 5 2 --CPTAC --add_p_to_file --add_perf_to_file  --mem_per_core=1000
#python3 ../analysis/ML/prep_ML_model_scripts.py --cancer_types gbm_cptac --scATAC_cell_number_filter 100 --annotation_dir finalized_annotation --datasets Bingren Bingren_adult_brain Shendure Greenleaf_brain Greenleaf_pbmc_bm Greenleaf_colon Tsankov Yang_kidney --seed_interval=1-10 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --CPTAC --mem_per_core 1000 --tissues_to_consider=all --add_p_to_file --add_perf_to_file 


### Sup fig 10d ###

# python3 ../analysis/ML/prep_ML_model_scripts.py --cancer_types HTCFR --datasets Bingren Greenleaf_pbmc_bm Greenleaf_colon Shendure Tsankov Yang_kidney --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --dataset_abbrev D2 \
#     --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 \
#     --cores=10 --feature_importance_method=permutation_importance --test_set_perf_num_features all --top_features_to_plot_feat_imp 5 \
#     --tissues_to_consider=all --custom_mutations --add_p_to_file --add_perf_to_file --mem_per_core=1000 --submit_jobs



# python3 ../analysis/ML/prep_ML_model_scripts.py --cancer_types HTCJP --datasets Bingren Greenleaf_pbmc_bm Greenleaf_colon Shendure Tsankov Yang_kidney --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --dataset_abbrev D2 \
#     --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 \
#     --cores=10 --feature_importance_method=permutation_importance --test_set_perf_num_features all --top_features_to_plot_feat_imp 5 \
#     --tissues_to_consider=all --custom_mutations --add_p_to_file --add_perf_to_file --mem_per_core=1000 --submit_jobs

# python3 ../analysis/ML/prep_ML_model_scripts.py --cancer_types HTCUS --datasets Bingren Greenleaf_pbmc_bm Greenleaf_colon Shendure Tsankov Yang_kidney --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --dataset_abbrev D2 \
#     --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 \
#     --cores=10 --feature_importance_method=permutation_importance --test_set_perf_num_features all --top_features_to_plot_feat_imp 5 \
#     --tissues_to_consider=all --custom_mutations --add_p_to_file --add_perf_to_file --mem_per_core=1000 --submit_jobs





