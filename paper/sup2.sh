### Extended Data Figure 1A ###
#python3 ../analysis/ML/prep_ML_model_scripts.py --add_perf_to_file --cancer_types Lung-AdenoCA --scATAC_cell_number_filter 100 --annotation_dir finalized_annotation --datasets Bingren Shendure Greenleaf_colon Greenleaf_pbmc_bm Tsankov Rawlins_fetal_lung Yang_kidney --seed_interval=1-10 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method=permutation_importance --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --woo_pcawg --mem_per_core 1 --cell_types_keep="lung Neuroendocrine-Tsankov" --submit_jobs --top_features_to_plot_feat_imp 10 5 2
## run robustness plotting as explained in README

#python3 ../analysis/ML/prep_ML_model_scripts.py --add_perf_to_file --cancer_types Lung-SCC --scATAC_cell_number_filter 100 --annotation_dir finalized_annotation --datasets Bingren Shendure Greenleaf_colon Greenleaf_pbmc_bm Tsankov Rawlins_fetal_lung Yang_kidney --seed_interval=1-10 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method=permutation_importance --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --woo_pcawg --mem_per_core 1 --cell_types_keep="lung Neuroendocrine-Tsankov" --submit_jobs --top_features_to_plot_feat_imp 10 5 2 
## run robustness plotting as explained in README

#python3 ../analysis/ML/prep_ML_model_scripts.py --add_perf_to_file --cancer_types epithelioid_waddell --scATAC_cell_number_filter 100 --annotation_dir finalized_annotation --datasets Bingren Shendure Greenleaf_colon Greenleaf_pbmc_bm Tsankov Rawlins_fetal_lung Yang_kidney --seed_interval=1-10 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method=permutation_importance --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --meso --mem_per_core 1 --cell_types_keep="lung Neuroendocrine-Tsankov" --top_features_to_plot_feat_imp 10 5 2 --submit_jobs
## run robustness plotting as explained in README

#python3 ../analysis/ML/prep_ML_model_scripts.py --add_perf_to_file --cancer_types SCLC --scATAC_cell_number_filter 100 --annotation_dir finalized_annotation --datasets Bingren Shendure Greenleaf_colon Greenleaf_pbmc_bm Tsankov Rawlins_fetal_lung Yang_kidney --seed_interval=1-10 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method=permutation_importance --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --SCLC --mem_per_core 1 --cell_types_keep="lung Neuroendocrine-Tsankov" --top_features_to_plot_feat_imp 10 5 2 --submit_jobs
## run robustness plotting as explained in README

### Extended Data Figure 1B ###
TODO
