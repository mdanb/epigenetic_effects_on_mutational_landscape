### Extended Data Figure 4 ###

#python3 ../analysis/ML/prep_ML_model_scripts.py --add_perf_to_file --cancer_types Kidney-ChRCC --datasets Bingren Shendure Greenleaf_colon Greenleaf_pbmc_bm Tsankov Yang_kidney --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=8 --feature_importance_method=permutation_importance --test_set_perf_num_features all --top_features_to_plot_feat_imp 10 5 2 --woo_pcawg --submit_jobs
## run robustness plotting as explained in README

#python3 ../analysis/ML/prep_ML_model_scripts.py --add_perf_to_file --cancer_types kidney_clear_cell --datasets Bingren Shendure Greenleaf_colon Greenleaf_pbmc_bm Tsankov Yang_kidney --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=8 --feature_importance_method=permutation_importance --test_set_perf_num_features all --top_features_to_plot_feat_imp 10 5 2 --histologically_subtyped_mutations --submit_jobs
## run robustness plotting as explained in README

#python3 ../analysis/ML/prep_ML_model_scripts.py --add_perf_to_file --cancer_types kidney_papillary --datasets Bingren Shendure Greenleaf_colon Greenleaf_pbmc_bm Tsankov Yang_kidney --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=8 --feature_importance_method=permutation_importance --test_set_perf_num_features all --top_features_to_plot_feat_imp 10 5 2 --histologically_subtyped_mutations --submit_jobs
## run robustness plotting as explained in README

#python3 ../analysis/ML/prep_ML_model_scripts.py --add_perf_to_file --cancer_types Panc-AdenoCA --datasets Bingren Shendure Greenleaf_colon Greenleaf_pbmc_bm Tsankov Yang_kidney --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=8 --feature_importance_method=permutation_importance --test_set_perf_num_features all --top_features_to_plot_feat_imp 10 5 2 --woo_pcawg --submit_jobs
## run robustness plotting as explained in README

#python3 ../analysis/ML/prep_ML_model_scripts.py --add_perf_to_file --cancer_types Panc-Endocrine --datasets Bingren Shendure Greenleaf_colon Greenleaf_pbmc_bm Tsankov Yang_kidney --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=8 --feature_importance_method=permutation_importance --test_set_perf_num_features all --top_features_to_plot_feat_imp 10 5 2 --woo_pcawg --submit_jobs
## run robustness plotting as explained in README

#python3 ../analysis/ML/prep_ML_model_scripts.py --add_perf_to_file --cancer_types Panc-Endocrine --datasets Bingren Shendure --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=8 --feature_importance_method=permutation_importance --test_set_perf_num_features all --top_features_to_plot_feat_imp 10 5 2 --woo_pcawg --submit_jobs
## run robustness plotting as explained in README

#python3 ../analysis/ML/prep_ML_model_scripts.py --add_perf_to_file --cancer_types msi_high --datasets Bingren Shendure Greenleaf_colon Greenleaf_pbmc_bm Tsankov Yang_kidney --scATAC_cell_number_filter=100 --annotation_dir=Greenleaf_colon_remove_cancer_merge_normal_unaffected --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=8 --feature_importance_method=permutation_importance --test_set_perf_num_features all --top_features_to_plot_feat_imp 10 5 2 --msi_high --submit_jobs
## run robustness plotting as explained in README

