### Fig 4A ###
#python3 prep_ML_model_scripts.py --cancer_types CNS-Medullo --scATAC_cell_number_filter 100 --annotation_dir finalized_annotation --datasets Bingren Bingren_adult_brain Shendure Greenleaf_brain --seed_interval=1-10 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --woo_pcawg --mem_per_core 1 --tissues_to_consider adult_brain frontal_cortex cerebrum brain cerebellum
## run robustness plotting as explained in README
#python3 prep_ML_model_scripts.py --cancer_types CNS-PiloAstro --scATAC_cell_number_filter 100 --annotation_dir finalized_annotation --datasets Bingren Bingren_adult_brain Shendure Greenleaf_brain --seed_interval=1-10 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --woo_pcawg --mem_per_core 1 --tissues_to_consider adult_brain frontal_cortex cerebrum brain cerebellum --top_features_to_plot_feat_imp 2 5 10
## run robustness plotting as explained in README
#python3 prep_ML_model_scripts.py --cancer_types CNS-GBM --scATAC_cell_number_filter 100 --annotation_dir finalized_annotation --datasets Bingren Bingren_adult_brain Shendure Greenleaf_brain --seed_interval=1-10 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --woo_pcawg --mem_per_core 1 --tissues_to_consider adult_brain frontal_cortex cerebrum brain cerebellum
## run robustness plotting as explained in README
#python3 prep_ML_model_scripts.py --cancer_types CNS-Oligo --scATAC_cell_number_filter 100 --annotation_dir finalized_annotation --datasets Bingren Bingren_adult_brain Shendure Greenleaf_brain --seed_interval=1-10 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --woo_pcawg --mem_per_core 1 --tissues_to_consider adult_brain frontal_cortex cerebrum brain cerebellum
## run robustness plotting as explained in README
### Fig 4B ###


### Fig 4C ###
#python3 prep_ML_model_scripts.py --cancer_types Lymph-BNHL --scATAC_cell_number_filter 100 --annotation_dir new_intermediate_blood_bm_annotation --datasets Greenleaf_pbmc_bm Shendure Bingren --seed_interval=1-10 --top_features_to_plot 1 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method=permutation_importance --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --woo_pcawg --mem_per_core 1 --top_features_to_plot_feat_imp 2 5 10 --create_bash_scripts --grid_analysis --grid_cell_types "bonemarrow B GL_BlBm,stomach Stromal cells SH,placenta PAEP_MECOM positive cells SH,thyroid Thyroid Follicular Cell BR"
## run robustness plotting as explained in README
#python3 prep_ML_model_scripts.py --cancer_types Bone-Leiomyo --scATAC_cell_number_filter 100 --annotation_dir new_intermediate_blood_bm_annotation --datasets Greenleaf_pbmc_bm Shendure Bingren --seed_interval=1-10 --top_features_to_plot 1 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method=permutation_importance --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --woo_pcawg --mem_per_core 1 --top_features_to_plot_feat_imp 2 5 10 --create_bash_scripts --grid_analysis --grid_cell_types "bonemarrow B GL_BlBm,stomach Stromal cells SH,placenta PAEP_MECOM positive cells SH,thyroid Thyroid Follicular Cell BR"
## run robustness plotting as explained in README
#python3 prep_ML_model_scripts.py --cancer_types Thy-AdenoCA --scATAC_cell_number_filter 100 --annotation_dir new_intermediate_blood_bm_annotation --datasets Greenleaf_pbmc_bm Shendure Bingren --seed_interval=1-10 --top_features_to_plot 1 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method=permutation_importance --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --woo_pcawg --mem_per_core 1 --top_features_to_plot_feat_imp 2 5 10 --create_bash_scripts --grid_analysis --grid_cell_types "bonemarrow B GL_BlBm,stomach Stromal cells SH,placenta PAEP_MECOM positive cells SH,thyroid Thyroid Follicular Cell BR"
## run robustness plotting as explained in README
#python3 prep_ML_model_scripts.py --cancer_types Uterus-AdenoCA --scATAC_cell_number_filter 100 --annotation_dir new_intermediate_blood_bm_annotation --datasets Greenleaf_pbmc_bm Shendure Bingren --seed_interval=1-10 --top_features_to_plot 1 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method=permutation_importance --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --woo_pcawg --mem_per_core 1 --top_features_to_plot_feat_imp 2 5 10 --create_bash_scripts --grid_analysis --grid_cell_types "bonemarrow B GL_BlBm,stomach Stromal cells SH,placenta PAEP_MECOM positive cells SH,thyroid Thyroid Follicular Cell BR"
## run robustness plotting as explained in README


### Fig 4D ###
#python3 ../analysis/ML/prep_ML_model_scripts.py --add_perf_to_file --cancer_types Biliary-AdenoCA --scATAC_cell_number_filter 100 --annotation_dir finalized_annotation --datasets Bingren Shendure Greenleaf_colon Greenleaf_pbmc_bm Tsankov Yang_kidney --seed_interval=1-10 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method=permutation_importance --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --woo_pcawg --mem_per_core 1 --top_features_to_plot_feat_imp 10 5 2 --submit_jobs 
## run robustness plotting as explained in README

#python3 ../analysis/ML/prep_ML_model_scripts.py --add_perf_to_file --cancer_types Bladder-TCC --scATAC_cell_number_filter 100 --annotation_dir finalized_annotation --datasets Bingren Shendure Greenleaf_colon Greenleaf_pbmc_bm Tsankov Yang_kidney --seed_interval=1-10 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method=permutation_importance --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --woo_pcawg --mem_per_core 1 --top_features_to_plot_feat_imp 10 5 2 --submit_jobs 
## run robustness plotting as explained in README

#python3 ../analysis/ML/prep_ML_model_scripts.py --add_perf_to_file --cancer_types Eso-AdenoCa --scATAC_cell_number_filter 100 --annotation_dir finalized_annotation --datasets Bingren Shendure Greenleaf_colon Greenleaf_pbmc_bm Tsankov Yang_kidney --seed_interval=1-10 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method=permutation_importance --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --woo_pcawg --mem_per_core 1 --top_features_to_plot_feat_imp 10 5 2 --submit_jobs 
## run robustness plotting as explained in README

#python3 ../analysis/ML/prep_ML_model_scripts.py --add_perf_to_file --cancer_types Stomach-AdenoCA --scATAC_cell_number_filter 100 --annotation_dir finalized_annotation --datasets Bingren Shendure Greenleaf_colon Greenleaf_pbmc_bm Tsankov Yang_kidney --seed_interval=1-10 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method=permutation_importance --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --woo_pcawg --mem_per_core 1 --top_features_to_plot_feat_imp 10 5 2 --submit_jobs 
## run robustness plotting as explained in README


### Fig 4E ###
TODO: Instructions for creating data

