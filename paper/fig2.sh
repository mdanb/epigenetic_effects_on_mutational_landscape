
### Fig 2A ###
cd ../data/scripts/
Rscript create_arrow_files_and_tss.R --dataset=Greenleaf_colon --cores=1
Rscript create_arrow_files_and_tss.R --dataset=Shendure --cores=1
Rscript create_arrow_files_and_tss.R --dataset=Greenleaf_brain --cores=1
Rscript create_ArchR_project.R --cores=8
cd ../../paper

### Fig 2B ###
#python3 ../analysis/ML/prep_ML_model_scripts.py --cancer_types mss --datasets Bingren Shendure Greenleaf_colon --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=8 --feature_importance_method=permutation_importance --test_set_perf_num_features all --top_features_to_plot_feat_imp 5 --tissues_to_consider colon_sigmoid,colon_transverse,intestine,normal_colon --msi_high --submit_jobs
## run robustness plotting as explained in README

### Fig 2C ###

### Fig 2D ###
#python3 ../analysis/ML/prep_ML_model_scripts.py --cancer_types Lymph-CLL --datasets Greenleaf_pbmc_bm --scATAC_cell_number_filter=100 --annotation_dir=new_intermediate_blood_bm_annotation --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=8 --feature_importance_method=permutation_importance --test_set_perf_num_features all --top_features_to_plot_feat_imp 5 --tissues_to_consider all --woo_pcawg --submit_jobs

#python3 ../analysis/ML/prep_ML_model_scripts.py --cancer_types Myeloid-AML --datasets Greenleaf_pbmc_bm --scATAC_cell_number_filter=100 --annotation_dir=new_intermediate_blood_bm_annotation --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=8 --feature_importance_method=permutation_importance --test_set_perf_num_features all --top_features_to_plot_feat_imp 5 --tissues_to_consider all --woo_pcawg --submit_jobs

## run robustness plotting as explained in README

