#python3 prep_ML_model_scripts.py --cancer_types mesomics_top_r1_90th_perc --datasets Tsankov --scATAC_cell_number_filter=30 --annotation_dir=finalized_annotation --seed_interval=1-100 --seed_interval_step=5 --iters_dont_skip 17 15 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --meso --cores=8 
#python3 prep_ML_model_scripts.py --cancer_types mesomics_top_r2_90th_perc --datasets Tsankov --scATAC_cell_number_filter=30 --annotation_dir=finalized_annotation --seed_interval=1-100 --seed_interval_step=5 --iters_dont_skip 17 15 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --meso --cores=8
#python3 prep_ML_model_scripts.py --cancer_types mesomics_top_r3_90th_perc --datasets Tsankov --scATAC_cell_number_filter=30 --annotation_dir=finalized_annotation --seed_interval=1-100 --seed_interval_step=5 --iters_dont_skip 17 15 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --meso --cores=8

#python3 prep_ML_model_scripts.py --cancer_types blum_top_10_perc --datasets Tsankov --meso --scATAC_cell_number_filter=30 --annotation_dir=finalized_annotation --seed_interval=1-100 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=2

#python3 prep_ML_model_scripts.py --cancer_types nmf_top_10_perc --datasets Tsankov --meso --scATAC_cell_number_filter=30 --annotation_dir=finalized_annotation --seed_interval=1-100 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=2

python3 prep_ML_model_scripts.py --make_plots --cancer_types Lung-SCC --datasets Tsankov --scATAC_cell_number_filter=30 --annotation_dir=finalized_annotation --seed_interval=1-100 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=2 --top_features_to_plot "19" "15" "10" "5" "2" --feature_importance_method=permutation_importance

#python3 prep_ML_model_scripts.py --cancer_types blum_top_10_perc nmf_top_10_perc blum_50_top_10_perc --datasets Tsankov --meso --scATAC_cell_number_filter=30 --annotation_dir=finalized_annotation --seed_interval=1-100 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=4

#python3 prep_ML_model_scripts.py --cancer_types blum_bottom_10_perc nmf_bottom_10_perc blum_50_bottom_10_perc --datasets Tsankov --meso --scATAC_cell_number_filter=30 --annotation_dir=finalized_annotation --seed_interval=1-100 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=4


#python3 prep_ML_model_scripts.py --cancer_types sarcomatoid_waddell_mesomics --datasets Tsankov --scATAC_cell_number_filter=30 --annotation_dir=finalized_annotation --seed_interval=1-100 --seed_interval_step=5 --iters_dont_skip 17 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --meso 
