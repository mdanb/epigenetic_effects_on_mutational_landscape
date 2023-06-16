#python3 build_ML_model.py --cancer_types sarcomatoid_waddell_mesomics --datasets Bingren Greenleaf_brain Shendure Yang_kidney --scATAC_cell_number_filter=30 --annotation_dir=finalized_annotation --meso --ML_model=XGB --seed=42:42 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50

python3 build_ML_model.py --test_backward_selection_iter 17 --cancer_types sarcomatoid_waddell_mesomics --datasets Tsankov --scATAC_cell_number_filter=30 --annotation_dir=finalized_annotation --meso --ML_model=XGB --seed_range=1-100 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --iters_dont_skip=17

#python3 build_ML_model.py --cancer_types SCLC --datasets Greenleaf_brain Bingren Shendure Tsankov --scATAC_cell_number_filter=30 --annotation_dir=finalized_annotation --ML_model=XGB --seed=42:42 --SCLC --n_optuna_trials_prebackward_selection=0 --n_optuna_trials_backward_selection=0




