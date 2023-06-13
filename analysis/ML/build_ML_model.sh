#python3 build_ML_model.py --cancer_types sarcomatoid_waddell_mesomics --datasets Bingren Greenleaf_brain Shendure Yang_kidney --scATAC_cell_number_filter=30 --annotation_dir=finalized_annotation --meso --ML_model=XGB --seed=42:42 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50

python3 build_ML_model.py --cancer_types sarcomatoid_waddell_mesomics --datasets Tsankov --scATAC_cell_number_filter=30 --annotation_dir=finalized_annotation --meso --ML_model=XGB --seed=42:42 --n_optuna_trials_prebackward_selection=250 --n_optuna_trials_backward_selection=250

#python3 build_ML_model.py --cancer_types SCLC --datasets Greenleaf_brain Bingren Shendure Tsankov --scATAC_cell_number_filter=30 --annotation_dir=finalized_annotation --ML_model=XGB --seed=42:42 --SCLC --n_optuna_trials_prebackward_selection=0 --n_optuna_trials_backward_selection=0


