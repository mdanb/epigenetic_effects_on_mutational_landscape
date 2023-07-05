#python3 build_ML_model.py --cancer_types sarcomatoid_waddell_mesomics --datasets Bingren Greenleaf_brain Shendure Yang_kidney --scATAC_cell_number_filter=30 --annotation_dir=finalized_annotation --meso --ML_model=XGB --seed=42:42 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50

#python3 build_ML_model.py --test_backward_selection_iter 17 --cancer_types sarcomatoid_waddell_mesomics --datasets Tsankov --scATAC_cell_number_filter=30 --annotation_dir=finalized_annotation --meso --ML_model=XGB --seed_range=1-100 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --iters_dont_skip=17

#python3 build_ML_model.py --cancer_types sarcomatoid_waddell_mesomics --datasets Tsankov --scATAC_cell_number_filter=30 --annotation_dir=finalized_annotation --meso --ML_model=XGB --seed_range=101-101 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --iters_dont_skip=17

#python3 build_ML_model.py --test_backward_selection_iter 15 --cancer_types sarcomatoid_waddell_mesomics --datasets Tsankov --scATAC_cell_number_filter=30 --annotation_dir=finalized_annotation --meso --ML_model=XGB --seed_range=1-100 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --iters_dont_skip=17

#python3 build_ML_model.py --cancer_types SCLC --datasets Greenleaf_brain Bingren Shendure Tsankov --scATAC_cell_number_filter=30 --annotation_dir=finalized_annotation --ML_model=XGB --seed=42:42 --SCLC --n_optuna_trials_prebackward_selection=0 --n_optuna_trials_backward_selection=0

#python3 build_ML_model.py --cancer_types Lung-AdenoCA Lung-SCC --datasets Tsankov --scATAC_cell_number_filter=30 --annotation_dir=finalized_annotation --ML_model=XGB --seed_range=1-1 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --iters_dont_skip=17

#python3 /broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML/build_ML_model.py --cancer_types Lung-AdenoCA --datasets Tsankov --scATAC_cell_number_filter 30 --annotation_dir finalized_annotation --seed_range 6-10 --iters_dont_skip 17 15 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50

#python3 /broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML/build_ML_model.py --cancer_types Lung-SCC --datasets Tsankov --scATAC_cell_number_filter 30 --annotation_dir finalized_annotation --seed_range 6-10 --iters_dont_skip 17 15 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50

#python3 /broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML/build_ML_model.py --cancer_types Lung-AdenoCA --datasets Tsankov --scATAC_cell_number_filter 30 --annotation_dir finalized_annotation --seed_range 6-10 --iters_dont_skip 17 15 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50

#python3 /broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML/build_ML_model.py --cancer_types Lung-AdenoCA --datasets Tsankov --scATAC_cell_number_filter 30 --annotation_dir finalized_annotation --seed_range 1-5 --iters_dont_skip 17 15 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50

python3 /broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML/build_ML_model.py --cancer_types Lung-SCC --datasets Tsankov --scATAC_cell_number_filter 30 --annotation_dir finalized_annotation --seed_range 1-5 --iters_dont_skip 17 15 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50
