#!/bin/bash
#$ -cwd
#$ -l h_rt=24:00:00
#$ -l os=RedHat7
#$ -pe smp 8
#$ -l h_vmem=1G
#$ -N multiple_myeloma_seed_range_6-10_fold_for_test_10
#$ -o /broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML
#$ -e /broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML
#$ -binding linear:8
python3 /broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML/build_ML_model.py --cancer_types multiple_myeloma --datasets Bingren Greenleaf_colon Greenleaf_pbmc_bm Shendure Tsankov Yang_kidney --scATAC_cell_number_filter 100 --annotation_dir new_intermediate_blood_bm_annotation --seed_range 6-10 --top_features_to_plot 1 2 5 10 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method permutation_importance --fold_for_test_set 10 --mm --tissues_to_consider all --test_set_perf_num_features all