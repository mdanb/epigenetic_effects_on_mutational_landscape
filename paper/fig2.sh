
### Fig 2A ###
cd ../analysis/ArchR_analysis/

Rscript reannotate_datasets.R --cores=8 --dataset=Greenleaf_colon --metadata_for_celltype_fn=greenleaf_colon_metadata.csv --sep_for_metadata=, --cell_type_col_in_metadata=general_cell_type --tissue=all --nfrags_filter=10000 --tss_filter=0 --cell_types=epithelial --min_cells_per_cell_type=1 --filter_per_cell_type --plot_custom_column --plus_to_add_to_metadata=GrossPathology,CellType --plus_filters="Normal|Unaffected" --color_embedding_by=CellType --harmonize --fig2_colon #--get_metacells 

### Fig 2B ###
#python3 ../analysis/ML/prep_ML_model_scripts.py --cancer_types mss --datasets Bingren Shendure Greenleaf_colon --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=8 --feature_importance_method=permutation_importance --test_set_perf_num_features all --top_features_to_plot_feat_imp 5 --tissues_to_consider colon_sigmoid colon_transverse intestine normal_colon --msi_high --submit_jobs
## run robustness plotting as explained in README

### Fig 2C ###
#cd ../analysis/ArchR_analysis

#Rscript reannotate_datasets.R --cores=8 --dataset=Greenleaf_pbmc_bm --metadata_for_celltype_fn=intermediate_blood_bm_annotation_metadata.csv --sep_for_metadata=, --cell_type_col_in_metadata=cell_type --tissue=all --nfrags_filter=1 --tss_filter=0 --min_cells_per_cell_type=100 --filter_per_cell_type --plot_custom_column  --color_embedding_by=cell_type --get_metacells
### Fig 2D ###
#python3 ../analysis/ML/prep_ML_model_scripts.py --cancer_types Lymph-CLL --datasets Greenleaf_pbmc_bm --scATAC_cell_number_filter=100 --annotation_dir=new_intermediate_blood_bm_annotation --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=8 --feature_importance_method=permutation_importance --test_set_perf_num_features all --top_features_to_plot_feat_imp 5 --tissues_to_consider all --woo_pcawg --submit_jobs
## run robustness plotting as explained in README

#python3 ../analysis/ML/prep_ML_model_scripts.py --cancer_types Myeloid-AML --datasets Greenleaf_pbmc_bm --scATAC_cell_number_filter=100 --annotation_dir=new_intermediate_blood_bm_annotation --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=8 --feature_importance_method=permutation_importance --test_set_perf_num_features all --top_features_to_plot_feat_imp 5 --tissues_to_consider all --woo_pcawg
## run robustness plotting as explained in README

