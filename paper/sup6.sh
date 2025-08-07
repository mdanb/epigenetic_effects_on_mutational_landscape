### Supplementary Figure 6 ###
## Supp Fig 6A

#python3 ../analysis/ML/prep_ML_model_scripts.py --cancer_types Biliary-AdenoCA --datasets Bingren Shendure Greenleaf_colon Greenleaf_pbmc_bm Tsankov Yang_kidney --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=8 --feature_importance_method=permutation_importance --test_set_perf_num_features all --top_features_to_plot_feat_imp 10 5 2 --tissues_to_consider all --woo_pcawg --add_p_to_file --add_perf_to_file --submit_jobs
## run robustness plotting as explained in README

#python3 ../analysis/ML/prep_ML_model_scripts.py --cancer_types Eso-AdenoCa --datasets Bingren Shendure Greenleaf_colon Greenleaf_pbmc_bm Tsankov Yang_kidney --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=8 --feature_importance_method=permutation_importance --test_set_perf_num_features all --top_features_to_plot_feat_imp 10 5 2 --tissues_to_consider all --woo_pcawg --add_p_to_file --add_perf_to_file --submit_jobs
## run robustness plotting as explained in README

#python3 ../analysis/ML/prep_ML_model_scripts.py --cancer_types Stomach-AdenoCA --datasets Bingren Shendure Greenleaf_colon Greenleaf_pbmc_bm Tsankov Yang_kidney --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 --cores=8 --feature_importance_method=permutation_importance --test_set_perf_num_features all --top_features_to_plot_feat_imp 10 5 2 --tissues_to_consider all --woo_pcawg --add_p_to_file --add_perf_to_file --submit_jobs
## run robustness plotting as explained in README


#python3 ../analysis/ML/prep_ML_model_scripts.py --cancer_types Metaplasia --datasets Bingren Greenleaf_pbmc_bm Greenleaf_colon Shendure Tsankov Yang_kidney --scATAC_cell_number_filter=100 --annotation_dir=finalized_annotation --dataset_abbrev D2 \
    # --seed_interval=1-10 --fold_for_test_set_range=1-10 --seed_interval_step=5 --n_optuna_trials_prebackward_selection=50 --n_optuna_trials_backward_selection=50 \
    # --cores=10 --feature_importance_method=permutation_importance --test_set_perf_num_features all --top_features_to_plot_feat_imp 5 \
    # --tissues_to_consider=all --custom_mutations --add_p_to_file --add_perf_to_file --mem_per_core=1000 --submit_jobs




## Supp Fig 6B
# cd ../analysis/ArchR_analysis
# Rscript reannotate_datasets.R \
#--cores=8
#--dataset=Shendure
#--metadata_for_celltype_fn=GSE149683_File_S2.Metadata_of_high_quality_cells.txt
#--sep_for_metadata=\t
#--cell_type_col_in_metadata=cell_type
#--tissue=all
#--nfrags_filter=1
#--tss_filter=0
#--min_cells_per_cell_type=1
#--filter_per_cell_type
#--marker_genes=ACTA2,TAGLN,MYH11,COL1A2,COL3A1,DCN





