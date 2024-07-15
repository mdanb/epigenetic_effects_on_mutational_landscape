### Fig 4A ###
#python3 ../analysis/ML/prep_ML_model_scripts.py  --add_p_to_file --seed_interval_step=5 --add_perf_to_file --cancer_types CNS-Medullo --scATAC_cell_number_filter 100 --annotation_dir finalized_annotation --datasets Bingren Bingren_adult_brain Shendure Greenleaf_brain --seed_interval=1-10 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --woo_pcawg --mem_per_core 1 --tissues_to_consider adult_brain brain frontal_cortex cerebrum cerebellum --top_features_to_plot_feat_imp 5 #-submit_jobs

## run robustness plotting as explained in README
#python3 ../analysis/ML/prep_ML_model_scripts.py --add_p_to_file --seed_interval_step=5 --cancer_types CNS-PiloAstro --scATAC_cell_number_filter 100 --annotation_dir finalized_annotation --datasets Bingren Bingren_adult_brain Shendure Greenleaf_brain --seed_interval=1-10 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --woo_pcawg --mem_per_core 1 --tissues_to_consider adult_brain brain frontal_cortex cerebrum cerebellum --top_features_to_plot_feat_imp 5 #--submit_jobs
## run robustness plotting as explained in README
#python3 ../analysis/ML/prep_ML_model_scripts.py --add_p_to_file --seed_interval_step=5 --cancer_types CNS-GBM --scATAC_cell_number_filter 100 --annotation_dir finalized_annotation --datasets Bingren Bingren_adult_brain Shendure Greenleaf_brain --seed_interval=1-10 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --woo_pcawg --mem_per_core 1 --tissues_to_consider adult_brain brain frontal_cortex cerebrum cerebellum --top_features_to_plot_feat_imp 5 #--submit_jobs

#python3 ../analysis/ML/prep_ML_model_scripts.py --add_p_to_file --seed_interval_step=5 --cancer_types CNS-GBM --scATAC_cell_number_filter 100 --annotation_dir finalized_annotation --datasets Bingren Bingren_adult_brain Shendure Greenleaf_brain --seed_interval=1-10 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --woo_pcawg --mem_per_core 1 --tissues_to_consider adult_brain frontal_cortex cerebrum brain cerebellum --top_features_to_plot_feat_imp 5
## run robustness plotting as explained in README

#python3 ../analysis/ML/prep_ML_model_scripts.py --add_p_to_file --seed_interval_step=5 --cancer_types CNS-Oligo --scATAC_cell_number_filter 100 --annotation_dir finalized_annotation --datasets Bingren Bingren_adult_brain Shendure Greenleaf_brain --seed_interval=1-10 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --woo_pcawg --mem_per_core 1 --tissues_to_consider adult_brain brain frontal_cortex cerebrum cerebellum --top_features_to_plot_feat_imp 5
## run robustness plotting as explained in README

### Fig 4B ###
#Rscript ../analysis/ArchR_analysis/reannotate_datasets.R --cores=4 --dataset=Greenleaf_brain --metadata_for_celltype_fn=GSE162170_atac_cell_metadata.txt.gz --sep_for_metadata=\t --cell_type_col_in_metadata=cell_type --tissue=all --nfrags_filter=1 --tss_filter=0 --min_cells_per_cell_type=1 --filter_per_cell_type --harmonize --fig4 --plot_custom_column --get_metacells

#Rscript ../data/scripts/create_count_overlaps.R --dataset=Greenleaf_brain --cores=4 --annotation=Greenleaf_brain_lowest_level_annotation --overlaps_per_cell

#Rscript ../analysis/correlation_analysis_clean.R --fig4_oligo
#Rscript ../analysis/correlation_analysis_clean.R --fig4_astro
#Rscript ../analysis/correlation_analysis_clean.R --fig4_gbm


### Fig 4C ###
#python3 ../analysis/ML/prep_ML_model_scripts.py --add_perf_to_file_grid --cancer_types Lymph-BNHL --scATAC_cell_number_filter 100 --annotation_dir new_intermediate_blood_bm_annotation --datasets Greenleaf_pbmc_bm Shendure Bingren --seed_interval=1-10 --top_features_to_plot 1 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method=permutation_importance --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --woo_pcawg --mem_per_core 1 --top_features_to_plot_feat_imp 2 5 10 --grid_analysis --grid_cell_types "bonemarrow B GL_BlBm,stomach Stromal cells SH,placenta PAEP_MECOM positive cells SH,thyroid Thyroid Follicular Cell BR"
## run robustness plotting as explained in README
#python3 ../analysis/ML/prep_ML_model_scripts.py --add_perf_to_file_grid --cancer_types Bone-Leiomyo --scATAC_cell_number_filter 100 --annotation_dir new_intermediate_blood_bm_annotation --datasets Greenleaf_pbmc_bm Shendure Bingren --seed_interval=1-10 --top_features_to_plot 1 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method=permutation_importance --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --woo_pcawg --mem_per_core 1 --top_features_to_plot_feat_imp 2 5 10 --grid_analysis --grid_cell_types "bonemarrow B GL_BlBm,stomach Stromal cells SH,placenta PAEP_MECOM positive cells SH,thyroid Thyroid Follicular Cell BR"
## run robustness plotting as explained in README
#python3 ../analysis/ML/prep_ML_model_scripts.py --add_perf_to_file_grid --cancer_types Thy-AdenoCA --scATAC_cell_number_filter 100 --annotation_dir new_intermediate_blood_bm_annotation --datasets Greenleaf_pbmc_bm Shendure Bingren --seed_interval=1-10 --top_features_to_plot 1 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method=permutation_importance --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --woo_pcawg --mem_per_core 1 --top_features_to_plot_feat_imp 2 5 10 --grid_analysis --grid_cell_types "bonemarrow B GL_BlBm,stomach Stromal cells SH,placenta PAEP_MECOM positive cells SH,thyroid Thyroid Follicular Cell BR"
## run robustness plotting as explained in README
#python3 ../analysis/ML/prep_ML_model_scripts.py --add_perf_to_file_grid --cancer_types Uterus-AdenoCA --scATAC_cell_number_filter 100 --annotation_dir new_intermediate_blood_bm_annotation --datasets Greenleaf_pbmc_bm Shendure Bingren --seed_interval=1-10 --top_features_to_plot 1 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method=permutation_importance --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --woo_pcawg --mem_per_core 1 --top_features_to_plot_feat_imp 2 5 10 --grid_analysis --grid_cell_types "bonemarrow B GL_BlBm,stomach Stromal cells SH,placenta PAEP_MECOM positive cells SH,thyroid Thyroid Follicular Cell BR"

#Rscript ../analysis/ML/plot_top_features.R --cancer_types="Lymph-BNHL,Bone-Leiomyo,Thy-AdenoCA,Uterus-AdenoCA" --datasets="Bingren,Shendure,Greenleaf_colon,Greenleaf_blood_bm,Tsankov" --ML_model=XGB --annotation="finalized_annotation,new_intermediate_blood_bm_annotation,finalized_annotation,finalized_annotation" --robustness_analysis --seed_range=1-10 --feature_importance_method=permutation_importance --folds_for_test_set=1-10 --grid_analysis --top_features_to_plot=1 --grid_cell_types="bonemarrow B GL_BlBm,stomach Stromal cells SH,thyroid Thyroid Follicular Cell BR,placenta PAEP_MECOM positive cells SH"
## run robustness plotting as explained in README


### Fig 4D ###
#python3 ../analysis/ML/prep_ML_model_scripts.py --add_p_to_file --cancer_types Biliary-AdenoCA --scATAC_cell_number_filter 100 --annotation_dir finalized_annotation --datasets Bingren Shendure Greenleaf_colon Greenleaf_pbmc_bm Tsankov Yang_kidney --seed_interval=1-10 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method=permutation_importance --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --woo_pcawg --mem_per_core 1 --top_features_to_plot_feat_imp 10 5 2 #--submit_jobs 
## run robustness plotting as explained in README

#python3 ../analysis/ML/prep_ML_model_scripts.py --add_p_to_file --cancer_types Bladder-TCC --scATAC_cell_number_filter 100 --annotation_dir finalized_annotation --datasets Bingren Shendure Greenleaf_colon Greenleaf_pbmc_bm Tsankov Yang_kidney --seed_interval=1-10 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method=permutation_importance --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --woo_pcawg --mem_per_core 1 --top_features_to_plot_feat_imp 10 5 2
## run robustness plotting as explained in README

#python3 ../analysis/ML/prep_ML_model_scripts.py --add_p_to_file --cancer_types Eso-AdenoCa --scATAC_cell_number_filter 100 --annotation_dir finalized_annotation --datasets Bingren Shendure Greenleaf_colon Greenleaf_pbmc_bm Tsankov Yang_kidney --seed_interval=1-10 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method=permutation_importance --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --woo_pcawg --mem_per_core 1 --top_features_to_plot_feat_imp 10 5 2 --top_features_to_plot 1
## run robustness plotting as explained in README

#python3 ../analysis/ML/prep_ML_model_scripts.py --add_p_to_file --cancer_types Stomach-AdenoCA --scATAC_cell_number_filter 100 --annotation_dir finalized_annotation --datasets Bingren Shendure Greenleaf_colon Greenleaf_pbmc_bm Tsankov Yang_kidney --seed_interval=1-10 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method=permutation_importance --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --woo_pcawg --mem_per_core 1 --top_features_to_plot_feat_imp 10 5 2
## run robustness plotting as explained in README


### Fig 4E ###
# Uncomment lines in script below per cancer type
# sh ../data/scripts/create_subsampled_data.sh 

#for i in 1 $(seq 5 5 35); do
#        python3 prep_ML_model_scripts.py --tissues_to_consider lung fetal_lung --cancer_types Lung-SCC_n_$i --scATAC_cell_number_filter 100 --annotation_dir finalized_annotation --datasets Bingren Shendure Rawlins_fetal_lung Tsankov --seed_interval=1-10 --top_features_to_plot 1 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method=permutation_importance --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --subsampled_mutations --mem_per_core 1 --test_set_perf_num_features all --top_features_to_plot_feat_imp 1 2 5 10 --cell_types_keep="lung Neuroendocrine-Tsankov" --submit_jobs
#done

#for i in 1 $(seq 5 5 35); do
#        python3 prep_ML_model_scripts.py --tissues_to_consider all --cancer_types Kidney-ChRCC_n_$i --scATAC_cell_number_filter 100 --annotation_dir finalized_annotation --datasets Bingren Greenleaf_pbmc_bm Greenleaf_colon Shendure Tsankov Yang_kidney --seed_interval=1-10 --top_features_to_plot 1 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method=permutation_importance --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --subsampled_mutations --mem_per_core 1 --test_set_perf_num_features all --top_features_to_plot_feat_imp 1 2 5 10 --submit_jobs
#done

#for i in 1 $(seq 5 5 35); do
#       python3 prep_ML_model_scripts.py --tissues_to_consider adult_brain frontal_cortex cerebrum brain cerebellum --cancer_types CNS-GBM_n_$i --scATAC_cell_number_filter 100 --annotation_dir finalized_annotation --datasets Bingren Bingren_adult_brain Greenleaf_brain Shendure --seed_interval=1-10 --top_features_to_plot 1 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method=permutation_importance --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --subsampled_mutations --mem_per_core 1 --test_set_perf_num_features all --top_features_to_plot_feat_imp 1 2 5 10 --submit_jobs
#done

#for i in 1 $(seq 5 5 35); do
#       python3 prep_ML_model_scripts.py --tissues_to_consider adult_brain frontal_cortex cerebrum brain cerebellum --cancer_types CNS-Medullo_n_$i --scATAC_cell_number_filter 100 --annotation_dir finalized_annotation --datasets Bingren Bingren_adult_brain Greenleaf_brain Shendure --seed_interval=1-10 --top_features_to_plot 1 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method=permutation_importance --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --subsampled_mutations --mem_per_core 1 --test_set_perf_num_features all --top_features_to_plot_feat_imp 1 2 5 10 --submit_jobs
#done

#for i in 1 $(seq 5 5 35); do
#        python3 prep_ML_model_scripts.py --tissues_to_consider all --cancer_types Skin-Melanoma_n_$i --scATAC_cell_number_filter 100 --annotation_dir finalized_annotation --datasets Bingren Greenleaf_pbmc_bm Greenleaf_colon Shendure Tsankov Yang_kidney --seed_interval=1-10 --top_features_to_plot 1 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method=permutation_importance --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --subsampled_mutations --mem_per_core 1 --test_set_perf_num_features all --top_features_to_plot_feat_imp 1 2 5 10 --submit_jobs
#done

#for i in 1 $(seq 5 5 35); do
#        python3 prep_ML_model_scripts.py --tissues_to_consider all --cancer_types Liver-HCC_n_$i --scATAC_cell_number_filter 100 --annotation_dir finalized_annotation --datasets Bingren Greenleaf_pbmc_bm Greenleaf_colon Shendure Tsankov Yang_kidney --seed_interval=1-10 --top_features_to_plot 1 --n_optuna_trials_prebackward_selection 50 --n_optuna_trials_backward_selection 50 --feature_importance_method=permutation_importance --fold_for_test_set_range 1-10 --test_set_perf_num_features all --cores=8 --seed_interval_step=5 --subsampled_mutations --mem_per_core 1 --test_set_perf_num_features all --top_features_to_plot_feat_imp 1 2 5 10 --submit_jobs
#done


