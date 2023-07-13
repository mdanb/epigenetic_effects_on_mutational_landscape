Rscript plot_top_features.R --robustness_analysis --datasets="Tsankov" --cancer_types=Lung-SCC --cell_number_filter=30 --ML_model=XGB --annotation=finalized_annotation --top_features_to_plot=2,5,10,15,19 --robustness_test_perf_boxplot

Rscript plot_top_features.R --robustness_analysis --datasets="Tsankov" --cancer_types=Lung-AdenoCA --cell_number_filter=30 --ML_model=XGB --annotation=finalized_annotation --top_features_to_plot=2,5,10,15,19 --robustness_test_perf_boxplot
