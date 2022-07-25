import argparse
import pickle
from utils import load_scATAC, load_agg_mutations, filter_agg_data, filter_mutations_by_cancer
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure

parser = argparse.ArgumentParser()

parser.add_argument('--cancer_types', nargs="+", type=str,
                    help='which cancer types to analyze', required=True)
parser.add_argument('--all_cells', action="store_true",
                    help='perform error analysis on model trained with '
                         'all cells', default=False)
parser.add_argument('--combined_datasets', action="store_true",
                    help='perform error analysis on models trained with all scATACseq',
                    default=False)

config = parser.parse_args()
cancer_types = config.cancer_types
run_all_cells = config.all_cells
combined_datasets = config.combined_datasets
scATAC_sources = []

if (combined_datasets):
    scATAC_sources.append("combined_datasets")

scATAC_df = load_scATAC("processed_data/count_overlap_data/combined_count_overlaps" \
                        "/count_filter_100_combined_count_overlaps.rds")
scATAC_df_shendure = load_scATAC("processed_data/count_overlap_data/combined_count_overlaps" \
                                 "/shendure_count_filter_100_combined_count_overlaps.rds")
scATAC_df_tsankov = load_scATAC("processed_data/count_overlap_data/combined_count_overlaps" \
                                 "/tsankov_count_filter_100_combined_count_overlaps.rds")
combined_scATAC_df = pd.concat((scATAC_df, scATAC_df_shendure, scATAC_df_tsankov), axis=1)
mutations_df = load_agg_mutations()
combined_scATAC_df, mutations_df = filter_agg_data(combined_scATAC_df, mutations_df)

for scATAC_source in scATAC_sources:
    for cancer_type in cancer_types:
        mutations = filter_mutations_by_cancer(mutations_df, cancer_type)
        path = f"models/{cancer_type}/scATAC_source_{scATAC_source}/grid_search_results.pkl"
        gs_results = pickle.load(open(path, "rb"))
        best_model = gs_results.best_estimator_.get_params()['regressor__selected_model']
        errors = []
        for train_index, val_index in gs_results.cv.split(combined_scATAC_df, mutations):
            val_actual = mutations.iloc[val_index]
            val_set_preds = best_model.predict(combined_scATAC_df.iloc[val_index, :])
            error = np.absolute(val_actual - val_set_preds) / val_actual
            errors.append(error)
        errors = [err for err_list in errors for err in err_list]
        figure(figsize=(16, 12))
        x = np.arange(len(errors))
        plt.bar(x, np.asarray(errors, dtype=np.float32), width=4)
        for i in range(0, 2200, 213):
            plt.axvline(x=i, c="red")
        plt.xticks(np.arange(min(x), max(x)+1, 100))
        plt.xlabel("Bin")
        plt.ylabel("Percent error")
        # plt.title("Percent error versus chromosomal bins")
        plt.savefig('figures/bin_vs_percent_error.png')