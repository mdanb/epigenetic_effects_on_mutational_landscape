#### Imports ####
import pandas as pd
import numpy as np
import pickle
import os
import time
from pipelinehelper import PipelineHelper
from sklearn.pipeline import Pipeline
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score
import argparse
import glob
from ML_utils import load_scATAC, load_agg_mutations, filter_agg_data, \
                     filter_mutations_by_cancer, load_meso_mutations

from itertools import chain

parser = argparse.ArgumentParser()

parser.add_argument('--cancer_types', nargs="+", type=str,
                    help='which cancer types to analyze', required=True)
parser.add_argument('--all_cells', action="store_true",
                    help='run model on all cells', default=False)
parser.add_argument('--tissue_spec_cells', action="store_true",
                     help='run model on tissue specific cells', default=False)
parser.add_argument('--clustered_mutations', action="store_true",
                    help='run model on hierarchically clustered mutations', default=False)
# parser.add_argument('--bing_ren', action="store_true",
#                     help='use Bing Ren ATACseq', default=False)
# parser.add_argument('--shendure', action="store_true",
#                     help='use Shendure ATACseq', default=False)
# parser.add_argument('--tsankov', action="store_true",
#                     help='use Tsankov ATACseq', default=False)
parser.add_argument('--datasets', nargs="+", type=str,
                    help='which sc-ATACseq datasets to analyze', required=True)
parser.add_argument('--scATAC_cell_number_filter', type=int,
                    help='minimum number of cells per cell type in scATAC', default=100)
# parser.add_argument('--combined_datasets', action="store_true",
#                     help='combine all scATACseq', default=False)
group = parser.add_mutually_exclusive_group()
# group.add_argument('--meso_waddell_and_biphasic', action="store_true",
#                     default=False)
# group.add_argument('--meso_waddell_only', action="store_true", default=False)
# group.add_argument('--meso_waddell_and_broad_only', action="store_true", default=False)
# group.add_argument('--meso_waddell_biph_786_846', action="store_true", default=False)

# group.add_argument('--waddell_biph', action="store_true",
#                     default=False)
# group.add_argument('--waddell_sarc', action="store_true", default=False)
# group.add_argument('--waddell_sarc_tsankov_sarc', action="store_true", default=False)
# group.add_argument('--waddell_sarc_biph_tsankov_sarc_biph_waddell_epith', action="store_true", default=False)

group.add_argument('--waddell_sarc_biph', action="store_true",
                    default=False)
group.add_argument('--waddell_sarc', action="store_true", default=False)
group.add_argument('--waddell_sarc_tsankov_sarc', action="store_true", default=False)
group.add_argument('--waddell_sarc_biph_tsankov_sarc_biph', action="store_true", default=False)


# parser.add_argument('--tss_filtered', action="store_true",
#                     help='Use TSS filtered data', default=False)
parser.add_argument('--tss_fragment_filter', type=int, default=-1)
# parser.add_argument('--bioRxiv_method', action="store_true",
#                     help='Use method from bioRxiv paper. Different from tissue_spec_cells by fact that here ' \
#                          'dont a priori know which tissue to train on, so trains on all of them', default=False)

config = parser.parse_args()
cancer_types = config.cancer_types
run_all_cells = config.all_cells
run_tissue_spec_cells = config.tissue_spec_cells
run_clustered_mutations = config.clustered_mutations
# bing_ren = config.bing_ren
# shendure = config.shendure
# tsankov = config.tsankov
datasets = sorted(config.datasets)
scATAC_cell_number_filter = config.scATAC_cell_number_filter
# meso_waddell_and_biphasic = config.meso_waddell_and_biphasic
# meso_waddell_only = config.meso_waddell_only
# meso_waddell_and_broad_only = config.meso_waddell_and_broad_only
# meso_waddell_biph_786_846 = config.meso_waddell_biph_786_846

waddell_sarc_biph = config.waddell_sarc_biph
waddell_sarc = config.waddell_sarc
waddell_sarc_tsankov_sarc = config.waddell_sarc_tsankov_sarc
waddell_sarc_biph_tsankov_sarc_biph = config.waddell_sarc_biph_tsankov_sarc_biph

# tss_filtered = config.tss_filtered
tss_fragment_filter = config.tss_fragment_filter
# bioRxiv_method = config.bioRxiv_method

#### Helpers ####
# Dataframe curation
def add_dataset_origin_to_cell_types(df, dataset):
    if (dataset == "Bingren"):
        df.columns = [c + " BR" for c in df.columns]
    elif (dataset == "Shendure"):
        df.columns = [c + " SH" for c in df.columns]
    elif (dataset == "Tsankov"):
        df.columns = [c + " TS" for c in df.columns]
    elif (dataset == "Greenleaf_brain"):
        df.columns = [c + " GL_Br" for c in df.columns]
    elif (dataset == "Greenleaf_pbmc_bm"):
        df.columns = [c + " GL_BlBm" for c in df.columns]
    elif (dataset == "Yang_kidney"):
        df.columns = [c + " Y_K" for c in df.columns]
    return(df)

# Filter Data helpers
def filter_clustered_data(scATAC_df, mutations_df):
    return pd.DataFrame(scATAC_df, index=mutations_df.index)
# Split Train/Test helpers
def get_train_test_split(X, y, test_size):
    X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                        test_size=test_size, random_state=42)
    return X_train, X_test, y_train, y_test

# Cross Validation helpers
def grid_search(X_train, y_train, pipe, params, num_k_folds):
    start_time = time.time()
    grid_search_object = GridSearchCV(pipe, params, scoring="r2",
                                      cv=KFold(num_k_folds), n_jobs=-1, verbose=100)
    grid_search_object.fit(X_train, y_train)
    print(f"--- {time.time() - start_time} seconds ---")
    return grid_search_object

# Feature Selection helpers
def get_top_n_features(clf, n, features):
    feature_importances = clf.feature_importances_
    feat_importance_idx = np.argsort(feature_importances)[::-1]
    top_n_feats = features[feat_importance_idx][:n]
    return top_n_feats

def print_and_save_top_features(top_features, filepath):
    n = len(top_features)
    print(f"Top {n} features")
    f = open(filepath, "w")
    for idx, feature in enumerate(top_features):
        print(f"{idx+1}. {feature}")
        f.write(f"{idx+1}. {feature}\n")

def backward_eliminate_features(X_train, y_train, starting_clf, starting_n,
                                params, num_k_folds, backwards_elim_dir):
    os.makedirs(backwards_elim_dir, exist_ok=True)
    top_n_feats = get_top_n_features(starting_clf, starting_n, X_train.columns.values)
    all_feature_rankings = get_top_n_features(starting_clf, len(X_train.columns.values), X_train.columns.values)
    print_and_save_top_features(all_feature_rankings, filepath=f"{backwards_elim_dir}/all_features_rankings.txt")
    X_train = X_train.loc[:, top_n_feats]
    print_and_save_top_features(top_n_feats,
                                filepath=f"{backwards_elim_dir}/starter_model_top_features.txt")
    for idx in range(1, len(top_n_feats)):
        filepath = f"{backwards_elim_dir}/top_features_iteration_{idx}.txt"
        if (not os.path.exists(filepath)):
            pipe = Pipeline([
                ('regressor', PipelineHelper([
                    ('rf', RandomForestRegressor(random_state=idx)),
                ])),
            ])
            grid_search_results = grid_search(X_train, y_train, pipe, params, num_k_folds)
            pickle.dump(grid_search_results, open(f"{backwards_elim_dir}/model_iteration_{idx}.pkl", 'wb'))
        grid_search_results = pickle.load(open(f"{backwards_elim_dir}/model_iteration_{idx}.pkl", 'rb'))
        best_model = grid_search_results.best_estimator_.get_params()['regressor__selected_model']
        top_n_feats = get_top_n_features(best_model, len(X_train.columns) - 1, X_train.columns)
        X_train = X_train.loc[:, top_n_feats]
        if (not os.path.exists(filepath)):
            print_and_save_top_features(top_n_feats, filepath=filepath)


# Test set performance helpers
def print_and_save_test_set_perf(X_test, y_test, model, filename):
    test_preds = model.predict(X_test)
    test_set_performance = r2_score(y_test, test_preds)
    with open(filename, "w") as f:
        f.write(str(test_set_performance))
    print(f"Test set performance: {test_set_performance}")

def train_val_test(scATAC_df, mutations, cv_filename, backwards_elim_dir,
                   test_set_perf_filename):
    X_train, X_test, y_train, y_test = get_train_test_split(scATAC_df, mutations,
                                                            0.10)
    pipe = Pipeline([
        ('regressor', PipelineHelper([
            ('rf', RandomForestRegressor(random_state=0)),
        ])),
    ])

    params = {
        'regressor__selected_model': pipe.named_steps['regressor'].generate({
            'rf__n_estimators':[10, 100, 1000],
        })
    }

    if (not os.path.exists(cv_filename)):
        grid_search_results = grid_search(X_train, y_train, pipe, params, 10)
        pickle.dump(grid_search_results, open(cv_filename, 'wb'))

    grid_search_results = pickle.load(open(cv_filename, 'rb'))
    print(f"Best Score: {grid_search_results.best_score_}\n")
    best_model = grid_search_results.best_estimator_.get_params()['regressor__selected_model']
    n = 20
    # top_n_feats = get_top_n_features(best_model, n, X_train.columns.values)
    # top_n_feats_tissue_spec = get_top_n_features(best_model_tissue_spec, n, X_train_tissue_spec.columns.values)
    #print_top_features(top_n_feats)

    backward_eliminate_features(X_train, y_train, best_model, n, params, 10, backwards_elim_dir)

    #### Test Set Performance ####
    print_and_save_test_set_perf(X_test, y_test, best_model, test_set_perf_filename)

# def run_unclustered_data_analysis(scATAC_df, run_all_cells, run_tissue_spec, cancer_types,
#                                   meso_waddell_and_biphasic, meso_waddell_only, meso_waddell_and_broad_only,
#                                   meso_waddell_biph_786_846, tss_fragment_filter, scATAC_source="Bingren"):


def run_unclustered_data_analysis(scATAC_df, run_all_cells, run_tissue_spec, cancer_types,
                                  waddell_sarc_biph, waddell_sarc,
                                  waddell_sarc_tsankov_sarc, waddell_sarc_biph_tsankov_sarc_biph,
                                  tss_fragment_filter, scATAC_source="Bingren"):
    # waddell_sarc_biph_waddell_epith = config.waddell_sarc_biph_waddell_epith
    # waddell_sarc_waddell_epith = config.waddell_sarc_waddell_epith
    # waddell_sarc_tsankov_sarc_waddell_epith = config.waddell_sarc_tsankov_sarc_waddell_epith
    # waddell_sarc_biph_tsankov_sarc_biph = config.waddell_sarc_biph_tsankov_sarc_biph
    if (waddell_sarc_biph or waddell_sarc or waddell_sarc_tsankov_sarc or
        waddell_sarc_biph_tsankov_sarc_biph):
        mutations_df = load_meso_mutations(waddell_sarc_biph, waddell_sarc,
                                           waddell_sarc_tsankov_sarc, waddell_sarc_biph_tsankov_sarc_biph)
    else:
        mutations_df = load_agg_mutations()

    #### Filter data ####
    scATAC_df, mutations_df = filter_agg_data(scATAC_df, mutations_df)

    scATAC_dir = f"scATAC_source_{scATAC_source}_cell_number_filter_{scATAC_cell_number_filter}"

    if (tss_fragment_filter != -1):
        scATAC_dir = scATAC_dir + "_tss_fragment_filter_" + str(tss_fragment_filter)

    if (waddell_sarc_biph):
        scATAC_dir = scATAC_dir + "_waddell_sarc_biph"
    elif (waddell_sarc):
        scATAC_dir = scATAC_dir + "_waddell_sarc"
    elif (waddell_sarc_tsankov_sarc):
        scATAC_dir = scATAC_dir + "_waddell_sarc_tsankov_sarc"
    elif (waddell_sarc_biph_tsankov_sarc_biph):
         scATAC_dir = scATAC_dir + "_waddell_sarc_biph_tsankov_sarc_biph"

    for cancer_type in cancer_types:
        os.makedirs(f"/broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML/models/{cancer_type}/{scATAC_dir}",
                    exist_ok=True)
        cancer_specific_mutations = filter_mutations_by_cancer(mutations_df, cancer_type)

        if (run_all_cells):
            backwards_elim_dir=f"/broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML/models/{cancer_type}/{scATAC_dir}/backwards_elimination_results"
            grid_search_filename = f"/broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML/models/{cancer_type}/{scATAC_dir}/grid_search_results.pkl"
            test_set_perf_filename = f"/broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML/models/{cancer_type}/{scATAC_dir}/test_set_performance.txt"

            # All Cells
            train_val_test(scATAC_df, cancer_specific_mutations,
                           grid_search_filename,
                           backwards_elim_dir,
                           test_set_perf_filename)

        # Tissue Specific
        if (run_tissue_spec):
            backwards_elim_dir=f"/broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML/models/{cancer_type}/{scATAC_dir}/backwards_elimination_results_tissue_spec"
            grid_search_filename = f"/broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML/models/{cancer_type}/{scATAC_dir}/grid_search_results_tissue_specific.pkl"
            test_set_perf_filename = f"/broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML/models/{cancer_type}/{scATAC_dir}/test_set_performance_tissue_spec.txt"

            tissue = cancer_type.split("-")[0]
            tissue_specific_cell_types = [cell_type for cell_type in scATAC_df.columns.values if tissue in cell_type]
            per_tissue_df = scATAC_df.loc[:, tissue_specific_cell_types]
            train_val_test(per_tissue_df,
                           cancer_specific_mutations,
                           grid_search_filename,
                           backwards_elim_dir,
                           test_set_perf_filename)
        # if (bioRxiv_method):
        #     tissues = set(scATAC_df.columns.str.split().to_series().apply(lambda x: x[0]))
        #     for tissue in tissues:
        #         grid_search_filename = f"models/bioRxiv_method/{cancer_type}/{scATAC_dir}/grid_search_results.pkl"
        #         test_set_perf_filename = f"models/bioRxiv_method/{cancer_type}/{scATAC_dir}/test_set_performance.txt"
        #
        #         train_val_test(scATAC_df, cancer_specific_mutations,
        #                        grid_search_filename,
        #                        backwards_elim_dir,
        #                        test_set_perf_filename)


def run_clustered_data_analysis(scATAC_df, cancer_types):
    for cancer_type in cancer_types:
        os.makedirs(f"/broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML/models/{cancer_type}", exist_ok=True)
        cancer_hierarchical_dir = f"../../data/processed_data/hierarchically_clustered_mutations/{cancer_type}"
        for cluster_method_dir in os.listdir(cancer_hierarchical_dir):
            for threshold_dir in os.listdir(f"{cancer_hierarchical_dir}/{cluster_method_dir}"):
                run_per_cluster_models(scATAC_df, cancer_type, cancer_hierarchical_dir, cluster_method_dir,
                                       threshold_dir)

def run_per_cluster_models(scATAC_df, cancer_type, cancer_hierarchical_dir, cluster_method_dir,
                           threshold_dir):
    for cluster_file in glob.glob(f"{cancer_hierarchical_dir}/{cluster_method_dir}/{threshold_dir}/aggregated*"):
        mutations_df = pd.read_csv(cluster_file, index_col=0)
        cluster_dir = cluster_file.split("/")[-1].split(".")[0]
        cluster_dir = f"/broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML/models/{cancer_type}/{cluster_method_dir}/{threshold_dir}/{cluster_dir}/"
        os.makedirs(cluster_dir, exist_ok=True)
        backwards_elim_dir=f"{cluster_dir}/backwards_elimination_results"
        grid_search_filename = f"{cluster_dir}/grid_search_results.pkl"
        test_set_perf_filename = f"{cluster_dir}/test_set_performance.txt"

        scATAC_df = filter_clustered_data(scATAC_df, mutations_df)
        mutations = mutations_df.values.reshape(-1)
        train_val_test(scATAC_df, mutations,
                       grid_search_filename,
                       backwards_elim_dir,
                       test_set_perf_filename)

#### Load scATAC ####
tss_filtered_root = "../../data/processed_data/count_overlap_data/tsse_filtered"
datasets_combined_count_overlaps = []
if (tss_fragment_filter != -1):
    # TODO: Fix for other than Bing Ren
    #result = pyreadr.read_r("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/mutation_data/hg19.1Mb.ranges.Polak.Nature2015.RData") #
    chr_ranges = pd.read_csv("../../data/processed_data/chr_ranges.csv")
    if (bing_ren):
        scATAC_df_bingren = load_scATAC(f"{tss_filtered_root}/bing_ren/combined/" \
                               f"combined_{tss_fragment_filter}_fragments.rds").T
        scATAC_df_bingren.index = chr_ranges["x"].values

    if (shendure):
        scATAC_df_shendure = load_scATAC(f"{tss_filtered_root}/shendure/combined/" \
                                    f"combined_{tss_fragment_filter}_fragments.rds").T
        scATAC_df_shendure.index = chr_ranges["x"].values
    if (tsankov):
        scATAC_df_tsankov = load_scATAC(f"{tss_filtered_root}/tsankov/combined/" \
                                     f"combined_{tss_fragment_filter}_fragments.rds").T
        scATAC_df_tsankov.index = chr_ranges["x"].values
else:
    for dataset in datasets:
        datasets_combined_count_overlaps.append(load_scATAC("../../data/processed_data/count_overlap_data/combined_count_overlaps" \
                                f"/{dataset}_count_filter_{scATAC_cell_number_filter}_combined_count_overlaps.rds"))
        # scATAC_df_bingren = load_scATAC("processed_data/count_overlap_data/combined_count_overlaps" \
        #                         f"/count_filter_{scATAC_cell_number_filter}_combined_count_overlaps.rds")
        # scATAC_df_shendure = load_scATAC("processed_data/count_overlap_data/combined_count_overlaps" \
        #                              f"/shendure_count_filter_{scATAC_cell_number_filter}_combined_count_overlaps.rds")
        # scATAC_df_tsankov = load_scATAC("processed_data/count_overlap_data/combined_count_overlaps" \
        #                              f"/tsankov_count_filter_{scATAC_cell_number_filter}_combined_count_overlaps.rds")

scATAC_df = pd.DataFrame()
scATAC_sources = ""

for idx, dataset in enumerate(datasets):
    if (scATAC_sources == ""):
        scATAC_sources = dataset
    else:
        scATAC_sources = "_".join((scATAC_sources, dataset))
    datasets_combined_count_overlaps[idx] = [add_dataset_origin_to_cell_types(datasets_combined_count_overlaps[idx],
                                             dataset)]

scATAC_df = pd.concat(chain(*datasets_combined_count_overlaps), axis=1)

# scATAC_df = pd.concat((scATAC_df, scATAC_df_bingren), axis=1)
# scATAC_sources = scATAC_sources + "bing_ren"
# scATAC_df = pd.concat((scATAC_df, scATAC_df_shendure), axis=1)
# scATAC_sources = scATAC_sources + "shendure"
# scATAC_df = pd.concat((scATAC_df, scATAC_df_tsankov), axis=1)
# scATAC_sources = scATAC_sources + "tsankov"

# if (len(datasets) == 5):
#     scATAC_sources = "combined_datasets"

if (run_all_cells or run_tissue_spec_cells):
    run_unclustered_data_analysis(scATAC_df, run_all_cells, run_tissue_spec_cells,
                                  cancer_types, waddell_sarc_biph, waddell_sarc,
                                  waddell_sarc_tsankov_sarc, waddell_sarc_biph_tsankov_sarc_biph,
                                  tss_fragment_filter, scATAC_sources)


# if (run_all_cells and shendure):
#     run_unclustered_data_analysis(scATAC_df_shendure, run_all_cells, run_tissue_spec_cells, cancer_types, meso,
#                                   "shendure")


# if (run_clustered_mutations and bing_ren):
#     run_clustered_data_analysis(scATAC_df, cancer_types)

# if (bioRxiv_method):
#     run_unclustered_data_analysis(combined_scATAC_df, run_all_cells, run_tissue_spec_cells,
#                                   cancer_types, meso, "combined_datasets")