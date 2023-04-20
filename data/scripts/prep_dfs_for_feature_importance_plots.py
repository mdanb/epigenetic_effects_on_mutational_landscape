import pickle
import argparse
import glob
from natsort import natsorted
import pandas as pd
import os
from pathlib import Path

parser = argparse.ArgumentParser()

parser.add_argument('--cancer_types', nargs="+", type=str,
                    help='which cancer types to analyze', required=True,
                    )
parser.add_argument('--clustered_mutations', action="store_true",
                    help='obtain model run on hierarchically clustered mutations',
                    default=False)
parser.add_argument('--datasets', nargs="+", type=str,
                    help='which sc-ATACseq datasets to analyze', required=True)
parser.add_argument('--annotation', type=str, default="default_annotation")
parser.add_argument('--tissues_to_consider', nargs="+", type=str, default="all")
parser.add_argument('--iters_dont_skip', nargs="+", type=int)
parser.add_argument('--ML_model', type=str)

# parser.add_argument('--bing_ren', action="store_true",
#                     help='obtain model which used Bing Ren ATACseq', default=False)
# parser.add_argument('--shendure', action="store_true",
#                     help='obtain model which used Shendure ATACseq', default=False)
# parser.add_argument('--tsankov', action="store_true",
#                     help='obtain model which used tsankov ATACseq', default=False)
# parser.add_argument('--combined_datasets', action="store_true",
#                     help='obtain model which combined all scATACseq', default=False)
parser.add_argument('--cell_number_filter', type=int)
parser.add_argument('--num_iter_skips', type=int, default=5)
parser.add_argument('--tss_fragment_filter', nargs="+", type=int, default=[-1])
group = parser.add_mutually_exclusive_group()
group.add_argument('--waddell_sarc_biph', action="store_true",
                    default=False)
group.add_argument('--waddell_sarc', action="store_true", default=False)
group.add_argument('--waddell_sarc_tsankov_sarc', action="store_true", default=False)
group.add_argument('--waddell_sarc_biph_tsankov_sarc_biph', action="store_true", default=False)


def construct_backwards_elim_dir(cancer_type, scATAC_source, cell_number_filter,
                                 tss_fragment_filter, waddell_sarc_biph,
                                 waddell_sarc, waddell_sarc_tsankov_sarc,
                                 waddell_sarc_biph_tsankov_sarc_biph,
                                 annotation, tissues_to_consider, ML_model):
    dir = f"/broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML/models/{ML_model}/{cancer_type}/scATAC_source_" \
          f"{scATAC_source}_cell_number_filter_{cell_number_filter}"
    if (tss_fragment_filter != -1):
        dir = dir + f"_tss_fragment_filter_{tss_fragment_filter}"

    if (waddell_sarc_biph):
        dir = dir + "_waddell_sarc_biph"
    elif (waddell_sarc):
        dir = dir + "_waddell_sarc"
    elif (waddell_sarc_tsankov_sarc):
        dir = dir + "_waddell_sarc_tsankov_sarc"
    elif (waddell_sarc_biph_tsankov_sarc_biph):
        dir = dir + "_waddell_sarc_biph_tsankov_sarc_biph"

    dir = dir + f"_annotation_{annotation}"

    if (tissues_to_consider == "all"):
        return f"{dir}/backwards_elimination_results/"
    else:
        return f"{dir}/backwards_elimination_results_{tissues_to_consider}"

def get_relevant_backwards_elim_dirs(config):
    print("cancer types")
    cancer_types = config.cancer_types
    print(cancer_types)

    # tissue_spec_cells = config.tissue_spec_cells
    # clustered_mutations = config.clustered_mutations
    # bing_ren = config.bing_ren
    # shendure = config.shendure
    # tsankov = config.tsankov
    datasets = config.datasets
    waddell_sarc_biph = config.waddell_sarc_biph
    waddell_sarc = config.waddell_sarc
    waddell_sarc_tsankov_sarc = config.waddell_sarc_tsankov_sarc
    waddell_sarc_biph_tsankov_sarc_biph = config.waddell_sarc_biph_tsankov_sarc_biph
    # combined_datasets = config.combined_datasets
    cell_number_filter = config.cell_number_filter
    tss_fragment_filter = config.tss_fragment_filter
    tissues_to_consider = "_".join(config.tissues_to_consider)
    annotation = config.annotation
    ML_model = config.ML_model
    backward_elim_dirs = []

    scATAC_sources = ""
    for idx, dataset in enumerate(datasets):
        if (scATAC_sources == ""):
            scATAC_sources = dataset
        else:
            scATAC_sources = "_".join((scATAC_sources, dataset))

    # scATAC_sources = ""
    # if (bing_ren):
    #     scATAC_sources = scATAC_sources + "bing_ren"
    # if (shendure):
    #     scATAC_sources = scATAC_sources + "shendure"
    # if (tsankov):
    #     scATAC_sources = scATAC_sources + "tsankov"

    # if (bing_ren and shendure and tsankov):
    #     scATAC_sources = "combined_datasets"

    for cancer_type in cancer_types:
        for tss_filter in tss_fragment_filter:
            backward_elim_dirs.append(construct_backwards_elim_dir(cancer_type, scATAC_sources, cell_number_filter,
                                                                   tss_filter, waddell_sarc_biph,
                                                                   waddell_sarc,
                                                                   waddell_sarc_tsankov_sarc,
                                                                   waddell_sarc_biph_tsankov_sarc_biph,
                                                                   annotation,
                                                                   tissues_to_consider,
                                                                   ML_model))
    print(backward_elim_dirs)
    return backward_elim_dirs
        # if (run_tissue_spec):
        #     backwards_elim_dir=f"models/{cancer_type}/scATAC_source_{scATAC_source}/" \
        #                        f"backwards_elimination_results_tissue_spec"
        #     grid_search_filename = f"models/{cancer_type}/scATAC_source_{scATAC_source}/" \
        #                            f"grid_search_results_tissue_specific.pkl"
        #     test_set_perf_filename = f"models/{cancer_type}/scATAC_source_{scATAC_source}/" \
        #                              f"test_set_performance_tissue_spec.txt"
        # if (clustered_mutations):
        # cancer_hierarchical_dir = f"processed_data/hierarchically_clustered_mutations/{cancer_type}"
        #     for cluster_method_dir in os.listdir(cancer_hierarchical_dir):
        #         for threshold_dir in os.listdir(f"{cancer_hierarchical_dir}/{cluster_method_dir}"):
        #             run_per_cluster_models(scATAC_df, cancer_type, cancer_hierarchical_dir, cluster_method_dir,
        #                                    threshold_dir)

def prep_df_for_feat_importance_plots(backwards_elim_dirs, num_iter_skips, iters_dont_skip, ML_model):
    for backwards_elim_dir in backwards_elim_dirs:
        df = pd.DataFrame(columns = ["features", "importance", "num_features", "score"])

        temp = backwards_elim_dir.split("/")
        cancer_type = temp[-4]
        cancer_type_dir = temp[-3]

        figure_path = os.path.join("../../figures/models",
                                   ML_model,
                                   cancer_type,
                                   cancer_type_dir,
                                   "backwards_elimination_results")
        #os.makedirs(figure_path, exist_ok=True)
        files = natsorted(glob.glob(f"{backwards_elim_dir}/*pkl"))
        for idx, file in enumerate(files):
            if (idx % num_iter_skips == 0 or idx in iters_dont_skip):
                gs = pickle.load(open(file, "rb"))
                model = gs.best_estimator_.get_params()['regressor__selected_model']
                cv_score = gs.best_score_
                features = model.feature_names_in_
                feature_importances = model.feature_importances_
                df_curr = pd.DataFrame((features, feature_importances, [len(features)] * len(features),
                                        [cv_score] * len(features))).T
                df_curr.columns = df.columns
                df = pd.concat((df, df_curr))
        Path(figure_path).mkdir(parents=True, exist_ok=True)
        df.to_csv(os.path.join(figure_path, "df_for_feature_importance_plots.csv"), index=False)
        # import matplotlib.pyplot as plt
        # ax = df.plot.bar(x = "num_features", y = "importance", stacked=True)
        # fig = ax.get_figure()
        # fig.savefig(os.path.join(figure_path, "stacked_feature_importances.png"))
        # plot = (
        #         ggplot(df, aes(y="importance", fill="features")) +
        #             geom_bar(stat="identity", width=1) +
        #             facet_wrap("~num_features")
        #             #coord_polar("y", start=0)
        #         )
        # plot = (
        #         ggplot(df, aes(x="num_features", fill="features"))
        #               + geom_bar(fill="position")
        #               + xlab("Number of Features")
        #               + ylab("")
        #         )
        # plot = (
        #         ggplot(df, aes(x="num_features", y="importance", fill="features"))
        #               + geom_bar(stat="identity")
        #               + xlab("Number of Features")
        #               + ylab("")
        #         )
        # ggsave(plot, os.path.join(figure_path, "stacked_feature_importances.png"))
            # feat_importance_idx = np.argsort(feature_importances)[::-1]
            # sorted_features = features[feat_importance_idx]
            # sorted_feature_importances = sorted(feature_importances, reverse=True)

config = parser.parse_args()
backwards_elim_dirs = get_relevant_backwards_elim_dirs(config)
num_iter_skips = config.num_iter_skips
iters_dont_skip = config.iters_dont_skip
ML_model = config.ML_model
prep_df_for_feat_importance_plots(backwards_elim_dirs, num_iter_skips, iters_dont_skip, ML_model)
