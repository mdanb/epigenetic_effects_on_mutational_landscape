import pickle
import argparse
import glob
from natsort import natsorted
import pandas as pd
import os
from pathlib import Path
import optuna
import sys

# sys.path.insert(0, '/broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML')
sys.path.insert(0, '/home/mdanb/research/mount_sinai/epigenetic_effects_on_mutational_landscape')
sys.path.insert(0, '/broad/hptmp/bgiotti/BingRen_scATAC_atlas/')

from analysis.ML.ML_utils import construct_scATAC_dir, construct_scATAC_sources
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
parser.add_argument('--tissues_to_consider', nargs="+", type=str, default=["all"])
parser.add_argument('--iters_dont_skip', nargs="+", type=int)
parser.add_argument('--ML_model', type=str)
parser.add_argument('--cell_number_filter', type=int)
parser.add_argument('--num_iter_skips', type=int, default=5)
parser.add_argument('--seed', type=int, default=42)
parser.add_argument('--tss_fragment_filter', nargs="+", type=int, default=[-1])


def construct_backwards_elim_dir(cancer_type, scATAC_source, cell_number_filter,
                                 tss_fragment_filter, annotation, tissues_to_consider,
                                 ML_model, seed):
    dir = f"/broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML/" \
          f"models/{ML_model}/{cancer_type}/scATAC_source_" \
          f"{scATAC_source}_cell_number_filter_{cell_number_filter}"
    if (tss_fragment_filter != -1):
        dir = dir + f"_tss_fragment_filter_{tss_fragment_filter}"
    dir = dir + f"_annotation_{annotation}_seed_{seed}"

    if (tissues_to_consider == "all"):
        return f"{dir}/backwards_elimination_results/"
    else:
        return f"{dir}/backwards_elimination_results_{tissues_to_consider}"

def get_relevant_backwards_elim_dirs(config):
    print("cancer types")
    cancer_types = config.cancer_types
    print(cancer_types)
    datasets = config.datasets
    cell_number_filter = config.cell_number_filter
    tss_fragment_filter = config.tss_fragment_filter
    tissues_to_consider = "_".join(config.tissues_to_consider)
    annotation = config.annotation
    ML_model = config.ML_model
    seed = config.seed
    backward_elim_dirs = []

    scATAC_sources = ""
    for idx, dataset in enumerate(datasets):
        if (scATAC_sources == ""):
            scATAC_sources = dataset
        else:
            scATAC_sources = "_".join((scATAC_sources, dataset))

    for cancer_type in cancer_types:
        for tss_filter in tss_fragment_filter:
            backward_elim_dirs.append(construct_backwards_elim_dir(cancer_type,
                                                                   scATAC_sources,
                                                                   cell_number_filter,
                                                                   tss_filter,
                                                                   annotation,
                                                                   tissues_to_consider,
                                                                   ML_model,
                                                                   seed))
    print(backward_elim_dirs)
    return backward_elim_dirs

def prep_df_for_feat_importance_plots(backwards_elim_dirs, num_iter_skips, iters_dont_skip,
                                      ML_model, optuna_base_dir):
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
        files = natsorted(glob.glob(f"{backwards_elim_dir}/*pkl"))
        for idx, file in enumerate(files):
            if (idx % num_iter_skips == 0 or idx in iters_dont_skip):
                model = pickle.load(open(file, "rb"))
                study_name = f"iter_{idx}_{optuna_base_dir}"
                study = optuna.load_study(study_name=study_name, storage="sqlite:///db.sqlite3")
                best_trial = study.best_trial
                best_cv_score = best_trial.value
                features = model.feature_names_in_
                feature_importances = model.feature_importances_
                df_curr = pd.DataFrame((features, feature_importances,
                                        [len(features)] * len(features),
                                        [best_cv_score] * len(features))).T
                df_curr.columns = df.columns
                df = pd.concat((df, df_curr))
        Path(figure_path).mkdir(parents=True, exist_ok=True)
        df.to_csv(os.path.join(figure_path, "df_for_feature_importance_plots.csv"),
                  index=False)

config = parser.parse_args(["--cancer_types", "CNS-GBM", "--datasets", "Bingren", "Greenleaf_brain", "Shendure",
                            "Tsankov", "--cell_number_filter", "1", "--annotation", "finalized_annotation",
                            "--ML_model", "XGB", "--seed", "42"])
backwards_elim_dirs = get_relevant_backwards_elim_dirs(config)
print(backwards_elim_dirs)
num_iter_skips = config.num_iter_skips
iters_dont_skip = config.iters_dont_skip
ML_model = config.ML_model
optuna_base_dir = construct_scATAC_dir(scATAC_sources=construct_scATAC_sources(config.datasets),
                                       scATAC_cell_number_filter=config.cell_number_filter,
                                       tss_filter=None,
                                       annotation_dir=config.annotation,
                                       seed=config.seed)
prep_df_for_feat_importance_plots(backwards_elim_dirs, num_iter_skips,
                                  iters_dont_skip, ML_model, optuna_base_dir)
