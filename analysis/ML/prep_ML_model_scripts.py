import argparse
import subprocess
import re
import os
from config import int_or_all

parser = argparse.ArgumentParser()
parser.add_argument('--cancer_types', nargs="+", type=str,
                    help='which cancer types to analyze', default=None)
parser.add_argument('--datasets', nargs="+", type=str,
                    help='which sc-ATACseq datasets to analyze', required=True)
parser.add_argument('--scATAC_cell_number_filter', type=str,
                    help='minimum number of cells per cell type in scATAC', default="100")
parser.add_argument('--annotation_dir', type=str,
                    help='name of annotation directory', default="default_annotation")
parser.add_argument('--seed_interval', type=str)
parser.add_argument('--seed_interval_step', type=int)
parser.add_argument('--top_features_to_plot', nargs="+", type=str,
                    default=["18"])
parser.add_argument('--top_features_to_plot_feat_imp', nargs="+", type=str)
parser.add_argument('--fold_for_test_set_range', type=str)
parser.add_argument('--n_optuna_trials_prebackward_selection', type=str)
parser.add_argument('--n_optuna_trials_backward_selection', type=str)
# parser.add_argument("--save_test_set_perf", action="store_true", default=False)
parser.add_argument('--test_set_perf_num_features', nargs="+", type=int_or_all)
parser.add_argument('--cores', type=str, default="8")
parser.add_argument('--mem_per_core', type=str, default="8")
parser.add_argument('--time', type=str, default="24:00:00")
parser.add_argument('--feature_importance_method', type=str)
parser.add_argument("--make_plots", action="store_true", default=False)
parser.add_argument("--create_bash_scripts", action="store_true", default=False)
parser.add_argument("--cleanup", action="store_true", default=False)
parser.add_argument('--feat_imp_min_n_robustness', type=int)
parser.add_argument('--grid_cell_types', nargs="+", type=str, default=None)
parser.add_argument('--cell_types_keep', nargs="+", type=str, default=None)
parser.add_argument('--tss_fragment_filter', type=str,
                    help='tss fragment filter to consider', default=None)
parser.add_argument('--which_interval_ranges', type=str, default=None)

group = parser.add_mutually_exclusive_group()
group.add_argument("--meso", action="store_true", default=False)
group.add_argument("--SCLC", action="store_true", default=False)
group.add_argument("--de_novo_seurat_clustering", action="store_true", default=False)
group.add_argument("--histologically_subtyped_mutations", action="store_true", default=False)
group.add_argument("--hierarchically_subtyped_mutations", action="store_true", default=False)
group.add_argument("--hundred_kb", action="store_true", default=False)
group.add_argument("--expanded_hundred_kb", action="store_true", default=False)
group.add_argument("--mm", action="store_true", default=False)
group.add_argument("--woo_pcawg", action="store_true", default=False)
group.add_argument("--msi_high", action="store_true", default=False)
group.add_argument("--subsampled_mutations", action="store_true", default=False)
group.add_argument("--custom_mutations", action="store_true", default=False)

group2 = parser.add_mutually_exclusive_group()
group2.add_argument('--grid_analysis', action="store_true", default=False)
group2.add_argument('--tissues_to_consider', nargs="+", type=str, default=["all"])

group_run_loc = parser.add_mutually_exclusive_group()
group_run_loc.add_argument("--submit_jobs", action="store_true", default=False)
group_run_loc.add_argument("--run_locally", action="store_true", default=False)

config = parser.parse_args()
scATAC_cell_number_filter = config.scATAC_cell_number_filter
annotation_dir = config.annotation_dir
# iters_dont_skip = " ".join(config.iters_dont_skip)
seed_interval = config.seed_interval
start, end = map(int, seed_interval.split('-'))
seed_interval_step = config.seed_interval_step
seed_ranges = list(range(start - 1, end + 1, seed_interval_step))
seed_ranges = [f"{seed_ranges[i] + 1}-{seed_ranges[i+1]}" for i in range(len(seed_ranges)-1)]
start, end = map(int, config.fold_for_test_set_range.split('-'))
fold_for_test_set_range = range(start, end + 1)
create_bash_scripts = config.create_bash_scripts
n_optuna_trials_prebackward_selection = config.n_optuna_trials_prebackward_selection
n_optuna_trials_backward_selection = config.n_optuna_trials_backward_selection
meso = config.meso
SCLC = config.SCLC
de_novo_seurat_clustering = config.de_novo_seurat_clustering
histologically_subtyped_mutations = config.histologically_subtyped_mutations
hierarchically_subtyped_mutations = config.hierarchically_subtyped_mutations
mm = config.mm
woo_pcawg = config.woo_pcawg
hundred_kb = config.hundred_kb
cores = config.cores
mem_per_core = config.mem_per_core
time = config.time
top_features_to_plot = config.top_features_to_plot
feature_importance_method = config.feature_importance_method
make_plots = config.make_plots
# save_test_set_perf = config.save_test_set_perf
test_set_perf_num_features = config.test_set_perf_num_features
submit_jobs = config.submit_jobs
cleanup = config.cleanup
run_locally = config.run_locally
expanded_hundred_kb = config.expanded_hundred_kb
tissues_to_consider = config.tissues_to_consider
top_features_to_plot_feat_imp = config.top_features_to_plot_feat_imp
feat_imp_min_n_robustness = config.feat_imp_min_n_robustness
grid_analysis = config.grid_analysis
grid_cell_types = config.grid_cell_types
msi_high = config.msi_high
cell_types_keep = config.cell_types_keep
tss_fragment_filter = config.tss_fragment_filter
subsampled_mutations = config.subsampled_mutations
custom_mutations = config.custom_mutations
which_interval_ranges = config.which_interval_ranges

for fold in fold_for_test_set_range:
    for seed_range in seed_ranges:
        # cancer_types = config.cancer_types
        # match = re.search(r'cluster_[0-9]+', config.cancer_types[0])
        # if match:
        #     cancer_types = match.group()
        script_filename = "_".join(["cancer_types", "_".join(config.cancer_types),
                                    "datasets", "_".join(sorted(config.datasets)),
                                    "scATAC_cell_number_filter", scATAC_cell_number_filter,
                                    "annotation_dir", annotation_dir,
                                    "tissues_to_consider", "_".join(tissues_to_consider)])
                                    # "feature_importance_method", feature_importance_method])
                                        # "top_features_to_plot", "_".join(config.top_features_to_plot),

        command_args = " ".join(["--cancer_types", " ".join(config.cancer_types),
                                 "--datasets", " ".join(sorted(config.datasets)),
                                 "--scATAC_cell_number_filter", scATAC_cell_number_filter,
                                 "--annotation_dir", annotation_dir,
                                 "--seed_range", seed_range,
                                 "--top_features_to_plot", " ".join(config.top_features_to_plot),
                                 "--n_optuna_trials_prebackward_selection", n_optuna_trials_prebackward_selection,
                                 "--n_optuna_trials_backward_selection", n_optuna_trials_backward_selection,
                                 "--feature_importance_method", feature_importance_method,
                                 "--fold_for_test_set", str(fold)])

        if cell_types_keep:
            ctk = " ".join(cell_types_keep)
            command_args = command_args + " " + "--cell_types_keep" + " " + f"'{ctk}'"

        if meso:
            script_filename = script_filename + "_" + "meso"
            command_args = command_args + " " + "--meso"
        elif SCLC:
            script_filename = script_filename + "_" + "SCLC"
            command_args = command_args + " " + "--SCLC"
        elif de_novo_seurat_clustering:
            script_filename = script_filename + "_" + "de_novo_seurat_clustering"
            command_args = command_args + " " + "--de_novo_seurat_clustering"
        elif histologically_subtyped_mutations:
            # script_filename = script_filename + "_" + "histologically_subtyped_mutations"
            command_args = command_args + " " + "--histologically_subtyped_mutations"
        elif hierarchically_subtyped_mutations:
            command_args = command_args + " " + "--hierarchically_subtyped_mutations"
        elif hundred_kb:
            script_filename = script_filename + "_" + "hundred_kb"
            command_args = command_args + " " + "--hundred_kb"
        elif expanded_hundred_kb:
            script_filename = script_filename + "_" + "expanded_hundred_kb"
            command_args = command_args + " " + "--expanded_hundred_kb"
        elif mm:
            script_filename = script_filename + "_" + "mm"
            command_args = command_args + " " + "--mm"
        elif woo_pcawg:
            script_filename = script_filename + "_" + "woo_pcawg"
            command_args = command_args + " " + "--woo_pcawg"
        elif msi_high:
            script_filename = script_filename + "_" + "msi_high"
            command_args = command_args + " " + "--msi_high"
        elif subsampled_mutations:
            command_args = command_args + " " + "--subsampled_mutations"
        elif custom_mutations:
            command_args = command_args + " " + "--custom_mutations"

        if grid_analysis:
            command_args = command_args + " " + "--grid_analysis" + " " + "--grid_cell_types" + " " + \
                           " ".join([f'"{item}"' for item in grid_cell_types])
        else:
            command_args = command_args + " " + "--tissues_to_consider" + " " + " ".join(tissues_to_consider)

        if tss_fragment_filter:
            command_args = command_args + " " + "--tss_fragment_filter" + " " + tss_fragment_filter
            script_filename = script_filename + "_" + "tss_filter" + "_" + tss_fragment_filter

        robustness_filename = script_filename
        script_filename = "_".join([script_filename, "seed_range", seed_range, "fold_for_test_set", str(fold)])

        script_filename = f"{script_filename}.sh"

        if make_plots:
            command_args = command_args + " " + "--make_plots"
        # if save_test_set_perf:
        #     command_args = command_args + " " + "--save_test_set_perf"
        command_args = command_args + " " + "--test_set_perf_num_features" + " " + " ".join(map(str,
                                                                                    test_set_perf_num_features))
        if which_interval_ranges:
            command_args = command_args + " " + "--which_interval_ranges" + " " + which_interval_ranges

        python_command = "python3 /broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML/build_ML_model.py " + \
                         command_args

        job_script = "\n".join(["#!/bin/bash",
                                "#$ -cwd",
                                f"#$ -l h_rt={time}",
                                "#$ -l os=RedHat7",
                                f"#$ -pe smp {cores}",
                                f"#$ -l h_vmem={mem_per_core}G",
                                f"#$ -N {config.cancer_types[0]}_seed_range_{seed_range}_fold_for_test_{fold}",
                                "#$ -o /broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML",
                                "#$ -e /broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML",
                                f"#$ -binding linear:{cores}"])
                                # "",
                                # "source /broad/software/scripts/useuse",
                                # "use Anaconda",
                                # "source activate /ahg/regevdata/projects/ICA_Lung/Wooseung/conda/coo",
                                # "",
                                # python_command])

        if submit_jobs:
            job_script = '\n'.join([job_script,
                                   "",
                                   "source /broad/software/scripts/useuse",
                                   "use Anaconda",
                                   "source activate /ahg/regevdata/projects/ICA_Lung/Wooseung/conda/coo",
                                   ""])
        job_script = '\n'.join([job_script,
                                python_command])

        try:
            print("Creating scripts...")
            with open(script_filename, "w") as f:
                f.write(job_script)
        except:
            if "cluster" in script_filename:
                name = re.search("cluster_[0-9]+", script_filename).group()
                script_filename = name + "_" + re.search("seed_range_.*", script_filename).group()
            else:
            	script_filename = "temp.sh"

            with open(script_filename, "w") as f:
                 print("Creating scripts...")
                 f.write(job_script)

        if submit_jobs:
            subprocess.run(["qsub", f"{script_filename}"])
        elif run_locally:
            subprocess.run(["sh", f"{script_filename}"])

        if cleanup:
            subprocess.run(["rm", f"{script_filename}"])

os.makedirs("robustness_scripts", exist_ok=True)
robustness_filename = f"robustness_{robustness_filename}.sh"
robustness_fp = os.path.join("robustness_scripts", robustness_filename)
command_args = " ".join(["--cancer_types", ",".join(config.cancer_types),
                         "--datasets", ",".join(sorted(config.datasets)),
                         "--cell_number_filter", scATAC_cell_number_filter,
                         "--annotation", annotation_dir,
                         "--seed_range 1-10",
                         "--top_features_to_plot_feat_imp", ",".join(config.top_features_to_plot_feat_imp),
                         "--top_features_to_plot", ",".join(config.top_features_to_plot),
                         "--feature_importance_method", feature_importance_method,
                         "--folds_for_test_set 1-10",
                         "--robustness_analysis",
                         "--feat_imp_min_n_robustness", str(feat_imp_min_n_robustness)])

if grid_analysis:
    command_args = command_args + " " + "--grid_analysis" + " " + "--grid_cell_types" + " " + \
                   " ".join([f'"{item}"' for item in grid_cell_types])
else:
    command_args = command_args + " " + "--tissues_to_consider" + " " + ",".join(tissues_to_consider)

if cell_types_keep:
    cell_types_keep = ",".join(cell_types_keep)
    command_args = command_args + " " + "--cell_types_keep" + " " + f"'{cell_types_keep}'"

if tss_fragment_filter:
    command_args = command_args + " " + "--tss_fragment_filter" + " " + tss_fragment_filter

if hundred_kb:
    command_args = command_args + " " + "--hundred_kb"

Rscript_command = "Rscript /broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML/plot_top_features.R " + \
                  command_args

try:
    with open(robustness_fp, "w") as f:
        f.write(Rscript_command)
except:
    print("exception")
    robustness_fp = "robustness_scripts/temp.sh"
    print(robustness_fp)
    with open(robustness_fp, "w") as f:
        f.write(Rscript_command)




