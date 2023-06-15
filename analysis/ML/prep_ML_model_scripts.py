import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('--cancer_types', nargs="+", type=str,
                    help='which cancer types to analyze', default=None)
parser.add_argument('--datasets', nargs="+", type=str,
                    help='which sc-ATACseq datasets to analyze', required=True)
parser.add_argument('--scATAC_cell_number_filter', type=str,
                    help='minimum number of cells per cell type in scATAC', default=100)
parser.add_argument('--annotation_dir', type=str,
                    help='name of annotation directory', default="default_annotation")
parser.add_argument('--seed_interval', type=str)
parser.add_argument('--seed_interval_step', type=int)
parser.add_argument('--iters_dont_skip', nargs="+", type=str,
                    default=["18"])
parser.add_argument('--n_optuna_trials_prebackward_selection', type=str)
parser.add_argument('--n_optuna_trials_backward_selection', type=str)
group = parser.add_mutually_exclusive_group()
group.add_argument("--meso", action="store_true", default=False)

config = parser.parse_args()
cancer_types = " ".join(config.cancer_types)
datasets = " ".join(sorted(config.datasets))
scATAC_cell_number_filter = config.scATAC_cell_number_filter
annotation_dir = config.annotation_dir
iters_dont_skip = " ".join(config.iters_dont_skip)
seed_interval = config.seed_interval
start, end = map(int, seed_interval.split('-'))
seed_interval_step = config.seed_interval_step
seed_ranges = list(range(start - 1, end + 1, seed_interval_step))
seed_ranges = [f"{seed_ranges[i] + 1}-{seed_ranges[i+1]}" for i in range(len(seed_ranges)-1)]
n_optuna_trials_prebackward_selection = config.n_optuna_trials_prebackward_selection
n_optuna_trials_backward_selection = config.n_optuna_trials_backward_selection
meso = config.meso

for seed_range in seed_ranges:
    script_filename = "_".join(["cancer_types", cancer_types,
                                "datasets", datasets,
                                "scATAC_cell_number_filter", scATAC_cell_number_filter,
                                "annotation_dir", annotation_dir,
                                "seed_range", seed_range,
                                "iters_dont_skip", iters_dont_skip,
                                "n_optuna_trials_prebackward_selection", n_optuna_trials_prebackward_selection,
                                "n_optuna_trials_backward_selection", n_optuna_trials_backward_selection])
    if meso:
        script_filename = script_filename + "_" + "meso"
    script_filename = f"{script_filename}.sh"
    command_args = " ".join(["--cancer_types", cancer_types,
                             "--datasets", datasets,
                             "--scATAC_cell_number_filter", scATAC_cell_number_filter,
                             "--annotation_dir", annotation_dir,
                             "--seed_range", seed_range,
                             "--iters_dont_skip", iters_dont_skip,
                             "--n_optuna_trials_prebackward_selection", n_optuna_trials_prebackward_selection,
                             "--n_optuna_trials_backward_selection", n_optuna_trials_backward_selection])

    if meso:
        command_args = command_args + " " + "--meso"

    python_command = "python3 /broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML/build_ML_model.py " + \
                     command_args

    job_script = "\n".join(["#!/bin/bash",
                            "#$ -cwd",
                            "#$ -l h_rt=24:00:00",
                            "#$ -l os=RedHat7",
                            "#$ -pe smp 8",
                            "#$ -l h_vmem=8G",
                            f"#$ -N seed_range_{seed_range}",
                            "#$ -o /broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML",
                            "#$ -e /broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML",
                            "#$ -binding linear:8",
                            "#$ -M bgiotti@broadinstitute.org",
                            "#$ -m bea",
                            "",
                            "source /broad/software/scripts/useuse",
                            "use Anaconda",
                            "source activate /home/unix/bgiotti/conda/coo",
                            "",
                            python_command])
    with open(script_filename, "w") as f:
        f.write(job_script)

    # subprocess.run(["qsub", f"{script_filename}"])

