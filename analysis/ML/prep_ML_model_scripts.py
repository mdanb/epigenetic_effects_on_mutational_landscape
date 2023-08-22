import argparse
import subprocess

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
parser.add_argument('--fold_for_test_set_range', type=str)
parser.add_argument('--n_optuna_trials_prebackward_selection', type=str)
parser.add_argument('--n_optuna_trials_backward_selection', type=str)
parser.add_argument("--save_test_set_perf", action="store_true", default=False)
parser.add_argument('--test_set_perf_num_features', nargs="+", type=int)
parser.add_argument('--cores', type=str, default="8")
parser.add_argument('--time', type=str, default="24:00:00")
parser.add_argument('--feature_importance_method', type=str)
parser.add_argument("--make_plots", action="store_true", default=False)
parser.add_argument("--create_bash_scripts", action="store_true", default=False)
parser.add_argument("--cleanup", action="store_true", default=False)

group = parser.add_mutually_exclusive_group()
group.add_argument("--meso", action="store_true", default=False)
group.add_argument("--SCLC", action="store_true", default=False)
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
cores = config.cores
time = config.time
top_features_to_plot = config.top_features_to_plot
feature_importance_method = config.feature_importance_method
make_plots = config.make_plots
save_test_set_perf = config.save_test_set_perf
test_set_perf_num_features = config.test_set_perf_num_features
submit_jobs = config.submit_jobs
cleanup = config.cleanup
run_locally = config.run_locally

for fold in fold_for_test_set_range:
    for seed_range in seed_ranges:
        script_filename = "_".join(["cancer_types", "_".join(config.cancer_types),
                                    "datasets", "_".join(sorted(config.datasets)),
                                    "scATAC_cell_number_filter", scATAC_cell_number_filter,
                                    "annotation_dir", annotation_dir,
                                    "seed_range", seed_range, "fold_for_test_set", str(fold)])
                                    # "feature_importance_method", feature_importance_method])
                                        # "top_features_to_plot", "_".join(config.top_features_to_plot),

        if meso:
            script_filename = script_filename + "_" + "meso"
        elif SCLC:
            script_filename = script_filename + "_" + "SCLC"

        script_filename = f"{script_filename}.sh"
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


        if meso:
            command_args = command_args + " " + "--meso"
        elif SCLC:
            command_args = command_args + " " + "--SCLC"

        if make_plots:
            command_args = command_args + " " + "--make_plots"
        if save_test_set_perf:
            command_args = command_args + " " + "--save_test_set_perf"
            command_args = command_args + " " + "--test_set_perf_num_features" + " " + " ".join(map(str,
                                                                                        test_set_perf_num_features))

        print(command_args)
        python_command = "python3 /broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML/build_ML_model.py " + \
                         command_args

        job_script = "\n".join(["#!/bin/bash",
                                "#$ -cwd",
                                f"#$ -l h_rt={time}",
                                "#$ -l os=RedHat7",
                                f"#$ -pe smp {cores}",
                                "#$ -l h_vmem=8G",
                                f"#$ -N {config.cancer_types[0]}_seed_range_{seed_range}_fold_for_test_{fold}",
                                "#$ -o /broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML",
                                "#$ -e /broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML",
                                f"#$ -binding linear:{cores}",
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

        if submit_jobs:
            subprocess.run(["qsub", f"{script_filename}"])

        if run_locally:
            subprocess.run(["sh", f"{script_filename}"])

        if cleanup:
            subprocess.run(["rm", f"{script_filename}"])


