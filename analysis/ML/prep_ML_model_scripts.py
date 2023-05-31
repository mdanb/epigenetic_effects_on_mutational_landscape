import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('--cancer_types', nargs="+", type=str,
                    help='which cancer types to analyze', default=None)
parser.add_argument('--datasets', nargs="+", type=str,
                    help='which sc-ATACseq datasets to analyze', required=True)
parser.add_argument('--scATAC_cell_number_filter', type=int,
                    help='minimum number of cells per cell type in scATAC', default=100)
parser.add_argument('--annotation_dir', type=str,
                    help='name of annotation directory', default="default_annotation")
parser.add_argument('--seed_ranges', type=str)

config = parser.parse_args()
cancer_types = config.cancer_types
datasets = sorted(config.datasets)
scATAC_cell_number_filter = config.scATAC_cell_number_filter
annotation_dir = config.annotation_dir

ranges = config.seed_ranges.split(',')

for seed_range in ranges:
    script_filename = "_".join(["cancer_types", cancer_types,
                                "datasets", datasets,
                                "scATAC_cell_number_filter", scATAC_cell_number_filter,
                                "annotation_dir", annotation_dir,
                                "seed_range", seed_range])
    script_filename = f"{script_filename}.sh"
    command_args = " ".join(["--cancer_types", cancer_types,
                             "--datasets", datasets,
                             "--scATAC_cell_number_filter", scATAC_cell_number_filter,
                             "--annotation_dir", annotation_dir,
                             "--seed_range", seed_range])
    python_command = "python /broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML/build_ML_model.py " + \
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
                            "use .anaconda3-5.3.1",
                            "source activate /home/unix/bgiotti/conda/coo",
                            "",
                            python_command])
    with open(script_filename, "w") as f:
        f.write(job_script)

    # subprocess.run(["qsub", f"{script_filename}"])

