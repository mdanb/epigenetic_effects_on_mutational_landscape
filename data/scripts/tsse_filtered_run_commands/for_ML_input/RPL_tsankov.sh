#!/bin/bash
#$ -cwd
#$ -l h_rt=48:00:00
#$ -l os=RedHat7
#$ -pe smp 1
#$ -l h_vmem=128G
#$ -N RPL_tsankov
#$ -o /ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/data/scripts/tsse_filtered_run_commands/for_ML_input/log_dir/RPL_tsankov.out
#$ -e /ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/data/scripts/tsse_filtered_run_commands/for_ML_input/log_dir/RPL_tsankov.err
#$ -binding linear:1
#$ -M bgiotti@broadinstitute.org
#$ -m bea

source /broad/software/scripts/useuse
use .anaconda3-5.3.1
source activate /home/unix/bgiotti/conda/coo

Rscript /ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/data/scripts/create_tsse_filtered_count_overlaps.R --cell_types='SmoothMuscle;Fibroblasts;AT2;Immune;AT1;Endothelial;Ciliated;Mesothelium;Secretory;B_cells' --dataset='tsankov' --top_tsse_fragment_count_range='100000,150000,200000,300000,400000,500000,1000000' --files_pattern='RPL' --cores=1
