#!/bin/bash
#$ -cwd
#$ -l h_rt=24:00:00
#$ -l os=RedHat7
#$ -pe smp 1
#$ -l h_vmem=320G
#$ -N placenta_shendure
#$ -o /ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/tsse_filtered_run_commands/for_ML_input/log_dir/placenta_shendure.out
#$ -e /ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/tsse_filtered_run_commands/for_ML_input/log_dir/placenta_shendure.err
#$ -binding linear:1
#$ -M bgiotti@broadinstitute.org
#$ -m bea

source /broad/software/scripts/useuse
use .anaconda3-5.3.1
source activate /home/unix/bgiotti/conda/coo

Rscript /ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/create_tsse_filtered_count_overlaps.R --cell_types='Trophoblast giant cells;Lymphoid cells;Syncytiotrophoblasts and villous cytotrophoblasts?;Erythroblasts;Syncytiotrophoblasts and villous cytotrophoblasts;Smooth muscle cells;IGFBP1_DKK1 positive cells;Stromal cells;Myeloid cells;Vascular endothelial cells;Extravillous trophoblasts;PAEP_MECOM positive cells' --dataset='shendure' --top_tsse_fragment_count_range='100000,150000,200000' --files_pattern='placenta' --cores=1
