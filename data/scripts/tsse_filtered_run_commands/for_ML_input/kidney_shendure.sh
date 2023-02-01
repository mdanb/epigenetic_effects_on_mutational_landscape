#!/bin/bash
#$ -cwd
#$ -l h_rt=24:00:00
#$ -l os=RedHat7
#$ -pe smp 1
#$ -l h_vmem=512G
#$ -N kidney_shendure
#$ -o /ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/tsse_filtered_run_commands/for_ML_input/log_dir/kidney_shendure.out
#$ -e /ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/tsse_filtered_run_commands/for_ML_input/log_dir/kidney_shendure.err
#$ -binding linear:1
#$ -M bgiotti@broadinstitute.org
#$ -m bea

source /broad/software/scripts/useuse
use .anaconda3-5.3.1
source activate /home/unix/bgiotti/conda/coo

Rscript /ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/create_tsse_filtered_count_overlaps.R --cell_types='Metanephric cells;Ureteric bud cells;Unknown.7;Mesangial cells;Vascular endothelial cells;Mesangial cells?;Myeloid cells;Stromal cells?;Unknown.14' --dataset='shendure' --top_tsse_fragment_count_range='100000,150000,200000' --files_pattern='kidney' --cores=1
