#!/bin/bash
#$ -cwd
#$ -l h_rt=24:00:00
#$ -l os=RedHat7
#$ -pe smp 1
#$ -l h_vmem=32G
#$ -N liver_SM_bing_ren
#$ -o /ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/tsse_filtered_run_commands/for_ML_input/log_dir/liver_SM_bing_ren.out
#$ -e /ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/tsse_filtered_run_commands/for_ML_input/log_dir/liver_SM_bing_ren.err
#$ -binding linear:1
#$ -M bgiotti@broadinstitute.org
#$ -m bea

source /broad/software/scripts/useuse
use .anaconda3-5.3.1
source activate /home/unix/bgiotti/conda/coo

Rscript /ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/create_tsse_filtered_count_overlaps.R --cell_types='Fibroblast (General);Fibroblast (Liver Adrenal);T Lymphocyte 1 (CD8+);Naive T cell;Pericyte (General) 1;Pericyte (General) 3;Pericyte (General) 4;Cardiac Pericyte 3;Cardiac Pericyte 4;Mammary Luminal Epithelial Cell 2;Schwann Cell (General);Macrophage (General);Macrophage (General,Alveolar);Hepatocyte;Ductal Cell (Pancreatic);Mast Cell;Endothelial Cell (General) 1;Endothelial Cell (Myocardial);Lymphatic Endothelial Cell;Endothelial Cell (General) 2;Endothelial Cell (General) 3;Endocardial Cell;Endothelial (Exocrine Tissues);Cardiac Fibroblasts' --dataset='bing_ren' --top_tsse_fragment_count_range='100000,150000,200000' --files_pattern='liver_SM' --cores=1
