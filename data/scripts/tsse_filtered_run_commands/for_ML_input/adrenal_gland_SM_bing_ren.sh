#!/bin/bash
#$ -cwd
#$ -l h_rt=48:00:00
#$ -l os=RedHat7
#$ -pe smp 1
#$ -l h_vmem=32G
#$ -N adrenal_gland_SM_bing_ren
#$ -o /ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/tsse_filtered_run_commands/for_ML_input/log_dir/adrenal_gland_SM_bing_ren.out
#$ -e /ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/tsse_filtered_run_commands/for_ML_input/log_dir/adrenal_gland_SM_bing_ren.err
#$ -binding linear:1
#$ -M bgiotti@broadinstitute.org
#$ -m bea

source /broad/software/scripts/useuse
use .anaconda3-5.3.1
source activate /home/unix/bgiotti/conda/coo

Rscript /ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/create_tsse_filtered_count_overlaps.R --cell_types='Fibroblast (General);Adipocyte;Fibroblast (Peripheral Nerve);Fibroblast (Liver Adrenal);T Lymphocyte 1 (CD8+);T lymphocyte 2 (CD4+);Natural Killer T Cell;Naive T cell;Pericyte (General) 1;Pericyte (General) 3;Pericyte (General) 4;Cardiac Pericyte 3;Cardiac Pericyte 4;Schwann Cell (General);Macrophage (General);Macrophage (General,Alveolar);CNS,Enteric Neuron;Transitional Zone Cortical Cell;Zona Fasciculata Cortical Cell;Zona Glomerulosa Cortical Cell;Cortical Epithelial-like;Luteal Cell (Ovarian);Memory B Cell;Mast Cell;Endothelial Cell (General) 1;Lymphatic Endothelial Cell;Endothelial Cell (General) 2;Endothelial (Exocrine Tissues);Smooth Muscle (General);Smooth Muscle (Vaginal);Vascular Smooth Muscle 2' --dataset='bing_ren' --top_tsse_fragment_count_range='100000,150000,200000,300000,400000,500000,1000000' --cores=1 --files_pattern='adrenal_gland_SM'
