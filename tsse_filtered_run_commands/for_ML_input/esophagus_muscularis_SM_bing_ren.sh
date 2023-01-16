#!/bin/bash
#$ -cwd
#$ -l h_rt=24:00:00
#$ -l os=RedHat7
#$ -pe smp 1
#$ -l h_vmem=128G
#$ -N esophagus_muscularis_SM_bing_ren
#$ -o /ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/tsse_filtered_run_commands/for_ML_input/log_dir/esophagus_muscularis_SM_bing_ren.out
#$ -e /ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/tsse_filtered_run_commands/for_ML_input/log_dir/esophagus_muscularis_SM_bing_ren.err
#$ -binding linear:1
#$ -M bgiotti@broadinstitute.org
#$ -m bea

source /broad/software/scripts/useuse
use .anaconda3-5.3.1
source activate /home/unix/bgiotti/conda/coo

Rscript /ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/create_tsse_filtered_count_overlaps.R --cell_types='Fibroblast (General);Fibroblast (Epithelial);Fibroblast (Gastrointestinal);Mesothelial Cell;Fibroblast (Peripheral Nerve);Fibroblast (Sk Muscle Associated);Satellite Cell;Fibroblast (Liver Adrenal);Thyroid Follicular Cell;Pancreatic Acinar Cell;T Lymphocyte 1 (CD8+);T lymphocyte 2 (CD4+);Natural Killer T Cell;Naive T cell;Pericyte (General) 1;Pericyte (General) 2;Pericyte (General) 3;Cardiac Pericyte 2;Pericyte (Esophageal Muscularis);Pericyte (General) 4;Cardiac Pericyte 3;Foveolar Cell;Schwann Cell (General);Macrophage (General);Macrophage (General,Alveolar);GABAergic Neuron 1;Luteal Cell (Ovarian);Plasma Cell;Mast Cell;Endothelial Cell (General) 1;Endothelial Cell (Myocardial);Lymphatic Endothelial Cell;Endothelial Cell (General) 2;Alveolar Capillary Endothelial Cell;Smooth Muscle (Esophageal Muscularis) 1;Smooth Muscle (Colon) 1;Smooth Muscle (Esophageal Muscularis) 2;Smooth Muscle (Uterine);Smooth Muscle (General);Smooth Muscle (GE Junction);Smooth Muscle (Colon) 2;Smooth Muscle (Esophageal Muscularis) 3;Smooth Muscle (General Gastrointestinal);Smooth Muscle (Vaginal);Smooth Muscle (Esophageal Mucosal);Type I Skeletal Myocyte;Type II Skeletal Myocyte;Vascular Smooth Muscle 2;Cardiac Fibroblasts;Peripheral Nerve Stromal;Small Intestinal Enterocyte' --dataset='bing_ren' --top_tsse_fragment_count_range='100000,150000,200000,300000,400000,500000,1000000' --cores=1
