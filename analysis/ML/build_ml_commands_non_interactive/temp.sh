#!/bin/bash
#$ -cwd
#$ -l h_rt=24:00:00
#$ -l os=RedHat7
#$ -pe smp 1 # SET THIS
#$ -l h_vmem=8G # SET THIS
#$ -N temp
#$ -o /ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/build_ml_models_non_interactive/temp.out # SET THIS
#$ -e /ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/build_ml_models_non_interactive/temp.err # SET THIS
#$ -binding linear:8 # SET THIS
#$ -M bgiotti@broadinstitute.org 
#$ -m bea

source /broad/software/scripts/useuse
use .anaconda3-5.3.1
source activate /home/unix/bgiotti/conda/coo

sh /ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/build_ml_models_non_interactive/temp.sh

