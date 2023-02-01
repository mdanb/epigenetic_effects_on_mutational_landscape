#!/bin/bash
#$ -cwd
#$ -l h_rt=24:00:00
#$ -l os=RedHat7
#$ -pe smp 1 # SET THIS
#$ -l h_vmem=32G # SET THIS
#$ -N distal_lung_tsankov
#$ -o /ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/distal_lung_tsankov_job.out # SET THIS
#$ -e /ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/distal_lung_tsankov_job.err # SET THIS
#$ -binding linear:4 # SET THIS
#$ -M bgiotti@broadinstitute.org 
#$ -m bea

source /broad/software/scripts/useuse
use .anaconda3-5.3.1
source activate /home/unix/bgiotti/conda/coo

sh /ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/tsse_filtered_run_commands/distal_lung_tsankov.sh
