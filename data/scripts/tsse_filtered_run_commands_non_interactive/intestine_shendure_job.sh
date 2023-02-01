#!/bin/bash
#$ -cwd
#$ -l h_rt=48:00:00
#$ -l os=RedHat7
#$ -pe smp 1 # SET THIS
#$ -l h_vmem=192G # SET THIS
#$ -N intestine_shendure
#$ -o /ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/intestine_shendure_job.out # SET THIS
#$ -e /ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/intestine_shendure_job.err # SET THIS
#$ -binding linear:1 # SET THIS
#$ -M bgiotti@broadinstitute.org 
#$ -m bea

source /broad/software/scripts/useuse
use .anaconda3-5.3.1
source activate /home/unix/bgiotti/conda/coo

sh /ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/tsse_filtered_run_commands/intestine_shendure.sh
