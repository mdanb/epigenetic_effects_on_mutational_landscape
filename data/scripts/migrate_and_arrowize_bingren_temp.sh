#!/bin/bash
#$ -cwd
#$ -l h_rt=48:00:00
#$ -l os=RedHat7
#$ -pe smp 4 
#$ -l h_vmem=96G
#$ -N bingren_migrate_and_arrowize
#$ -o /broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/migrate_and_arrowize_bingren.out
#$ -e /broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/migrate_and_arrowize_bingren.err
#$ -binding linear:4
#$ -M bgiotti@broadinstitute.org
#$ -m bea

source /broad/software/scripts/useuse
use .anaconda3-5.3.1
source activate /home/unix/bgiotti/conda/coo

Rscript /ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/data/scripts/migrate_and_save_fragments.R --dataset=Bingren --cores=4
Rscript /broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/create_arrow_files_and_tss.R --dataset=Bingren
