#!/bin/bash
#$ -cwd
#$ -l h_rt=36:00:00
#$ -l os=RedHat7
#$ -pe smp 4
#$ -l h_vmem=64G
#$ -N migrate_gl
#$ -o /broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/migrate_gl.out
#$ -e /broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/migrate_gl.err
#$ -binding linear:4
#$ -M bgiotti@broadinstitute.org
#$ -m bea

source /broad/software/scripts/useuse
use Anaconda 
source activate /home/unix/bgiotti/conda/coo

Rscript /broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/scripts/migrate_and_save_fragments.R --dataset=Greenleaf_colon --cores=4
