#!/bin/bash
#$ -cwd
#$ -l h_rt=336:00:00
#$ -l os=RedHat7
#$ -pe smp 1
#$ -l h_vmem=64G
#$ -N combine_ovs_per_cell
#$ -o /broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/scripts
#$ -e /broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/scripts
#$ -binding linear:8

source /broad/software/scripts/useuse
use Anaconda
source activate /home/unix/bgiotti/conda/coo

Rscript combine_overlaps.R --datasets=Greenleaf_colon --annotation=default_annotation --which_interval_ranges=polak --overlaps_per_cell
