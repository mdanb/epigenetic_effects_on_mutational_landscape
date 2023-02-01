#!/bin/bash

#$ -cwd
#$ -l h_rt=48:00:00
#$ -l os=RedHat7
#$ -pe smp 1 # SET THIS
#$ -l h_vmem=8G # SET THIS
#$ -N SRR13679157
#$ -o /ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/SRR13679157.out # SET THIS
#$ -e /ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/SRR13679157.err # SET THIS
#$ -binding linear:8 # SET THIS
#$ -M bgiotti@broadinstitute.org 
#$ -m bea

source /broad/software/scripts/useuse
use .anaconda3-5.3.1
source activate /home/unix/bgiotti/conda/coo

./software/cellranger-atac-2.1.0/cellranger-atac count \
--jobmode=local \
--localcores=8 \
--localmem=64 \
--id=SRR13679157 \
--reference=/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/GRCh37 \
--fastqs=/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/fastq/SRR13679157 \
