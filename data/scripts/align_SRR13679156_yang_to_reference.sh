#!/bin/bash

#$ -cwd
#$ -l h_rt=48:00:00
#$ -l os=RedHat7
#$ -pe smp 1 # SET THIS
#$ -l h_vmem=16G # SET THIS
#$ -N SRR13679156
#$ -o /broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/bed_files/yang_kidney_scATAC/SRR13679156.out # SET THIS
#$ -e /broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/bed_files/yang_kidney_scATAC/SRR13679156.err # SET THIS
#$ -binding linear:8 # SET THIS
#$ -M bgiotti@broadinstitute.org 
#$ -m bea

source /broad/software/scripts/useuse
use Anaconda
source activate /home/unix/bgiotti/conda/coo

cd ../bed_files/yang_kidney_scATAC

/seq/regev_genome_portal/SOFTWARE/10X/cellranger-atac-1.1.0/cellranger-atac count \
--jobmode=local \
--localcores=8 \
--localmem=16 \
--id=SRR13679156 \
--reference=/seq/regev_genome_portal/RESOURCES/10x_genomes/new/refdata-cellranger-atac-hg19-1.1.0 \
--fastqs=/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/fastq/SRR13679156 \
