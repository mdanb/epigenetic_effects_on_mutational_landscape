#!/bin/bash
#$ -cwd
#$ -M bgiotti@login.broadinstitute.org
#$ -l h_rt=96:00:00
#$ -l os=RedHat7
#$ -pe smp 8
#$ -l h_vmem=8G
#$ -N CD45_pos_neg_count_hg38_atac_1.2
#$ -o stdout_CD45_pos_neg_count_hg38_atac_1.2
#$ -e stderr_CD45_pos_neg_count_hg38_atac_1.2
#$ -binding linear:8
#$ -M bgiotti@broadinstitute.org
#$ -m bea source /broad/software/scripts/useuse
use .bcl2fastq2-v2.20.0
cd /ahg/regevdata/projects/lungCancerBueno/10x/200128/scATAC_Pt811_mesothelioma_CD45pos_neg_cellranger_atac_v1.2
/seq/regev_genome_portal/SOFTWARE/10X/cellranger-atac-1.2.0/cellranger-atac count --jobmode=local --localcores=8 --localmem=64 --id=161_ATACseq_Pt811_mesothelioma_CD45pos_CD45neg_ATAC \
--reference=/seq/regev_genome_portal/RESOURCES/10x_genomes/refdata-cellranger-atac-GRCh38-1.2.0 \
--fastqs=/ahg/regevdata/projects/lungCancerBueno/10x/200128/wangy33.u.hpc.mssm.edu/10X_Single_Cell_RNA/TD01218_Robby_AlexTsankov/fastqs/161_ATACseq_Pt811_mesothelioma_CD45pos_CD45neg_ATAC \
--sample=161_ATACseq_Pt811_mesothelioma_CD45pos_CD45neg_ATAC
