#!/bin/bash

#$ -cwd
#$ -M bgiotti@login.broadinstitute.org
#$ -l h_rt=48:00:00
#$ -l os=RedHat7
#$ -pe smp 8
#$ -l h_vmem=8G
#$ -N CD45_pos_neg_count_hg38_atac_1.2
#$ -o stdout_CD45_pos_neg_count_hg38_atac_1.2
#$ -e stderr_CD45_pos_neg_count_hg38_atac_1.2
#$ -binding linear:8
#$ -M bgiotti@broadinstitute.org
#$ -m bea

