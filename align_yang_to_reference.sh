#!/bin/bash

/seq/regev_genome_portal/SOFTWARE/10X/cellranger-atac-1.2.0/cellranger-atac count --jobmode=local --localcores=8 --localmem=64 --id=SRR13679156 \
--reference=/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/hg19 \
--fastqs=/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/bed_files/yang_kidney_scATAC/SRR13679156.fq
--sample=SRR13679156
