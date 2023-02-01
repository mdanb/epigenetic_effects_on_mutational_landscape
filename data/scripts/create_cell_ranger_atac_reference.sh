#!/bin/bash


/seq/regev_genome_portal/SOFTWARE/10X/cellranger-3.1.0/cellranger mkgtf /broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/Homo_sapiens.GRCh37.87.gtf /broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/Homo_sapiens.GRCh37.87.filtered.gtf \
                 --attribute=gene_biotype:protein_coding \
                 --attribute=gene_biotype:lincRNA \
                 --attribute=gene_biotype:antisense

cd /broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/
/seq/regev_genome_portal/SOFTWARE/10X/cellranger-3.1.0/cellranger mkref --genome=hg19 \
                 --fasta=/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/Homo_sapiens.GRCh37.dna.primary_assembly.fa \
                 --genes=/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/Homo_sapiens.GRCh37.87.filtered.gtf \
                 --ref-version=3.0.0
