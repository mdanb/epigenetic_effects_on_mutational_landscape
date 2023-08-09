cd .. 
cellranger mkgtf Homo_sapiens.GRCh37.87.gtf Homo_sapiens.GRCh37.87.filtered.gtf \
                 --attribute=gene_biotype:protein_coding \
                 --attribute=gene_biotype:lincRNA \
                 --attribute=gene_biotype:antisense
cellranger mkref --genome=hg19 \
                 --fasta=Homo_sapiens.GRCh37.dna.primary_assembly.fa \
                 --genes=Homo_sapiens.GRCh37.87.filtered.gtf \
                 --ref-version=3.0.0
