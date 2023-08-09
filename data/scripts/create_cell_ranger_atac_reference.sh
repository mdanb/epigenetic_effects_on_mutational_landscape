cd ..
/seq/regev_genome_portal/SOFTWARE/10X/cellranger-3.1.0/cellranger mkgtf Homo_sapiens.GRCh37.87.gtf Homo_sapiens.GRCh37.87.filtered.gtf \
	--attribute=gene_biotype:protein_coding \
	--attribute=gene_biotype:lincRNA \
	--attribute=gene_biotype:antisense

/seq/regev_genome_portal/SOFTWARE/10X/cellranger-3.1.0/cellranger mkref --genome=hg19 \
                 --fasta=Homo_sapiens.GRCh37.dna.primary_assembly.fa \
                 --genes=Homo_sapiens.GRCh37.87.filtered.gtf \
                 --ref-version=3.1.0
