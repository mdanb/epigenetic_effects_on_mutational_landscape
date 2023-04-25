import argparse

def range_type(range_str):
    start, end = map(int, range_str.split('-'))
    if start > end:
        raise argparse.ArgumentTypeError('Invalid range: start value must be '
                                         'less than or equal to end value')
    return (start, end)

parser = argparse.ArgumentParser()

parser.add_argument('--cancer_types', nargs="+", type=str,
                    help='which cancer types to analyze', default=None)
parser.add_argument('--clustered_mutations', action="store_true",
                    help='run model on hierarchically clustered mutations', default=False)
parser.add_argument('--datasets', nargs="+", type=str,
                    help='which sc-ATACseq datasets to analyze', required=True)
parser.add_argument('--scATAC_cell_number_filter', type=int,
                    help='minimum number of cells per cell type in scATAC', default=100)
parser.add_argument('--annotation_dir', type=str,
                    help='name of annotation directory', default="default_annotation")
group = parser.add_mutually_exclusive_group()
group.add_argument('--waddell_sarc_biph', action="store_true",
                    default=False)
group.add_argument('--waddell_sarc', action="store_true", default=False)
group.add_argument('--waddell_sarc_tsankov_sarc', action="store_true", default=False)
group.add_argument('--waddell_sarc_biph_tsankov_sarc_biph', action="store_true", default=False)
group.add_argument("--SCLC", action="store_true", default=False)
group.add_argument("--lung_subtyped", action="store_true", default=False)
parser.add_argument("--woo_pcawg", action="store_true", default=False)
parser.add_argument("--histologically_subtyped_mutations", action="store_true", default=False)
parser.add_argument("--de_novo_seurat_clustering", action="store_true", default=False)
parser.add_argument("--per_donor", action="store_true", default=False)
parser.add_argument("--CPTAC", action="store_true", default=False)
parser.add_argument("--combined_CPTAC_ICGC", action="store_true", default=False)
parser.add_argument('--donor_range', type=range_type, help='Specify a range in the format start-end',
                    default=None)
parser.add_argument('--tss_fragment_filter', nargs="+", type=str,
                    help='tss fragment filters to consider', default="")
parser.add_argument('--tissues_to_consider', nargs="+", type=str, default="all")
parser.add_argument("--ML_model", type=str, default="RF")

config = parser.parse_args()

cancer_types = config.cancer_types
run_clustered_mutations = config.clustered_mutations
datasets = sorted(config.datasets)
scATAC_cell_number_filter = config.scATAC_cell_number_filter
waddell_sarc_biph = config.waddell_sarc_biph
waddell_sarc = config.waddell_sarc
waddell_sarc_tsankov_sarc = config.waddell_sarc_tsankov_sarc
waddell_sarc_biph_tsankov_sarc_biph = config.waddell_sarc_biph_tsankov_sarc_biph
annotation_dir = config.annotation_dir
SCLC = config.SCLC
CPTAC = config.CPTAC
combined_CPTAC_ICGC = config.combined_CPTAC_ICGC
tissues_to_consider = config.tissues_to_consider
# tss_filtered = config.tss_filtered
tss_fragment_filter = config.tss_fragment_filter
ML_model = config.ML_model
lung_subtyped = config.lung_subtyped
woo_pcawg = config.woo_pcawg
histologically_subtyped_mutations = config.histologically_subtyped_mutations
de_novo_seurat_clustering = config.de_novo_seurat_clustering
per_donor = config.per_donor
donor_range = config.donor_range
