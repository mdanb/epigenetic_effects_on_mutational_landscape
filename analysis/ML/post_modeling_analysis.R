library(optparse)
source("ML_utils.R")

parser <- OptionParser()
parser <- add_option(parser, c("--datasets"), type="character")
parser <- add_option(parser, c("--cancer_types"), type="character")
parser <- add_option(parser, c("--cell_number_filter"), type="integer")
parser <- add_option(parser, c("--tss_fragment_filter"),
                     type="character", default="-1")
parser <- add_option(parser, c("--tissues_to_consider"), 
                     type="character", default="all")
parser <- add_option(parser, c("--ML_model"), type="character")
parser <- add_option(parser, c("--annotation"), 
                     type="character", default="default_annotation")
parser <- add_option(parser, c("--seed_range"), default="42-42")
parser <- add_option(parser, c("--feature_importance_method"), type="character")
parser <- add_option(parser, c("--folds_for_test_set"), type="character")

args = parse_args(parser)


args = parse_args(parser, args =
                    c("--datasets=Bingren,Greenleaf_brain,Greenleaf_pbmc_bm,Rawlins_fetal_lung,Shendure,Tsankov,Yang_kidney",
                      "--cancer_types=SCLC",
                      "--cell_number_filter=100",
                      "--ML_model=XGB",
                      "--annotation=finalized_annotation",
                      "--seed_range=1-1",
                      "--feature_importance_method=permutation_importance",
                      "--folds_for_test_set=1-10"))

cancer_types = args$cancer_types
cancer_types = unlist(strsplit(cancer_types, split = ","))
datasets = unlist(strsplit(args$datasets, split = ","))
datasets = sort(datasets)
cell_number_filter = args$cell_number_filter
tss_fragment_filter = unlist(strsplit(args$tss_fragment_filter, split = ","))
annotation = args$annotation
tissues_to_consider = strsplit(args$tissues_to_consider,  split=",")
ML_model = args$ML_model
seed_range = unlist(strsplit(args$seed_range, split = "-"))
seed_range = seq(seed_range[1], seed_range[2])
feature_importance_method = args$feature_importance_method

for (cancer_type in cancer_types) {
  # for (tss_filter in tss_fragment_filter) {
  scATAC_sources = construct_sources_string(datasets)
  scATAC_source = paste("scATAC_source", scATAC_sources, "cell_number_filter", 
                        cell_number_filter, sep="_")
  if (tss_fragment_filter != "-1") {
    scATAC_source = paste(scATAC_source, "tss_fragment_filter",
                          tss_fragment_filter, sep="_")
  }
  
  # if (tss_filter != "-1") {
  #   scATAC_source = paste(scATAC_source, "tss_fragment_filter",
  #                         tss_filter, sep="_")
  # }
  
  scATAC_source = paste(scATAC_source, "annotation", annotation, sep="_")
  
  savepath = get_relevant_backwards_elim_dirs(cancer_types=cancer_type, 
                                              # combined_datasets=combined_datasets,
                                              tissues_to_consider=tissues_to_consider,
                                              datasets=datasets,
                                              cell_number_filter=cell_number_filter,
                                              tss_fragment_filter=tss_fragment_filter,
                                              annotation=annotation,
                                              ML_model=ML_model)
}
