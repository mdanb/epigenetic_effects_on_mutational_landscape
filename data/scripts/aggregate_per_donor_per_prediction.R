library(optparse)
library(GenomicRanges)
source("../../analysis/ML/ML_utils.R")

parser <- OptionParser()
parser <- add_option(parser, c("--datasets"), type="character")
parser <- add_option(parser, c("--cancer_type"), type="character")
parser <- add_option(parser, c("--cell_number_filter"), type="integer")
parser <- add_option(parser, c("--tss_fragment_filter"),
                     type="character", default="-1")
parser <- add_option(parser, c("--ML_model"), type="character")
parser <- add_option(parser, c("--annotation"), 
                     type="character", default="default_annotation")
parser <- add_option(parser, c("--seed"), default="1")
parser <- add_option(parser, c("--feature_importance_method"), type="character")
parser <- add_option(parser, c("--fold_for_test_set"), type="character")
parser <- add_option(parser, c("--hundred_kb"), action="store_true", default=F)
parser <- add_option(parser, c("--features_to_aggregate"), type="character")
parser <- add_option(parser, c("--agg_name"), type="character")


# args = parse_args(OptionParser(option_list=option_list))

args = parse_args(parser, args =
                    c("--datasets=Bingren,Shendure",
                      "--cancer_type=Panc-AdenoCA",
                      "--cell_number_filter=100",
                      "--ML_model=XGB",
                      "--annotation=finalized_annotation",
                      "--seed=1",
                      "--feature_importance_method=permutation_importance",
                      "--fold_for_test_set=1",
                      "--features_to_aggregate=esophagus_mucosa Esophageal Epithelial Cell BR,intestine Intestinal epithelial cells SH,pancreas Pancreatic Acinar Cell BR",
                      "--agg_name=epithelial"))

cancer_type = args$cancer_type
datasets = unlist(strsplit(args$datasets, split = ","))
datasets = sort(datasets)
cell_number_filter = args$cell_number_filter
tss_fragment_filter = unlist(strsplit(args$tss_fragment_filter, split = ","))
annotation = args$annotation
ML_model = args$ML_model
seed = args$seed
feature_importance_method = args$feature_importance_method
fold_for_test_set = args$fold_for_test_set
hundred_kb = args$hundred_kb
per_donor = args$per_donor
features_to_aggregate = args$features_to_aggregate
features_to_aggregate = unlist(strsplit(args$features_to_aggregate, split = "\\|"))
agg_name = args$agg_name
agg_name = unlist(strsplit(agg_name, split=","))

dir = get_relevant_backwards_elim_dirs(cancer_types=cancer_type, 
                                       datasets=datasets,
                                       cell_number_filter=cell_number_filter,
                                       tss_fragment_filter=tss_fragment_filter,
                                       annotation=annotation,
                                       ML_model=ML_model,
                                       hundred_kb=hundred_kb,
                                       fold_for_test_set=fold_for_test_set,
                                       seed=seed,
                                       accumulated_seeds=F,
                                       per_donor=F)

dir = paste("../../analysis/ML", dir, sep="/")
file = paste(cancer_type, "per_donor_predictions.csv", sep="_")
fp = paste(dir, file, sep="/")
df = read.csv(fp)
mutations = read.table(paste("../mutation_data", paste0(cancer_type, ".txt"), 
                             sep="/"), header=1)

idx = 1
# chr = read.csv("../processed_data/chr_ranges.csv")[["x"]]
load("../mutation_data/hg19.1Mb.ranges.Polak.Nature2015.RData")
ranges = mutations[, 1:3]
granges = GRanges(seqnames = ranges[["Chr"]], IRanges(start=ranges[["Start"]], 
                 end=ranges[["End"]]))
hits = findOverlaps(granges, interval.ranges)
hits = as.data.frame(hits)
ranges["Chr"] = names(interval.ranges[hits$subjectHits])
ranges = ranges["Chr"]
# sum(a$queryHits == a$subjectHits)
for (features in features_to_aggregate) {
  features = unlist(strsplit(features, split=","))
  mutations_new = data.frame(rep(0, length(ranges[["Chr"]])),
                             row.names = ranges[["Chr"]])
  for (feature in features) {
    donors = df[df[["features"]] == feature, "donor"]
    print(paste0("num donors with ", feature, ": ", length(donors)))
    mutations_new = mutations_new + rowSums(mutations[, donors])
  }
  name = paste(cancer_type, agg_name[idx], sep = "_")
  colnames(mutations_new) = name
  fn = paste0(name, ".csv")
  write.csv(mutations_new, paste("../processed_data", cancer_type, 
                                 fn, sep="/"))
  idx = idx + 1
}

