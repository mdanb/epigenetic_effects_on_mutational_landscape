library(optparse)
library(gridExtra)
library(tidyverse)
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
parser <- add_option(parser, c("--feature_importance_method"), type="character")
parser <- add_option(parser, c("--folds_for_test_set"), type="character")
parser <- add_option(parser, c("--plot_bin_counts"), action="store_true", 
                     default=F)
parser <- add_option(parser, c("--plot_bins_volcano"), action="store_true", 
                     default=F)
parser <- add_option(parser, c("--hundred_kb"), action="store_true", 
                     default=F)

args = parse_args(parser)


args = parse_args(parser, args =
                    c("--datasets=Bingren,Greenleaf_brain,Greenleaf_pbmc_bm,Rawlins_fetal_lung,Shendure,Tsankov,Yang_kidney",
                      "--cancer_types=Lung-SCC",
                      "--cell_number_filter=100",
                      "--ML_model=XGB",
                      "--annotation=finalized_annotation",
                      "--feature_importance_method=permutation_importance",
                      "--folds_for_test_set=1-10",
                      "--hundred_kb"))

cancer_types = args$cancer_types
cancer_types = unlist(strsplit(cancer_types, split = ","))
datasets = unlist(strsplit(args$datasets, split = ","))
datasets = sort(datasets)
cell_number_filter = args$cell_number_filter
tss_fragment_filter = unlist(strsplit(args$tss_fragment_filter, split = ","))
annotation = args$annotation
tissues_to_consider = strsplit(args$tissues_to_consider,  split=",")
ML_model = args$ML_model
feature_importance_method = args$feature_importance_method

folds_for_test_set = args$folds_for_test_set
folds_for_test_set = unlist(strsplit(args$folds_for_test_set, split = "-"))
folds_for_test_set = seq(folds_for_test_set[1], folds_for_test_set[2])
plot_bin_counts = args$plot_bin_counts
hundred_kb = args$hundred_kb

plot_dist <- function(X, bin_names) {
  df <- data.frame(index = 1:nrow(X), 
                          value = X)
  use_names = seq(1, length(df$index), by=10)
  use_names = append(use_names, nrow(X))
  p = ggplot(df) +
    geom_bar(aes(x=index, y=value),
             stat="identity") +
    xlab("Bin") +
    ylab("Counts") +
    geom_text(x=nrow(df) / 2, y=max(df[["value"]]), 
              aes(label=paste0("n=", sum(df[["value"]]))),
              size=10) +
    scale_x_continuous(name="Row Names", breaks=df$index[use_names], 
                       labels=bin_names[use_names]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  return(p)
}

scatac_counts_plots <- list()
mut_counts_plots <- list()

for (cancer_type in cancer_types) {
  scATAC_sources = construct_sources_string(datasets)
  if (plot_bins_volcano) {
    dir = construct_dir(scATAC_sources,
                        cell_number_filter,
                        tss_fragment_filter,
                        annotation,
                        "all_seeds",
                        -1,
                        ML_model,
                        cancer_type)
    fp = paste("../../figures", dir, "errors_df.csv", sep="/")
    
    if (hundred_kb) {
      load("../../data/100kb_interval_ranges.Rdata")
    }
    else {
      load("../../data/mutation_data/hg19.1Mb.ranges.Polak.Nature2015.RData")
    }
    
    df = read.csv(fp)
    colnames(df)[1] = "names"
    chr_to_range = data.frame(interval.ranges@ranges)
    df = left_join(df, chr_to_range)
    df["neglog_q"] = -log10(df["q.value"])
    df = df %>% 
          filter(rejected == "True")
    ggplot(df) + 
      geom_point(aes(x=normalized_percent_error, y=neglog_q)) +
      ylab("Mutation Enrichment, -log10 FDR")
    
    df[, c("names", "start", "end", "neglog_q")]
  }
  
  
  if (plot_bin_counts) {
    for (fold in folds_for_test_set) {
    # for (tss_filter in tss_fragment_filter) {
      dir = construct_dir(scATAC_sources,
                    cell_number_filter,
                    tss_fragment_filter,
                    annotation,
                    1,
                    fold,
                    ML_model,
                    cancer_type)
      X_test = read.csv(paste(dir, "X_test.csv", sep="/"), row.names = 1,
                        check.names = FALSE)
      y_test = read.csv(paste(dir, "y_test.csv", sep="/"), row.names = 1)
      backwards_elim_dir = construct_backwards_elim_dir(cancer_type, 
                                                        scATAC_sources, 
                                                        cell_number_filter,
                                                        tss_fragment_filter, 
                                                        annotation,
                                                        tissues_to_consider, 
                                                        ML_model,
                                                        1,
                                                        fold_for_test_set=fold,
                                                        test=T)
      fp = paste(backwards_elim_dir, 
              "top_features_iteration_19_by_permutation_importance.txt", sep="/")
      con = file(fp, "r")
      top_feature=readLines(con, n = 1)
      top_feature=substr(top_feature, 4, nchar(top_feature))
      X_test_top_feature=data.frame(X_test[, top_feature])
      colnames(X_test_top_feature) = "value"
      p1 = plot_dist(X_test_top_feature, rownames(X_test))
      scatac_counts_plots[[fold]] = p1
      colnames(y_test) = "value"
      p2 = plot_dist(y_test, rownames(X_test))
      mut_counts_plots[[fold]] = p2
    }
    pdf("../../figures/per_fold_scatac_counts.pdf", width = 50, height = 20)
    do.call(grid.arrange, c(scatac_counts_plots, nrow=2)) 
    dev.off()
    pdf("../../figures/per_fold_mut_counts.pdf", width = 50, height = 20)
    do.call(grid.arrange, c(mut_counts_plots, nrow=2))
    dev.off()
  }
}

