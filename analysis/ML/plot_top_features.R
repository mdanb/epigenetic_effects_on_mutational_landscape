library(ggplot2)
library(tidytext)
library(magrittr)
library(RColorBrewer)
library(dplyr)
library(optparse)
library(stringr)
library(tibble)
source("../../utils.R")

parser <- OptionParser()
parser <- add_option(parser, c("--datasets"), type="character")
parser <- add_option(parser, c("--cancer_types"), type="character")
parser <- add_option(parser, c("--cell_number_filter"), type="integer")
parser <- add_option(parser, c("--pie_chart"), action="store_true", default=F)
parser <- add_option(parser, c("--tss_fragment_filter"),
                     type="character", default="-1")
parser <- add_option(parser, c("--ML_model"), type="character")
# parser <- add_option(parser, c("--num_iter_skips"), 
#                      type="integer", default=5)
parser <- add_option(parser, c("--top_features_to_plot"),
                     type="character")
parser <- add_option(parser, c("--tissues_to_consider"), 
                     type="character", default="all")
parser <- add_option(parser, c("--annotation"), 
                     type="character", default="default_annotation")
# parser <- add_option(parser, c("--iters_dont_skip"), default="18")
parser <- add_option(parser, c("--seed"), default="42")
parser <- add_option(parser, c("--robustness_analysis"), action="store_true", 
                     default=F)
parser <- add_option(parser, c("--robustness_feature_importance_barplot"), 
                     action="store_true", 
                     default=F)
parser <- add_option(parser, c("--robustness_test_perf_boxplot"), 
                     action="store_true", 
                     default=F)
parser <- add_option(parser, c("--robustness_seed_range"), type="character",
                     default="1-100")
parser <- add_option(parser, c("--SCLC"), action="store_true", default=F)
parser <- add_option(parser, c("--lung_subtyped"), action="store_true", default=F)
parser <- add_option(parser, c("--woo_pcawg"), action="store_true", default=F)
parser <- add_option(parser, c("--histologically_subtyped_mutations"), 
                     action="store_true", default=F)
parser <- add_option(parser, c("--de_novo_seurat_clustering"), 
                     action="store_true", default=F)
parser <- add_option(parser, c("--per_donor"), action="store_true", default=F)
parser <- add_option(parser, c("--CPTAC"), action="store_true", default=F)
parser <- add_option(parser, c("--meso"), action="store_true", default=F)
parser <- add_option(parser, c("--combined_CPTAC_ICGC"), action="store_true", default=F)
parser <- add_option(parser, c("--donor_range"), type="character", default=NULL)
parser <- add_option(parser, c("--RNA_subtyped"), action="store_true", default=F)

# args = parse_args(parser, args = c("--cancer_types=Lung-AdenoCA",
#                                    "--robustness_analysis",
#                                    "--datasets=Tsankov",
#                                    "--cell_number_filter=30",
#                                    "--ML_model=XGB",
#                                    "--annotation=finalized_annotation",
#                                    "--iters_dont_skip=17",
#                                    "--robustness_top_ns=2,4"))
# args = parse_args(parser, args =
#                     c("--datasets=Tsankov",
#                       "--cancer_types=Lung-SCC",
#                       "--cell_number_filter=30",
#                       "--ML_model=XGB",
#                       "--annotation=finalized_annotation",
#                       "--robustness_analysis",
#                       "--robustness_seed_range=1-100",
#                       "--top_features_to_plot=15,10,5,2",
#                       "--robustness_test_perf_boxplot"))

args = parse_args(parser)

construct_backwards_elim_dir <- function(cancer_type, scATAC_source, 
                                         cell_number_filter,
                                         tss_fragment_filter, 
                                         annotation,
                                         tissues_to_consider, 
                                         ML_model,
                                         seed,
                                         test=F) {
  scATAC_source = paste("scATAC_source", scATAC_source, "cell_number_filter", 
                        cell_number_filter, sep="_")
  
  if (tss_fragment_filter != -1) {
    scATAC_source = paste(scATAC_source, "tss_fragment_filter", 
                          tss_fragment_filter, sep="_")
  }

  scATAC_source = paste(scATAC_source, "annotation", annotation, "seed", seed, 
                        sep="_")
  dir = paste("models", ML_model, cancer_type, scATAC_source,
              "backwards_elimination_results", sep="/")
  
  if (tissues_to_consider != "all") {
    dir = paste(dir, tissues_to_consider, sep="_")
  }
  
  if (!test) {
    dir = paste("../../figures", dir, sep="/")
  }
  
  return(dir)
}

construct_sources_string <- function(datasets) {
  scATAC_sources = ""
  
  for (i in 1:length(datasets)) {
    dataset = datasets[i]
    if (scATAC_sources == "") {
      scATAC_sources = dataset
    }
    else {
      scATAC_sources = paste(scATAC_sources, dataset, sep="_")
    }
  }
  return(scATAC_sources)
}

get_relevant_backwards_elim_dirs <- function(args, accumulated_seeds=F) {
    cancer_types = args$cancer_types
    cancer_types = unlist(strsplit(cancer_types, split = ","))
    combined_datasets = args$combined_datasets
    tissues_to_consider = paste(unlist(strsplit(args$tissues_to_consider, 
                                                split=","), 
                                       "_"))
    datasets = unlist(strsplit(args$datasets, split = ","))
    cell_number_filter = args$cell_number_filter
    tss_fragment_filter = unlist(strsplit(args$tss_fragment_filter, split = ","))
    annotation = args$annotation
    ML_model = args$ML_model
    if (!accumulated_seeds) {
      seed = args$seed
    } else {
      seed = "all_seeds"
    }
    scATAC_sources = construct_sources_string(datasets)

    backward_elim_dirs = list()
    for (cancer_type in cancer_types) {
      for (tss_filter in tss_fragment_filter) {
        backward_elim_dirs = append(backward_elim_dirs,
                                    construct_backwards_elim_dir(cancer_type, 
                                                                 scATAC_sources,
                                                                 cell_number_filter,
                                                                 tss_filter,
                                                                 annotation,
                                                                 tissues_to_consider,
                                                                 ML_model,
                                                                 seed))
      }
    }
    return(unlist(backward_elim_dirs))
}

# get_n_colors <- function(n) {
#   qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#   col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, 
#                              rownames(qual_col_pals)))
#   set.seed(1)
#   cols <- sample(col_vector, n)
#   return(cols)
# }

construct_pie_charts <- function(args) {
  dirs = get_relevant_backwards_elim_dirs(args)
  for (dir in dirs) {
    file = paste(dir, "df_for_feature_importance_plots.csv", sep="/")
    df = read.csv(file)
    df$num_features_f = factor(df$num_features, levels=unique(df$num_features))
    colors = get_n_colors(20, 1)
    from = as.character(unique(df$num_features))
    to = paste(paste(levels(df$num_features_f), "features"),
               paste("(R^2=", as.character(round(unique(df$score*100), 1)), 
                     ")", sep=""), sep=" ")
    
    names(to) <- from
    # for (num_feats in unique(df$num_features)) {
    #   score = df[df$num_features == num_feats, ]$score[1]
    #   labels = append(labels, num_feats = score)
    # }
    title = unlist(strsplit(dir, split ="/"))
    title = title[length(title)]
    print("title: ")
    print(title)
    plot = ggplot(df, aes(x="", y=importance, fill=features)) +
                  facet_wrap(~num_features_f, nrow=1, 
                             labeller = as_labeller(to)) +
                  geom_bar(stat="identity", width=1, color="white") +
                  coord_polar("y", start=0) +
                  xlab("") +
                  ylab("") +
                  theme(axis.text = element_blank(),
                        axis.ticks = element_blank(),
                        panel.grid  = element_blank()) +
                  scale_fill_manual(values=colors) +
                  ggtitle() +
                  # geom_text(aes(y = ypos, label = features), color = "white", size=1) +
                  theme(plot.title = element_text(hjust = 0.5))
    ggsave(paste(dir, "pie_chart.png", sep="/"), width = 20,
           height = 6, plot)
  }
}

ggplot_barplot_helper <- function(df, title, savepath, accumulated_imp=F) {
  df$num_features_f = factor(df$num_features, levels=unique(df$num_features))
  colors = get_n_colors(20, 1)
  
  from = as.character(unique(df$num_features))
  to = paste(levels(df$num_features_f), "features", sep=" ")
  if (!accumulated_imp) {
    to = paste(to, paste("(R^2=", 
                         as.character(round(unique(df$score*100), 1)), ")",
                         sep=""))
  }
  names(to) <- from
  plot = ggplot(df, aes(x=reorder_within(features, -importance, within=num_features_f,
                                         sep="."), 
                        y=importance, fill=features)) +
    facet_wrap(~num_features_f, nrow=1, 
               labeller = as_labeller(to), scales = "free") +
    geom_bar(stat="identity", width=1, color="white") +
    xlab("Cell type") +
    {if (!accumulated_imp) 
      ylab("Percent importance (%)")
    else {
      ylab("Accumulated percent importance")
    }} +
    # theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1),
    #       aspect.ratio = 1.1/1) +
    theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1,
                                     size=15)) +
    guides(fill="none") +
    scale_fill_manual(values=colors) +
    # ggtitle(paste0(unlist(strsplit(dir, split ="/"))[3], " (R^2=",
    #               as.character(round(unique(df$score*100), 1)), ")")) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste(savepath, "bar_plot.png", sep="/"), width = 20, height = 15, plot)
}

construct_bar_plots <- function(args) {
  dirs = get_relevant_backwards_elim_dirs(args)
  for (dir in dirs) {
    file = paste(dir, "df_for_feature_importance_plots.csv", sep="/")
    df = as_tibble(read.csv(file))
    title = unlist(strsplit(dir, split ="/"))
    title = title[length(title) - 2]
    ggplot_barplot_helper(df, title, savepath=dir)
  }
}

construct_test_boxplots <- function(df, title, savepath) {
  from = as.character(unique(df$top_n))
  to = paste("top", from, "features")
  
  names(to) <- from
  
  plot = ggplot(df) +
          geom_boxplot(aes(x = top_feature, y = test_set_perf, color=top_feature),
                       lwd=1.2) +
          geom_text(aes(x = top_feature, y = y_position, label = paste0("n=", n_top_feature)),
                    vjust = -0.5) +
          facet_wrap(~top_n, nrow=1, 
                     scales = "fixed",
                     labeller = as_labeller(to)) +
          xlab("Cell type") +
          ylab("Test set R^2") +
          ggtitle(title) +
          theme(plot.title = element_text(hjust = 0.5),
                axis.text.x = element_blank())
  ggsave(paste(savepath, "test_set_boxplots.png", sep="/"), 
         width = 20, height = 15, plot)
}

cancer_types = args$cancer_types
cancer_types = unlist(strsplit(cancer_types, split = ","))
top_features_to_plot = args$top_features_to_plot
print(top_features_to_plot)
top_features_to_plot = as.integer(unlist(strsplit(top_features_to_plot, 
                                                  split = ",")))
print(top_features_to_plot)
# iters_dont_skip = args$iters_dont_skip
# iters_dont_skip = as.integer(unlist(strsplit(iters_dont_skip, split = ",")))
datasets = unlist(strsplit(args$datasets, split = ","))
cell_number_filter = args$cell_number_filter
tss_fragment_filter = unlist(strsplit(args$tss_fragment_filter, split = ","))
# num_iter_skips = args$num_iter_skips
annotation = args$annotation
tissues_to_consider = strsplit(args$tissues_to_consider,  split=",")
ML_model = args$ML_model
seed = args$seed
robustness_analysis = args$robustness_analysis
robustness_seed_range = as.integer(unlist(strsplit(args$robustness_seed_range, split="-")))
robustness_test_perf_boxplot = args$robustness_test_perf_boxplot
robustness_feature_importance_barplot = args$robustness_feature_importance_barplot
cancer_types = paste(cancer_types, collapse = " ")
# print(top_features_to_plot)
# print(datasets)
# print(cell_number_filter)
# print(tss_fragment_filter)
# print(annotation)
# print(tissues_to_consider)
# print(ML_model)
# print(seed)
# print(robustness_analysis)
# print(robustness_seed_range)
# print(robustness_test_perf_boxplot)
# print(robustness_feature_importance_barplot)
SCLC = args$SCLC
lung_subtyped = args$lung_subtyped
woo_pcawg = args$woo_pcawg
histologically_subtyped_mutations = args$histologically_subtyped_mutations
de_novo_seurat_clustering = args$de_novo_seurat_clustering
per_donor = args$per_donor
CPTAC = args$CPTAC
meso = args$meso
combined_CPTAC_ICGC = args$combined_CPTAC_ICGC
donor_range = args$donor_range
RNA_subtyped = args$RNA_subtyped

print(robustness_analysis)
if (!robustness_analysis) {
  prep_dfs_command = paste("python3 ../../data/scripts/prep_dfs_for_feature_importance_plots.py", 
                           "--datasets", 
                           paste(datasets, collapse=" "), 
                           "--annotation", annotation, 
                           "--tissues_to_consider",paste(tissues_to_consider, 
                                                         collapse = " "), 
                           "--cell_number_filter", cell_number_filter, 
                           "--top_features_to_plot", top_features_to_plot,
                           "--tss_fragment_filter",  paste(tss_fragment_filter, 
                                                           collapse = " "),
                           "--ML_model", ML_model,
                           "--cancer_types", cancer_types, 
                           "--seed", seed)
  if (SCLC) {
    prep_dfs_command = paste(prep_dfs_command, "--SCLC")
  } else if (woo_pcawg) {
    prep_dfs_command = paste(prep_dfs_command, "--woo_pcawg")
  } else if (histologically_subtyped_mutations) {
    prep_dfs_command = paste(prep_dfs_command, "--histologically_subtyped_mutations")
  } else if (de_novo_seurat_clustering) {
    prep_dfs_command = paste(prep_dfs_command, "--de_novo_seurat_clustering")
  } else if (per_donor) {
    prep_dfs_command = paste(prep_dfs_command, "--per_donor")
  } else if (CPTAC) {
    prep_dfs_command = paste(prep_dfs_command, "--CPTAC")
  } else if (meso) {
    prep_dfs_command = paste(prep_dfs_command, "--meso")
  } else if (combined_CPTAC_ICGC) {
    prep_dfs_command = paste(prep_dfs_command, "--combined_CPTAC_ICGC")
  } else if (!is.null(donor_range)) {
    prep_dfs_command = paste(prep_dfs_command, "--donor_range")
  } else if (RNA_subtyped) {
    prep_dfs_command = paste(prep_dfs_command, "--RNA_subtyped")
  } else if (lung_subtyped) {
    prep_dfs_command = paste(prep_dfs_command, "--lung_subtyped")
  }
  
  print(paste0("Prepping feature importance dfs for seed ", seed, "..."))
  print(prep_dfs_command)
  system(prep_dfs_command)
  print(paste0("Done prepping feature importance dfs for seed ", seed, "!"))
  
  tissues_to_consider = paste(unlist(tissues_to_consider, "_"))
  
  if (args$pie_chart) {
    construct_pie_charts(args)
  } else {
    construct_bar_plots(args)
  }
} else {
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
    
    scATAC_source = paste(scATAC_source, "annotation", annotation, 
                          sep="_")
    savepath = get_relevant_backwards_elim_dirs(args, accumulated_seeds = T)
    if (robustness_feature_importance_barplot) {
      dirs = list.dirs(paste("../../figures", "models", ML_model, cancer_type, 
                             sep="/"), recursive = F)
      all_seeds_dirs = dirs[grepl(scATAC_source, dirs)]
      all_seeds_dirs = all_seeds_dirs[!grepl("all_seeds", all_seeds_dirs)]
      df_feature_importances_all_seeds = tibble()
      for (dir in all_seeds_dirs) {
        df_feature_importances = as_tibble(read.csv(paste(dir, "backwards_elimination_results", 
                                       "df_for_feature_importance_plots.csv",
                                       sep="/")))
        df_feature_importances = df_feature_importances %>%
                                  filter(num_features %in% top_features_to_plot)
        temp = unlist(strsplit(dir, split="_"))
        seed = temp[length(temp)]
        df_feature_importances["seed"] = seed
        if (nrow(df_feature_importances_all_seeds) == 0) {
          df_feature_importances_all_seeds = df_feature_importances
        } else {
          df_feature_importances_all_seeds = rbind(df_feature_importances_all_seeds,
                                                   df_feature_importances)
        }
      }
      # }
      df_accumulated_imp = df_feature_importances_all_seeds %>% 
                              group_by(features, num_features) %>%
                              summarise(sum(importance))
      colnames(df_accumulated_imp)[3] = "importance"
      dir.create(savepath, recursive = T)
      ggplot_barplot_helper(df=df_accumulated_imp, 
                            title=cancer_type, 
                            savepath=savepath,
                            accumulated_imp=T)
    }
    
    if (robustness_test_perf_boxplot) {
      df = tibble(top_feature = character(0),
                  top_n = integer(0),
                  test_set_perf = double(0),
                  seed = integer(0))
      
      for (seed in seq(robustness_seed_range[1], robustness_seed_range[2])) {
          test_dir = construct_backwards_elim_dir(cancer_type,
                                                  construct_sources_string(datasets),
                                                  cell_number_filter,
                                                  tss_fragment_filter, 
                                                  annotation,
                                                  tissues_to_consider, 
                                                  ML_model,
                                                  seed,
                                                  test=T)
          total_num_features = length(list.files(test_dir,
                                        pattern="model_iteration_[0-9]*\\.pkl")) + 1
          test_file_idx = total_num_features - top_features_to_plot + 1
          test_perf_filenames = paste("model_iteration", test_file_idx, "test_performance.txt",
                          sep = "_")
          test_set_perf_files = paste(test_dir, test_perf_filenames, sep="/")
          
          idx = 1
          for (file in test_set_perf_files) {
            tryCatch(
                {
                  top_feature_file = paste0("top_features_iteration_",
                                          test_file_idx[idx],
                                          ".txt")
                  top_feature_fp = paste(test_dir, top_feature_file, sep="/")
                  top_feature = readLines(top_feature_fp, n = 1)
                  top_feature = substring(top_feature, 4, nchar(top_feature))
                  suppressWarnings({
                    perf = read.table(file)[1,1]
                  })
                  df = df %>% add_row(top_feature = top_feature,
                                      top_n = top_features_to_plot[idx],
                                      test_set_perf = perf,
                                      seed = seed)
                  idx = idx + 1
                },
                error = function(e) {
                }
            )
          }
      }
      df = df %>% 
            group_by(top_n, top_feature) %>%
            mutate(n_top_feature = n(), y_position = max(test_set_perf))
      construct_test_boxplots(df, title=cancer_type, savepath=savepath)
      }
    }
}
