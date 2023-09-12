library(ggplot2)
library(tidytext)
library(magrittr)
library(RColorBrewer)
library(dplyr)
library(optparse)
library(stringr)
library(tibble)
source("../../utils.R")
source("ML_utils.R")

parser <- OptionParser()
parser <- add_option(parser, c("--datasets"), type="character")
parser <- add_option(parser, c("--cancer_types"), type="character")
parser <- add_option(parser, c("--cell_number_filter"), type="integer")
# parser <- add_option(parser, c("--pie_chart"), action="store_true", default=F)
parser <- add_option(parser, c("--tss_fragment_filter"),
                     type="character", default="-1")
parser <- add_option(parser, c("--ML_model"), type="character")
parser <- add_option(parser, c("--top_features_to_plot"),
                     type="character")
parser <- add_option(parser, c("--top_features_to_plot_feat_imp"),
                     type="character", default=NULL)
parser <- add_option(parser, c("--tissues_to_consider"), 
                     type="character", default="all")
parser <- add_option(parser, c("--annotation"), 
                     type="character", default="default_annotation")
parser <- add_option(parser, c("--seed_range"), default="42-42")
parser <- add_option(parser, c("--robustness_analysis"), action="store_true", 
                     default=F)
parser <- add_option(parser, c("--robustness_accumulated_feature_importance_barplot"), 
                     action="store_true", 
                     default=F)
parser <- add_option(parser, c("--feature_importance_method"), type="character")
parser <- add_option(parser, c("--skip_seeds_robustness"), default="")
parser <- add_option(parser, c("--folds_for_test_set"), type="character")

# args = parse_args(parser, args =
#                     c("--datasets=Bingren,Greenleaf_brain,Greenleaf_colon,Greenleaf_pbmc_bm,Rawlins_fetal_lung,Shendure,Tsankov,Yang_kidney",
#                       "--cancer_types=ColoRect-AdenoCA",
#                       "--cell_number_filter=100",
#                       "--top_features_to_plot=1,2,5,10",
#                       "--ML_model=XGB",
#                       "--annotation=finalized_annotation",
#                       "--robustness_analysis",
#                       "--seed_range=1-10",
#                       "--feature_importance_method=permutation_importance",
#                       "--top_features_to_plot_feat_imp=1,2,5,10",
#                       "--folds_for_test_set=1-10"))

args = parse_args(parser)


ggplot_barplot_helper <- function(df, title, savepath, y, ylab, 
                                  accumulated_imp=F) {
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
  plot = ggplot(df, aes(x=reorder_within(features, -!!sym(y), within=num_features_f,
                                         sep="."), 
                        y=!!sym(y), fill=features)) +
    facet_wrap(~num_features_f, nrow=1, 
               labeller = as_labeller(to), scales = "free") +
    geom_bar(stat="identity", width=1, color="white") +
    xlab("Cell type") +
    ylab(ylab) +
    theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1,
                                     size=15)) +
    guides(fill="none") +
    scale_fill_manual(values=colors) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste(savepath, paste(y, "bar_plot.png", sep="_"), sep="/"), 
         width = 20, height = 15, plot)
}

construct_bar_plots <- function(cancer_type, 
                                combined_datasets,
                                tissues_to_consider,
                                datasets,
                                cell_number_filter,
                                tss_fragment_filter,
                                annotation,
                                ML_model,
                                seed,
                                accumulated_seeds,
                                feature_importance_method, 
                                fold_for_test_set) {
  dirs = get_relevant_backwards_elim_dirs(cancer_types=cancer_type, 
                                          # combined_datasets=combined_datasets,
                                          tissues_to_consider=tissues_to_consider,
                                          datasets=datasets,
                                          cell_number_filter=cell_number_filter,
                                          tss_fragment_filter=tss_fragment_filter,
                                          annotation=annotation,
                                          ML_model=ML_model,
                                          fold_for_test_set=fold_for_test_set,
                                          seed=seed,
                                          accumulated_seeds=accumulated_seeds)
  for (dir in dirs) {
    file = paste(dir, "df_for_feature_importance_plots", sep="/")
    if (feature_importance_method != "default_importance") {
      file = paste(file, feature_importance_method, sep="_")
    }
    file = paste(file, "csv", sep=".") 
    df = as_tibble(read.csv(file))
    title = unlist(strsplit(dir, split ="/"))
    title = title[length(title) - 2]
    df = df[df$num_features %in% top_features_to_plot, ]
    ggplot_barplot_helper(df, title, savepath=dir, 
                          ylab=gsub("_", " ", feature_importance_method), 
                          y=feature_importance_method)
  }
}

construct_boxplots <- function(df, x, y, color, title, savepath, savefile,
                               n_name, facet_var, ylabel, sort_by) {
  from = as.character(unique(df[[facet_var]]))
  to = paste("top", from, "features")
  
  names(to) <- from
  df[[facet_var]] <- factor(df[[facet_var]], levels = unique(df[[facet_var]]))
  
  reorder_within <- function(x, by, within, fun = median, sep = "___") {
    new_x <- paste(x, within, sep = sep)
    ordered_factor <- stats::reorder(new_x, by, FUN = fun)
    levels_reversed <- rev(levels(ordered_factor))
    return(factor(ordered_factor, levels = levels_reversed))  
  }
  
  plot = ggplot(df) +
          geom_boxplot(aes(x = reorder_within(x=!!sym(x), 
                                              by=!!sym(y), 
                                              within=!!sym(facet_var), 
                                              median), y = !!sym(y), color=!!sym(x)),
                       lwd=1.2) +
          geom_text(aes(x = reorder_within(x=!!sym(x), 
                                           by=!!sym(y), 
                                           within=!!sym(facet_var), 
                                           median), y = y_position, label = paste0("n=",
          !!sym(n_name))), vjust = -0.5) +
          facet_wrap(as.formula(paste0("~", facet_var)), nrow=1, 
                     scales = "free_x",
                     labeller = as_labeller(to)) +
          ylab(ylabel) +
          xlab("") +
          ggtitle(title) +
          scale_y_continuous(breaks = seq(round(min(df[[y]]), 2),
                                          round(max(df[[y]]), 2), by = 0.05)) + 
          theme(plot.title = element_text(hjust = 0.5),
                axis.text.x = element_blank())
  print(plot)
  ggsave(paste(savepath, savefile, sep="/"), 
         width = 20, height = 15, plot)
}

construct_df_feature_importances_all_seeds <- function(all_seeds_dirs,
                                                       feature_importance_method,
                                                       top_features_to_plot) {
  df_feature_importances_all_seeds = tibble()
  for (dir in all_seeds_dirs) {
    tryCatch(
      {
        df_feature_importance_path = paste(dir, "backwards_elimination_results", 
                                           "df_for_feature_importance_plots",
                                           sep="/")
        if (feature_importance_method != "default_importance") {
          df_feature_importance_path = paste(df_feature_importance_path,
                                             feature_importance_method,
                                             sep = "_")
        }
        df_feature_importance_path = paste(df_feature_importance_path, 
                                           "csv",
                                           sep=".")
        df_feature_importances = as_tibble(read.csv(df_feature_importance_path))
        df_feature_importances = df_feature_importances %>%
          filter(num_features %in% top_features_to_plot)
        seed = str_extract(dir, pattern="seed_[0-9]*")
        seed = unlist(strsplit(seed, split="_"))
        seed = seed[length(seed)]
        fold_for_test_set = str_extract(dir, pattern="fold_for_test_set_[0-9]*")
        fold_for_test_set = unlist(strsplit(fold_for_test_set, split="_"))
        fold_for_test_set = fold_for_test_set[length(fold_for_test_set)]
        df_feature_importances["seed"] = seed
        df_feature_importances["fold_for_test_set"] = fold_for_test_set

        if (nrow(df_feature_importances_all_seeds) == 0) {
          df_feature_importances_all_seeds = df_feature_importances
        } else {
          df_feature_importances_all_seeds = rbind(df_feature_importances_all_seeds,
                                                   df_feature_importances)
        }
      },
      error = function(e) {
        print(paste("Problem with", dir))
      }
    )
  }
  return(df_feature_importances_all_seeds)
}

plot_accumulated_feature_importance <- function(df_feature_importances_all_seeds,
                                                savepath,
                                                cancer_type) {
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

construct_all_seeds_test_df <- function(top_features_to_plot,
                                        seed_range,
                                        skip_seeds_robustness,
                                        cancer_type,
                                        datasets,
                                        cell_number_filter,
                                        tss_fragment_filter, 
                                        annotation,
                                        tissues_to_consider, 
                                        ML_model,
                                        folds_for_test_set,
                                        feature_importance_method) {
  df = tibble(top_feature = character(0),
              top_n = integer(0),
              test_set_perf = double(0),
              seed = integer(0))
  
  top_features_to_plot = sort(top_features_to_plot, decreasing = T)
  for (seed in seed_range) {
    for (fold in folds_for_test_set) {
      if (!(seed %in% skip_seeds_robustness)) {
        test_dir = construct_backwards_elim_dir(cancer_type,
                                                construct_sources_string(datasets),
                                                cell_number_filter,
                                                tss_fragment_filter, 
                                                annotation,
                                                tissues_to_consider, 
                                                ML_model,
                                                seed,
                                                fold_for_test_set = fold,
                                                test=T)
        
        if (feature_importance_method != "default_importance") {
          model_pattern = paste("^model_iteration_[0-9]+_feature_importance",
                                feature_importance_method, sep="_")
          model_pattern = paste(model_pattern, "pkl", sep=".")
        } else {
          model_pattern = "^model_iteration_[0-9]+\\.pkl"
        }
        total_num_features = list.files(test_dir,
                                        pattern=model_pattern)
        total_num_features = length(total_num_features)
        # sort so that later the first test_file_idx can be checked for 
        # total_num_features == test_file_idx
        test_file_idx = sort(total_num_features - top_features_to_plot + 1)
        test_perf_filenames = paste("model_iteration", test_file_idx,
                                    sep = "_")
        if (feature_importance_method != "default_importance") {
          test_perf_filenames = paste(test_perf_filenames,
                                      "feature_importance", 
                                      feature_importance_method, 
                                      sep = "_")
        }
        test_perf_filenames = paste(test_perf_filenames, 
                                    "test_performance.txt", sep="_")
        test_set_perf_files = paste(test_dir, test_perf_filenames, sep="/")
        
        idx = 1
        for (file in test_set_perf_files) {
          if (!(idx == length(test_set_perf_files) &&
                test_file_idx[length(test_file_idx)] == total_num_features)) {
            top_feature_file = paste0("top_features_iteration_",
                                      test_file_idx[idx])
            
            if (feature_importance_method != "default_importance") {
              top_feature_file = paste(top_feature_file,
                                       "by", 
                                       feature_importance_method, 
                                       sep = "_")
            }
            
            top_feature_file = paste0(top_feature_file, ".txt")
            top_feature_fp = paste(test_dir, top_feature_file, sep="/")
            top_feature = readLines(top_feature_fp, n = 1)
            top_feature = substring(top_feature, 4, nchar(top_feature))
          }
          suppressWarnings({
            perf = read.table(file, header=F)[1,1]
          })
          df = df %>% add_row(top_feature = top_feature,
                              top_n = top_features_to_plot[idx],
                              test_set_perf = perf,
                              seed = seed)
          idx = idx + 1
        }
      }
    }
  }
  return(df)
}

cancer_types = args$cancer_types
cancer_types = unlist(strsplit(cancer_types, split = ","))
top_features_to_plot = args$top_features_to_plot
top_features_to_plot = as.integer(unlist(strsplit(top_features_to_plot, 
                                                  split = ",")))
datasets = unlist(strsplit(args$datasets, split = ","))
datasets = sort(datasets)
cell_number_filter = args$cell_number_filter
tss_fragment_filter = unlist(strsplit(args$tss_fragment_filter, split = ","))
annotation = args$annotation
tissues_to_consider = strsplit(args$tissues_to_consider,  split=",")
ML_model = args$ML_model
seed_range = unlist(strsplit(args$seed_range, split = "-"))
seed_range = seq(seed_range[1], seed_range[2])
robustness_analysis = args$robustness_analysis
robustness_accumulated_feature_importance_barplot = 
  args$robustness_accumulated_feature_importance_barplot
feature_importance_method = args$feature_importance_method
if (!is.null(args$top_features_to_plot_feat_imp)) {
  top_features_to_plot_feat_imp = as.numeric(unlist(strsplit(
                                              args$top_features_to_plot_feat_imp, 
                                              split=",")))
}
skip_seeds_robustness = args$skip_seeds_robustness
if (!is.null(skip_seeds_robustness)) {
  skip_seeds_robustness = unlist(strsplit(args$skip_seeds_robustness, split=","))
}

folds_for_test_set = args$folds_for_test_set
folds_for_test_set = unlist(strsplit(args$folds_for_test_set, split = "-"))
folds_for_test_set = seq(folds_for_test_set[1], folds_for_test_set[2])

if (!robustness_analysis) {
  tissues_to_consider = paste(unlist(tissues_to_consider, "_"))

  for (seed in seed_range) {
    for (fold in folds_for_test_set) {
      for (cancer_type in cancer_types) {
        construct_bar_plots(cancer_type, 
                            combined_datasets,
                            tissues_to_consider,
                            datasets,
                            cell_number_filter,
                            tss_fragment_filter,
                            annotation,
                            ML_model,
                            seed,
                            accumulated_seeds=F,
                            feature_importance_method=feature_importance_method,
                            fold_for_test_set=fold)
      }
    }
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
    
    scATAC_source = paste(scATAC_source, "annotation", annotation, sep="_")

    savepath = get_relevant_backwards_elim_dirs(cancer_types=cancer_type, 
                                                # combined_datasets=combined_datasets,
                                                tissues_to_consider=tissues_to_consider,
                                                datasets=datasets,
                                                cell_number_filter=cell_number_filter,
                                                tss_fragment_filter=tss_fragment_filter,
                                                annotation=annotation,
                                                ML_model=ML_model,
                                                accumulated_seeds=T)
    
    dirs = list.dirs(paste("../../figures", "models", ML_model, cancer_type, 
                           sep="/"), recursive = F)
    
    combos = expand.grid(seed = seed_range, fold = folds_for_test_set)
    seed_fold_for_test_combinations = apply(combos, 1, function(row) {
      paste(scATAC_source, "seed", row["seed"], "fold_for_test_set", row["fold"], 
            sep="_")
    })
  
    all_seeds_dirs = dirs[basename(dirs) %in% seed_fold_for_test_combinations]
    
    df_feature_importances_all_seeds = 
          construct_df_feature_importances_all_seeds(all_seeds_dirs,
                                                     feature_importance_method,
                                                     top_features_to_plot)

    if (robustness_accumulated_feature_importance_barplot) {
      plot_accumulated_feature_importance(df_feature_importances_all_seeds,
                                          savepath,
                                          cancer_type)
    }
    
    df = construct_all_seeds_test_df(top_features_to_plot=top_features_to_plot,
                                     seed_range=seed_range,
                                     skip_seeds_robustness=skip_seeds_robustness,
                                     cancer_type=cancer_type,
                                     datasets=datasets,
                                     cell_number_filter=cell_number_filter,
                                     tss_fragment_filter=tss_fragment_filter, 
                                     annotation=annotation,
                                     tissues_to_consider=tissues_to_consider, 
                                     ML_model=ML_model,
                                     folds_for_test_set=folds_for_test_set,
                                     feature_importance_method=feature_importance_method)
    
    # y_position is for plotting number of times feature appears at the top of
    # the boxplot. 
    df_feat_imp = df_feature_importances_all_seeds %>% 
          group_by(num_features, seed, fold_for_test_set) %>%
          group_by(num_features, features) %>%
          mutate(n_feature = n(), y_position = max(permutation_importance)) %>%
          filter(num_features %in% top_features_to_plot_feat_imp)
      
    construct_boxplots(df=df_feat_imp, x="features", y="permutation_importance", 
                       color="features", title=cancer_type, savepath=savepath,
                       savefile="feature_importance.png", 
                       n_name="n_feature", facet_var="num_features",
                       ylabel="Permutation Importance")
    
    df_feat_imp_at_least_50_n = df_feat_imp %>% 
      filter(n_feature >= 50)
    
    construct_boxplots(df=df_feat_imp_at_least_50_n, x="features", 
                       y="permutation_importance", 
                       color="features", title=cancer_type, savepath=savepath,
                       savefile="feature_importance_50_n_min.png", 
                       n_name="n_feature", facet_var="num_features",
                       ylabel="Permutation Importance")
    
    df_feat_imp = df_feature_importances_all_seeds %>% 
      group_by(num_features, seed, fold_for_test_set) %>%
      mutate(max_feature_importance=max(permutation_importance)) %>%
      filter(permutation_importance == max_feature_importance) %>%
      ungroup() %>%
      group_by(num_features, features) %>%
      mutate(n_feature = n(), y_position = max(permutation_importance)) %>%
      filter(num_features %in% top_features_to_plot_feat_imp)
    
    construct_boxplots(df=df_feat_imp, x="features", 
                       y="permutation_importance", 
                       color="features", title=cancer_type, savepath=savepath,
                       savefile="feature_importance_top_feat_only.png", 
                       n_name="n_feature", facet_var="num_features",
                       ylabel="Permutation Importance")
    
    df_test = df %>% 
      group_by(top_n, top_feature) %>%
      mutate(n_top_feature = n(), y_position = max(test_set_perf))
    construct_boxplots(df=df_test, x="top_feature", y="test_set_perf", 
                       color="top_feature", title=cancer_type, savepath=savepath,
                       savefile="test_set_boxplots.png",
                       n_name="n_top_feature", facet_var="top_n",
                       ylabel="Test set R^2",
                       sort_by="n")
    
    df_val = df_feature_importances_all_seeds %>% 
            group_by(num_features, seed, fold_for_test_set) %>%
            mutate(max_feature_importance=max(permutation_importance)) %>%
            filter(permutation_importance == max_feature_importance) %>%
            ungroup() %>%
            group_by(num_features, features) %>%
            mutate(n_feature = n(), y_position = max(score))

    construct_boxplots(df_val, x="features", y="score", 
                       color="features", title=cancer_type, savepath=savepath,
                       savefile="validation_boxplots.png", 
                       n_name="n_feature", facet_var="num_features",
                       ylabel="Validation R^2", 
                       sort_by="score")
    }
}
