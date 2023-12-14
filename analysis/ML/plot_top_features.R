library(ggplot2)
library(tidytext)
library(magrittr)
library(RColorBrewer)
library(dplyr)
library(optparse)
library(stringr)
library(tibble)
library(gridExtra)
library(patchwork)
library(hash)

source("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/utils.R")
source("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML/ML_utils.R")
# 
# source("/home/mdanb/research/mount_sinai/epigenetic_effects_on_mutational_landscape/utils.R")
# source("/home/mdanb/research/mount_sinai/epigenetic_effects_on_mutational_landscape/analysis/ML/ML_utils.R")

parser <- OptionParser()
parser <- add_option(parser, c("--datasets"), type="character",
                     default=NULL)
parser <- add_option(parser, c("--cancer_types"), type="character")
parser <- add_option(parser, c("--cell_number_filter"), type="integer")
# parser <- add_option(parser, c("--pie_chart"), action="store_true", default=F)
parser <- add_option(parser, c("--tss_fragment_filter"),
                     type="character", default="-1")
parser <- add_option(parser, c("--ML_model"), type="character", default="XGB")
parser <- add_option(parser, c("--top_features_to_plot"),
                     type="character", default="1,2,5,10")
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
parser <- add_option(parser, c("--folds_for_test_set"), type="character", 
                     default="1-10")
parser <- add_option(parser, c("--feat_imp_min_n_robustness"), type="integer")
parser <- add_option(parser, c("--plot_fold_on_test_set_plot"), action="store_true", 
                     default=F)
parser <- add_option(parser, c("--hundred_kb"), action="store_true", default=F)
parser <- add_option(parser, c("--per_donor"), action="store_true", default=F)
parser <- add_option(parser, c("--grid_analysis"), action="store_true", default=F)
parser <- add_option(parser, c("--grid_cell_types"), type="character")

# args = parse_args(parser, args =
#                       c("--datasets=Bingren_adult_brain",
#                         "--cancer_types=CNS-Oligo",
#                         "--cell_number_filter=100",
#                         "--ML_model=XGB",
#                         "--annotation=finalized_annotation",
#                         "--robustness_analysis",
#                         "--seed_range=1-10",
#                         "--feature_importance_method=permutation_importance",
#                         "--top_features_to_plot_feat_imp=10",
#                         "--folds_for_test_set=1-10",
#                         "--feat_imp_min_n_robustness=50"))

# args = parse_args(parser, args =
#                     c("--datasets=Bingren,Shendure",
#                       "--cancer_types=Panc-AdenoCA",
#                       "--cell_number_filter=100",
#                       "--top_features_to_plot=1",
#                       "--ML_model=XGB",
#                       "--annotation=finalized_annotation",
#                       "--seed_range=1-1",
#                       "--feature_importance_method=permutation_importance",
#                       "--top_features_to_plot_feat_imp=1",
#                       "--folds_for_test_set=1-1",
#                       "--per_donor"))

# args = parse_args(parser, args =
#                     c("--datasets=Greenleaf_colon",
#                       "--cancer_types=ColoRect-AdenoCA",
#                       "--cell_number_filter=1",
#                       "--top_features_to_plot=1",
#                       "--ML_model=XGB",
#                       "--annotation=Greenleaf_colon_remove_cancer_polyp_merge_normal_unaffected",
#                       "--seed_range=1-10",
#                       "--feature_importance_method=permutation_importance",
#                       "--top_features_to_plot_feat_imp=1",
#                       "--folds_for_test_set=1-10",
#                       "--robustness_analysis",
#                       "--hundred_kb"))

# args = parse_args(parser, args =
#                     c("--datasets=Tsankov,Rawlins_fetal_lung",
#                       "--cancer_types=SCLC",
#                       "--cell_number_filter=1",
#                       "--top_features_to_plot=1",
#                       "--ML_model=XGB",
#                       "--annotation=finalized_annotation",
#                       "--seed_range=1-1",
#                       "--feature_importance_method=permutation_importance",
#                       "--top_features_to_plot_feat_imp=1",
#                       "--folds_for_test_set=1-1",
#                       "--per_donor"))

# args = parse_args(parser, args =
#                     c("--datasets=Bingren,Greenleaf_brain,Greenleaf_pbmc_bm,Rawlins_fetal_lung,Shendure,Tsankov,Yang_kidney",
#                       "--cancer_types=ColoRect-AdenoCA_cluster_1",
#                       "--cell_number_filter=100",
#                       "--top_features_to_plot_feat_imp=10",
#                       "--ML_model=XGB",
#                       "--annotation=finalized_annotation",
#                       "--seed_range=1-10",
#                       "--feature_importance_method=permutation_importance",
#                       "--folds_for_test_set=1-10",
#                       "--tissues_to_consider=all",
#                       "--robustness_analysis",
#                       "--feat_imp_min_n_robustness=50"))

# args = parse_args(parser, args =
#                     c("--datasets=Bingren,Greenleaf_brain,Greenleaf_pbmc_bm,Rawlins_fetal_lung,Shendure,Tsankov,Yang_kidney",
#                       "--cancer_types=Skin-Melanoma",
#                       "--cell_number_filter=100",
#                       "--top_features_to_plot_feat_imp=10",
#                       "--ML_model=XGB",
#                       "--annotation=finalized_annotation",
#                       "--seed_range=1-10",
#                       "--feature_importance_method=permutation_importance",
#                       "--folds_for_test_set=1-10",
#                       "--tissues_to_consider=all",
#                       "--robustness_analysis",
#                       "--feat_imp_min_n_robustness=50"))

# ,CNS-GBM,ColoRect-AdenoCA,Eso-AdenoCA,Liver-HCC,Lung-AdenoCA,Lung-SCC,Lymph-BNHL,Lymph-CLL,Skin-Melanoma


# args = parse_args(parser, args =
#                     c("--datasets=Bingren,Greenleaf_colon,Greenleaf_pbmc_bm,Tsankov,Shendure",
#                       "--cancer_types=Skin-Melanoma",
#                       "--cell_number_filter=100",
#                       "--ML_model=XGB",
#                       "--annotation=finalized_annotation",
#                       "--robustness_analysis",
#                       "--seed_range=1-10",
#                       "--feature_importance_method=permutation_importance",
#                       "--folds_for_test_set=1",
#                       "--grid_analysis",
#                       "--top_features_to_plot=1",
#                       "--grid_cell_types=skin Melanocyte BR,lung AT2 TS,liver Hepatoblasts SH,stomach Foveolar Cell BR,lung Basal TS,normal_colon Stem GL_Co,bonemarrow B GL_BlBm,mammary_tissue Basal Epithelial (Mammary) BR,cerebrum Astrocytes-Oligodendrocytes SH"))

# args = parse_args(parser, args= c("--cancer_types=Thy-AdenoCA",
#                                   "--datasets=Bingren,Greenleaf_brain,Greenleaf_colon,Greenleaf_pbmc_bm,Rawlins_fetal_lung,Shendure,Tsankov,Yang_kidney",
#                                   "--cell_number_filter=100",
#                                   "--annotation=Greenleaf_colon_remove_cancer_polyp_merge_normal_unaffected",
#                                   "--seed_range=1-10",
#                                   "--top_features_to_plot_feat_imp=10",
#                                   "--feature_importance_method=permutation_importance",
#                                   "--folds_for_test_set=1-10",
#                                   "--tissues_to_consider=all",
#                                   "--robustness_analysis",
#                                   "--feat_imp_min_n_robustness=50"))

# args = parse_args(parser, args= c("--cancer_types=Thy-AdenoCA",
#                                   "--datasets=Bingren,Greenleaf_brain,Greenleaf_colon,Greenleaf_pbmc_bm,Rawlins_fetal_lung,Shendure,Tsankov,Yang_kidney",
#                                   "--cell_number_filter=100",
#                                   "--annotation=Greenleaf_colon_remove_cancer_polyp_merge_normal_unaffected",
#                                   "--seed_range=1-10",
#                                   "--top_features_to_plot_feat_imp=10",
#                                   "--feature_importance_method=permutation_importance",
#                                   "--folds_for_test_set=1-10",
#                                   "--tissues_to_consider=all",
#                                   "--robustness_analysis",
#                                   "--feat_imp_min_n_robustness=50"))

# args = parse_args(parser, args =
#                     c("--cancer_types=Skin-Melanoma,Liver-HCC,ColoRect-AdenoCA,multiple_myeloma,Eso-AdenoCA,CNS-GBM,Lung-AdenoCA,Lung-SCC",
#                       "--ML_model=XGB",
#                       "--annotation=finalized_annotation",
#                       "--robustness_analysis",
#                       "--seed_range=1-10",
#                       "--feature_importance_method=permutation_importance",
#                       "--folds_for_test_set=1-10",
#                       "--grid_analysis",
#                       "--top_features_to_plot=1",
#                       "--grid_cell_types=skin Melanocyte BR,liver Hepatoblasts SH,normal_colon Stem GL_Co,bonemarrow B GL_BlBm,stomach Goblet cells SH,cerebrum Astrocytes-Oligodendrocytes SH,lung AT2 TS,lung Basal TS"))
# args = parse_args(parser, args =
#                     c("--cancer_types=Breast-AdenoCA,Lymph-BNHL,Myeloid-AML,SoftTissue-Leiomyo,Thy-AdenoCA",
#                       "--ML_model=XGB",
#                       "--annotation=finalized_annotation",
#                       "--robustness_analysis",
#                       "--seed_range=1-10",
#                       "--feature_importance_method=permutation_importance",
#                       "--folds_for_test_set=1-10",
#                       "--grid_analysis",
#                       "--top_features_to_plot=1",
#                       "--grid_cell_types=mammary_tissue Basal Epithelial (Mammary) BR,bonemarrow B GL_BlBm,bonemarrow GMP GL_BlBm,stomach Stromal cells SH,thyroid Thyroid Follicular Cell BR"))

# args = parse_args(parser, args= c("--cancer_types=SCLC",
#                                   "--datasets=Bingren,Greenleaf_colon,Greenleaf_pbmc_bm,Shendure,Tsankov,Yang_kidney",
#                                   "--cell_number_filter=100",
#                                   "--annotation=finalized_annotation",
#                                   "--seed_range=1-10",
#                                   "--top_features_to_plot=10,5,2,1",
#                                   "--top_features_to_plot_feat_imp=10,5,2,1",
#                                   "--feature_importance_method=permutation_importance",
#                                   "--folds_for_test_set=1-10",
#                                   "--tissues_to_consider=all",
#                                   "--robustness_analysis",
#                                   "--feat_imp_min_n_robustness=50"))

# args = parse_args(parser, args= c("--cancer_types=msi_high",
#                                   "--datasets=Greenleaf_colon",
#                                   "--cell_number_filter=100",
#                                   "--annotation=Greenleaf_colon_normal_merge_goblet",
#                                   "--seed_range=1-10",
#                                   "--top_features_to_plot_feat_imp=5,2,1",
#                                   "--feature_importance_method=permutation_importance",
#                                   "--folds_for_test_set=1-10",
#                                   "--tissues_to_consider=all",
#                                   "--robustness_analysis",
#                                   "--feat_imp_min_n_robustness=50"))

args = parse_args(parser)

cancer_names = hash("Skin-Melanoma"="Melanoma",
                    "Liver-HCC"="Liver cancer",
                    "ColoRect-AdenoCA"="Colorectal cancer",
                    "Eso-AdenoCA"="Esophageal\ncancer",
                    "CNS-GBM"="Glioblastoma",
                    "Lung-AdenoCA"="Lung\nadenocarcinoma",
                    "Lung-SCC"="Lung squamous\ncell carcinoma",
                    "Breast-AdenoCA"="Breast\nadenocarcinoma",
                    "Lymph-CLL"="Lymphocytic\nleukemia",
                    "Lymph-BNHL"="Non-Hodgkin\nlymphoma",
                    "multiple_myeloma"="Multiple\nMyeloma",
                    "Myeloid-AML"="Acute myeloid\nleukemia",
                    "SoftTissue-Leiomyo"="Leiomyosarcoma",
                    "Thy-AdenoCA"="Thyroid\nadenocarcinoma")

# cell_types = c("lung Basal TS"="Basal, Lung",
#                 "esophagus_mucosa Esophageal Epithelial Cell BR"="Epithelial,\nEsophagus Mucosa", 
#                 "lung Ciliated epithelial cells SH"="Ciliated,\nLung",
#                 "colon_transverse Small Intestinal Enterocyte BR" = "Enterocyte, Colon\nTransverse",
#                 "artery_aorta Smooth Muscle (General) BR" = "Smooth Muscle,\nArtery Aorta",
#                 "lung Mesothelium TS" = "Mesothelium",
#                 "")

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
                                datasets,
                                cell_number_filter,
                                tss_fragment_filter,
                                annotation,
                                ML_model,
                                seed,
                                accumulated_seeds,
                                feature_importance_method, 
                                fold_for_test_set,
                                hundred_kb, 
                                per_donor,
                                n_filter=NULL) {
  dirs = get_relevant_backwards_elim_dirs(cancer_types=cancer_type, 
                                          datasets=datasets,
                                          cell_number_filter=cell_number_filter,
                                          tss_fragment_filter=tss_fragment_filter,
                                          annotation=annotation,
                                          ML_model=ML_model,
                                          hundred_kb=hundred_kb,
                                          fold_for_test_set=fold_for_test_set,
                                          seed=seed,
                                          accumulated_seeds=accumulated_seeds,
                                          per_donor=per_donor)
  
  if (per_donor) {
    counts=data.frame()
  }
  for (dir in dirs) {
    file = paste(dir, "df_for_feature_importance_plots", sep="/")
    if (feature_importance_method != "default_importance") {
      file = paste(file, feature_importance_method, sep="_")
    }
    file = paste(file, "csv", sep=".") 
    file = paste("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/figures", file, sep="/")
    df = as_tibble(read.csv(file))
    title = unlist(strsplit(dir, split ="/"))
    title = title[length(title) - 2]
    df = df[df$num_features %in% top_features_to_plot, ]
    if (per_donor) {
      df["donor"] = unlist(strsplit(title, split = "_"))[2]
      counts = rbind(counts, df)
    } else {
      ggplot_barplot_helper(df, title, savepath=paste("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/figures", dir, sep="/"), 
                          ylab=gsub("_", " ", feature_importance_method), 
                          y=feature_importance_method)
    }
  }
  
  if (per_donor) {
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
    fp = paste(dir, "per_donor_predictions.csv", sep="/")
    write.csv(counts, fp)
    bin_width = 0.05
    break_points <- seq(min(counts$score), max(counts$score) + bin_width, bin_width)
    labels <- sprintf("%.2f-%.2f", head(break_points, -1), tail(break_points, -1))
    if (!is.null(n_filter)) {
      dist = table(counts[["features"]])
      keep = names(dist)[dist >= n_filter]
      counts = counts %>% filter(features %in% keep)
    }
    df <- counts %>%
            mutate(bin = cut(score, breaks = break_points,labels=labels,
                             right = FALSE)) %>%
            group_by(bin, features) %>%
            summarise(count = n()) %>%
            arrange(desc(count))
    
    total_counts <- df %>%
                      group_by(features) %>%
                      summarise(total = sum(count)) %>%
                      arrange(desc(total))
    df = df %>% 
          left_join(total_counts) %>%
          mutate(features = paste0(features, " (n=", total, ")")) 
    ggplot(df, aes(x = bin, y = count, fill = features)) +
      geom_bar(stat = "identity", position = "stack") +
      geom_text(aes(label = count), position = position_stack(vjust = 0.5),
                size = 3) +
      labs(title = paste(cancer_type, "per donor"), x = "Validation R^2", 
           y = "Count") +
      scale_x_discrete(limits = levels(df$bin)) +
      theme(axis.text.x = element_text(size=6))
    
    # ggplot(counts, aes(x = score, fill = features)) +
    #   geom_histogram(position = "stack", bins = 10) +
    #   labs(title = paste(cancer_type, "per donor"), x = "Validation R^2", y = "Count")
  }
}

construct_robustness_boxplots <- function(df, x, y, title, savepath, savefile,
                                          facet_var, xlabel="", plot_fold=F, 
                                          n_name=NULL, width=12, height=8) {
  if (plot_fold) {
    outlier_shape = NA
  } else {
    outlier_shape = 19
  }
  
  plots <- list()
  # df = df %>% 
  #       ungroup() %>%
  #        mutate("{y}" := unname(cell_types[df %>% pull(!!sym(y))]))
  
  for (level in unique(df[[facet_var]])) {
    df_filtered <- df %>% 
                     filter(!!sym(facet_var) == level)
    
    if (x == "permutation_importance") {
      unique_combos = unique(df_filtered[, c("features", "n_feature", "med_imp")])
      sorted_features = unique_combos %>% 
        arrange(desc(n_feature), desc(med_imp)) %>%
        pull(features)
      df_filtered = df_filtered %>% 
          filter(features %in% unique(sorted_features)[1:5])
    }
    
    df_compressed = df_filtered %>%
                      group_by(!!sym(y)) %>%
                      summarise(med_x = median(!!sym(x)))
    df_filtered = df_filtered %>% 
                    left_join(df_compressed) %>%
                    arrange(!!sym(n_name), med_x)
    top = unique(df_filtered %>% pull(!!sym(y)))
    top = top[length(top)]
    
    df_filtered = df_filtered %>% 
                    ungroup() %>%
                    mutate(color=ifelse(!!sym(y)==top, "highlight", "other"))
    
    # ifelse(!!sym(n_name)==max(!!sym(n_name)) &
    #          med_x==max(med_x), "highlight", 
    #        "other")
    
    ordered_levels <- df_filtered %>% pull(!!sym(y)) %>% unique()
    df_filtered <- df_filtered %>% 
                      mutate(y_reordered = factor(!!sym(y), 
                                                  levels = ordered_levels))
    xlim_lower = min(df[[x]])
    xlim_upper = max(df[[x]])
    color = rep("#000000", length(unique(df_filtered[[y]])) - 1)
    color = append(color, "#EE4B2B")
    
    p <- ggplot(df_filtered) +
            geom_boxplot(aes(x = !!sym(x), y_reordered, fill=color), lwd = 1.0, 
                             outlier.shape = outlier_shape) +
            geom_text(aes(x = x_position + xlim_upper / 10,
                          y = y_reordered),
                          label = paste0("n=", df_filtered[[n_name]])) +
            ggtitle(paste(level, "features")) +
            scale_fill_manual(values = c("highlight" = "#EE4B2B",
                                         "other" = "#A9A9A9")) +
            theme_classic() +
            theme(
                legend.position="none",
                strip.background = element_blank(),
                strip.text.x = element_blank(),
                plot.title = element_text(hjust = 0.5),
                axis.text.y = element_text(size = 15, colour = color),
                axis.text.x = element_text(size = 14),
                axis.title.x=element_blank(),
                axis.title.y=element_blank()
              ) +
            xlim(xlim_lower - xlim_lower / 10, xlim_upper + xlim_upper / 10)
    plots[[as.character(level)]] <- p
  }
  
  if (length(plots) > 1) {
    plot <- wrap_plots(plots, ncol=length(unique(df[[facet_var]])), 
                                nrow=1)
    plot <- plot + 
              plot_annotation(title = title) +
              plot_layout(guides = 'collect') +
              plot_annotation(caption = xlabel)
  } else {
    plot <- plots[[1]] +
              xlab(xlabel) +
              theme(axis.title.x=element_text(size=15),
                    axis.text.x = element_text(size=12))
  }

  if (plot_fold) {
    plot = plot + 
      geom_point(aes(x = reorder_within(x=!!sym(x), 
                                      by=!!sym(y), 
                                      within=!!sym(facet_var), 
                                      median), 
                   y = !!sym(y), 
                   fill=fold,
                   alpha=scale(scatac_counts),
                   size=scale(mut_counts)), pch=21, stroke=NA) +
      scale_fill_manual(values = c("#800000", "#808000",
                                   "#469990", "#000000",
                                   "#e6194B", "#ffe119",
                                   "#3cb44b", "#4363d8",
                                   "#f032e6", "#a9a9a9",
                                   "#ffd8b1", "#aaffc3"))
  }
  print(plot)
  ggsave(paste(savepath, savefile, sep="/"), 
         width = width * length(unique(df[[facet_var]])), 
         height = height, plot)
}


construct_robustness_barplots <- function(df, x, y, title, add_to_pos) {
  df[x] = 100 * df[x]
  df = df %>% 
        group_by(!!sym(y)) %>%
        summarise(mean_se(test_set_perf)) %>%
        mutate(color=ifelse(y == max(y), "highlight","other"))
  
  xlim_lower = min(df[["ymin"]])
  xlim_upper = max(df[["ymax"]])
  x_breaks <- pretty(df$y, n = 3)
  plot <- ggplot(df) +
    geom_col(aes(x = y, y = factor(top_feature,
                                   levels=rev(c("skin-Melanocyte-BR",
                                               "liver-Hepatoblasts-SH",
                                               "normal_colon-Stem-GL_Co",
                                               "bonemarrow-B-GL_BlBm",
                                               "stomach-Goblet-cells-SH",
                                               "cerebrum-Astrocytes-Oligodendrocytes-SH",
                                               "lung-AT2-TS",
                                               "lung-Basal-TS"
                                               ))),
                 fill=color), lwd=1.2) +
    geom_errorbarh(aes(y = top_feature, xmin = ymin, xmax = ymax), linewidth=2) +  
    coord_cartesian(xlim = c(xlim_lower, xlim_upper)) +
    scale_fill_manual(values = c("highlight" = "#EE4B2B", "other" = "#A9A9A9")) +
    ylab("") +
    ggtitle(title) +
    scale_x_continuous(breaks = x_breaks) +   
    theme_bw() +
    theme(
      legend.position="none",
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 40),
      axis.text.x = element_text(vjust = 0.5, hjust=1, size=40),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x=element_blank(),
      panel.grid.major.y=element_blank(),
      panel.grid.minor.x=element_blank(),
      panel.grid.major.x = element_line(size = 1),
      plot.margin = margin(t=0,r=1,b=0,l=0, unit = "cm")
    ) 
  return(plot)
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

# plot_accumulated_feature_importance <- function(df_feature_importances_all_seeds,
#                                                 savepath,
#                                                 cancer_type) {
#   df_accumulated_imp = df_feature_importances_all_seeds %>% 
#     group_by(features, num_features) %>%
#     summarise(sum(importance))
#   colnames(df_accumulated_imp)[3] = "importance"
#   dir.create(savepath, recursive = T)
#   
#   ggplot_barplot_helper(df=df_accumulated_imp, 
#                         title=cancer_type, 
#                         savepath=savepath,
#                         accumulated_imp=T)
# }

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
                                        feature_importance_method,
                                        hundred_kb, 
                                        grid_analysis=F,
                                        grid_cell_type=NULL) {
  df = tibble(top_feature = character(0),
              top_n = integer(0),
              test_set_perf = double(0),
              seed = integer(0),
              fold = integer(0))
  
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
                                                hundred_kb = hundred_kb,
                                                grid_analysis=grid_analysis,
                                                grid_cell_type=grid_cell_type)
        
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

        # #### TEMP ####
        # total_num_features=1
        # ###############
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
          if (grid_analysis) {
            top_feature = grid_cell_type
          }
          df = df %>% add_row(top_feature = top_feature,
                              top_n = top_features_to_plot[idx],
                              test_set_perf = perf,
                              seed = seed, 
                              fold = fold)
          idx = idx + 1
        }
      }
    }
  }
  df["fold"] = as.factor(df[["fold"]])
  return(df)
}

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


get_and_plot_scatac_and_mutation_counts_per_fold <- function(cancer_type,
                                                    folds_for_test_set,
                                                    datasets,
                                                    cell_number_filter,
                                                    tss_fragment_filter,
                                                    annotation, 
                                                    ML_model,
                                                    hundred_kb,
                                                    tissues_to_consider) {
  scatac_counts_plots <- list()
  scatac_counts <- c()
  mut_counts_plots <- list()
  mut_counts <- c()
  for (fold in folds_for_test_set) {
    # for (tss_filter in tss_fragment_filter) {
    scATAC_sources = construct_sources_string(datasets)
    dir = construct_dir(scATAC_sources,
                        cell_number_filter,
                        tss_fragment_filter,
                        annotation,
                        1,
                        fold,
                        ML_model,
                        cancer_type,
                        hundred_kb,
                        tissues_to_consider)
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
                                                      hundred_kb = hundred_kb)
    if (feature_importance_method != "default_importance") {
      model_pattern = paste("^model_iteration_[0-9]+_feature_importance",
                            feature_importance_method, sep="_")
      model_pattern = paste(model_pattern, "pkl", sep=".")
    } else {
      model_pattern = "^model_iteration_[0-9]+\\.pkl"
    }
    total_num_features = list.files(backwards_elim_dir,
                                    pattern=model_pattern)
    total_num_features = length(total_num_features)
    fp = paste(backwards_elim_dir, 
               paste("top_features_iteration", total_num_features - 1,
               "by_permutation_importance.txt", sep="_"), sep="/")
    
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
    
    scatac_counts = append(scatac_counts, sum(X_test_top_feature))
    mut_counts = append(mut_counts, sum(y_test))
  }
  
  pdf("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/figures/per_fold_scatac_counts.pdf", width = 50, height = 20)
  do.call(grid.arrange, c(scatac_counts_plots, nrow=2)) 
  dev.off()
  pdf("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/figures/per_fold_mut_counts.pdf", width = 50, height = 20)
  do.call(grid.arrange, c(mut_counts_plots, nrow=2))
  dev.off()
  return(list(scatac_counts, mut_counts))
}


cancer_types = args$cancer_types
cancer_types = unlist(strsplit(cancer_types, split = ","))
top_features_to_plot = args$top_features_to_plot
top_features_to_plot = as.integer(unlist(strsplit(top_features_to_plot, 
                                                  split = ",")))
cell_number_filter = args$cell_number_filter
tss_fragment_filter = unlist(strsplit(args$tss_fragment_filter, split = ","))
annotation = args$annotation
tissues_to_consider = strsplit(args$tissues_to_consider,  split=",")
tissues_to_consider = paste(unlist(tissues_to_consider), collapse="_")
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
grid_analysis = args$grid_analysis
if (grid_analysis) {
  grid_cell_types = unlist(strsplit(args$grid_cell_types, split=","))
  grid_cell_types = gsub(" ", "-", grid_cell_types)
}
if (!grid_analysis) {
  datasets = unlist(strsplit(args$datasets, split = ","))
  datasets = sort(datasets)
} else {
  datasets=""
}

folds_for_test_set = args$folds_for_test_set
folds_for_test_set = unlist(strsplit(args$folds_for_test_set, split = "-"))
folds_for_test_set = seq(folds_for_test_set[1], folds_for_test_set[2])
feat_imp_min_n_robustness = args$feat_imp_min_n_robustness
plot_fold_on_test_set_plot = args$plot_fold_on_test_set_plot
hundred_kb = args$hundred_kb
per_donor = args$per_donor

if (!robustness_analysis) {
  for (seed in seed_range) {
    for (fold in folds_for_test_set) {
      for (cancer_type in cancer_types) {
        construct_bar_plots(cancer_type, 
                            datasets,
                            cell_number_filter,
                            tss_fragment_filter,
                            annotation,
                            ML_model,
                            seed,
                            accumulated_seeds=F,
                            feature_importance_method=feature_importance_method,
                            fold_for_test_set=fold,
                            hundred_kb=hundred_kb,
                            per_donor=per_donor)
      }
    }
  }
} else {
  grid_plots = list()
  i = 1
  for (cancer_type in cancer_types) {
    # for (tss_filter in tss_fragment_filter) {
    if (!grid_analysis) {
      scATAC_source = paste("cell_number_filter", cell_number_filter, sep="_")
      scATAC_sources = construct_sources_string(datasets)
      scATAC_source = paste("scATAC_source", scATAC_sources, 
                            scATAC_source, sep="_")
      
    } else {
      scATAC_source = ""
    }
    
    if (tss_fragment_filter != "-1") {
      scATAC_source = paste(scATAC_source, "tss_fragment_filter",
                            tss_fragment_filter, sep="_")
    }
    
    # if (tss_filter != "-1") {
    #   scATAC_source = paste(scATAC_source, "tss_fragment_filter",
    #                         tss_filter, sep="_")
    # }
    if (tissues_to_consider != "all") {
      scATAC_source = paste(scATAC_source, "tissues_to_consider", tissues_to_consider, 
                            sep="_")
    }
    if (scATAC_source == "") {
      scATAC_source = paste("annotation", annotation, sep="_")
    } else {
      scATAC_source = paste(scATAC_source, "annotation", annotation, sep="_")
    }
    if (!grid_analysis) {
      savepath = get_relevant_backwards_elim_dirs(cancer_types=cancer_type, 
                                                # combined_datasets=combined_datasets,
                                                tissues_to_consider=tissues_to_consider,
                                                datasets=datasets,
                                                cell_number_filter=cell_number_filter,
                                                tss_fragment_filter=tss_fragment_filter,
                                                annotation=annotation,
                                                ML_model=ML_model,
                                                hundred_kb=hundred_kb,
                                                accumulated_seeds=T)
      savepath = paste("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/figures", savepath, sep="/")
      # savepath = paste("/home/mdanb/research/mount_sinai/epigenetic_effects_on_mutational_landscape/figures", savepath, sep="/")
      
      dir.create(path=savepath, recursive = T)
      dirs = list.dirs(paste("../../figures", "models", ML_model, cancer_type,
                             sep="/"), recursive = F)
    } else {
      dirs = c()
      for (cell_type in grid_cell_types) {
        dirs = append(dirs, list.dirs(paste("../../figures", "models", 
                      ML_model, paste(cancer_type, cell_type, sep="_"), sep="/"),
                      recursive = F))
      }
    } 
    
    # dirs = list.dirs(paste("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/figures", "models", ML_model, cancer_type, 
    #                        sep="/"), recursive = F)
    combos = expand.grid(seed = seed_range, fold = folds_for_test_set)
    seed_fold_for_test_combinations = apply(combos, 1, function(row) {
      paste(scATAC_source, "seed", row["seed"], "fold_for_test_set", row["fold"], 
            sep="_")
    })
    
    if (hundred_kb) {
      seed_fold_for_test_combinations = paste("interval_ranges_100kb", 
                                              seed_fold_for_test_combinations,
                                              sep = "_")
    }
    all_seeds_dirs = dirs[basename(dirs) %in% seed_fold_for_test_combinations]

    if (grid_analysis) {
      df = data.frame()
      for (grid_cell_type in grid_cell_types) {
        df = rbind(df, construct_all_seeds_test_df(top_features_to_plot=top_features_to_plot,
                                                   seed_range=seed_range,
                                                   skip_seeds_robustness=skip_seeds_robustness,
                                                   cancer_type=paste(cancer_type, grid_cell_type, sep="_"),
                                                   datasets=datasets,
                                                   cell_number_filter=cell_number_filter,
                                                   tss_fragment_filter=tss_fragment_filter, 
                                                   annotation=annotation,
                                                   tissues_to_consider=tissues_to_consider,
                                                   ML_model=ML_model,
                                                   folds_for_test_set=folds_for_test_set,
                                                   feature_importance_method=feature_importance_method,
                                                   hundred_kb=hundred_kb,
                                                   grid_analysis=grid_analysis,
                                                   grid_cell_type=grid_cell_type))
      }
      plot = construct_robustness_barplots(df, x="test_set_perf", 
                                           y="top_feature",
                                           title=cancer_names[[cancer_type]])
      grid_plots[[i]] = plot
      i = i + 1
    } else {
      df_feature_importances_all_seeds = 
        construct_df_feature_importances_all_seeds(all_seeds_dirs,
                                                   feature_importance_method,
                                                   top_features_to_plot)
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
                                       feature_importance_method=feature_importance_method,
                                       hundred_kb=hundred_kb)
      # y_position is for plotting number of times feature appears at the top of
      # the boxplot. 
      df_feat_imp = df_feature_importances_all_seeds %>% 
        # group_by(num_features, seed, fold_for_test_set) %>%
        group_by(num_features, features) %>%
        mutate(n_feature = n(), 
               med_imp = median(permutation_importance), 
               x_position = max(permutation_importance)) %>%
        filter(num_features %in% top_features_to_plot_feat_imp) 
      # unique_combos = unique(df_feat_imp[, c("features", "n_feature", "med_imp")])
      # sorted_features = unique_combos %>% 
      #                         arrange(desc(n_feature), desc(med_imp)) %>%
      #                         pull(features)
      # print(sorted_features)
      # df_feat_imp_top_5 = df_feat_imp %>% 
      #   # filter(n_feature >= feat_imp_min_n_robustness)%>%
      #   filter(features %in% unique(sorted_features)[1:5])
      savefile = paste0(cancer_type, "_feature_importance_with_",
                        paste(top_features_to_plot_feat_imp, collapse="_"),
                        "_features_", "top_5_features.pdf")
      
      construct_robustness_boxplots(df=df_feat_imp, 
                                    x="permutation_importance", 
                                    y="features", 
                                    title=cancer_names[[cancer_type]], 
                                    savepath=savepath,
                                    savefile=savefile, 
                                    n_name="n_feature", 
                                    facet_var="num_features",
                                    xlabel="Feature Importance",
                                    width=11,
                                    height=7)
      df_test = df %>% 
        group_by(top_n, top_feature) %>%
        mutate(n_top_feature = n(), x_position = max(test_set_perf))
      
      if (plot_fold_on_test_set_plot) {
        l = get_and_plot_scatac_and_mutation_counts_per_fold(cancer_type,
                                                             folds_for_test_set,
                                                             datasets,
                                                             cell_number_filter,
                                                             tss_fragment_filter,
                                                             annotation, 
                                                             ML_model,
                                                             hundred_kb,
                                                             tissues_to_consider)
        scatac_counts = l[[1]]
        mut_counts = l[[2]]
        df_test["scatac_counts"] = scatac_counts[df[["fold"]]]
        df_test["mut_counts"] = mut_counts[df[["fold"]]]
      }
      
      savefile = paste0("test_set_boxplots_with_",
                        paste(top_features_to_plot_feat_imp, collapse="_"),
                        "_features.png")
      
      construct_robustness_boxplots(df=df_test, 
                                    x="test_set_perf", 
                                    y="top_feature", 
                                    title=cancer_type, 
                                    savepath=savepath,
                                    savefile=savefile, 
                                    n_name="n_top_feature",
                                    facet_var="top_n",
                                    xlabel="Variance Explained, Test Set",
                                    width=8, 
                                    height=5)
      
      df_val = df_feature_importances_all_seeds %>% 
        group_by(num_features, seed, fold_for_test_set) %>%
        mutate(max_feature_importance=max(permutation_importance)) %>%
        filter(permutation_importance == max_feature_importance) %>%
        ungroup() %>%
        group_by(num_features, features) %>%
        mutate(n_feature = n(), x_position = max(score))
      
      savefile = paste0("validation_boxplots_with_",
                        paste(top_features_to_plot_feat_imp, collapse="_"),
                        "_features.png")
      
      construct_robustness_boxplots(df=df_val, 
                                    x="score", 
                                    y="features", 
                                    title=cancer_type, 
                                    savepath=savepath,
                                    savefile=savefile, 
                                    n_name="n_feature",
                                    facet_var="num_features",
                                    xlabel="Variance Explained, Validation Set",
                                    width=12, 
                                    height=5)
      
    }
  }
}

if (grid_analysis) {
  plot <- wrap_plots(grid_plots, ncol=8, nrow=1)
  plot <- plot + 
    # plot_layout(guides = 'collect') +
    plot_annotation(caption = "Test Set Variance Explained (%)",
                    theme = theme(plot.caption = element_text(size = 40, hjust=0.5)))
  ggsave("../../figures/grid_analysis.pdf", 
         width = 6 * 8, height = 10, limitsize = FALSE)
}



# plot_list <- lapply(grid_plots, function(p) {
#   p + plot_spacer(width=0.2) 
# })


# construct_boxplots(df=df_test, 
#                    x="test_set_perf", 
#                    y="top_feature", 
#                    title=cancer_type,
#                    savepath=savepath,
#                    savefile=savefile,
#                    facet_var="top_n",
#                    xlabel="Variance Explained, Test Set")

# df_val = df_feature_importances_all_seeds %>% 
#           # group_by(num_features, seed, fold_for_test_set) %>%
#           group_by(num_features, features) %>%
#           mutate(n_feature = n(), 
#                  x_position = max(score)) %>%
#           filter(num_features %in% top_features_to_plot_feat_imp) 

# df_feat_imp = df_feature_importances_all_seeds %>% 
#   # group_by(num_features, seed, fold_for_test_set) %>%
#   group_by(num_features, features) %>%
#   mutate(n_feature = n(), 
#          x_position = max(permutation_importance)) %>%
#   filter(num_features %in% top_features_to_plot_feat_imp) 

# construct_robustness_boxplots(df=df_feat_imp_at_least_n, 
#                               x="permutation_importance", 
#                               y="features", 
#                               title=cancer_type, 
#                               savepath=savepath,
#                               savefile=savefile, 
#                               n_name="n_feature", 
#                               facet_var="num_features",
#                               xlabel="Permutation Importance")

# from = as.character(unique(df[[facet_var]]))
# to = paste("top", from, "features")
# 
# names(to) <- from
# df[[facet_var]] <- factor(df[[facet_var]], levels = unique(df[[facet_var]]))


# custom_label_function <- function(labels) {
#   cleaned_labels <- str_replace(labels, "___10", "")
#   cleaned_labels <- str_replace(cleaned_labels, "___1", "")
#   wrapped_labels <- label_wrap_gen(width = 1)(cleaned_labels)
#   return(wrapped_labels)
# }

# grayscale_palette <- gray(seq(0, 0.9, length.out = 5))


# } else if (test_set_boxplots) {
# #   reorder_within <- function(x, by, within, fun = median, sep = "___") {
# #     new_x <- paste(x, within, sep = sep)
# #     ordered_factor <- stats::reorder(new_x, by, FUN = fun)
# #     levels_reversed <- rev(levels(ordered_factor))
# #     return(factor(ordered_factor, levels = levels_reversed))
# # }
#   
#   # plot +
#   #   geom_boxplot(aes(x = reorder_within(x=!!sym(x), by=!!sym(y),
#   #                                       within=!!sym(facet_var), 
#   #                                       fun=median), y = !!sym(y)), 
#   #                lwd=1.2, outlier.shape = outlier_shape) 
#   # 
#   # df <- df %>%
#   #   group_by(num_features) %>%
#   #   mutate(reordered_features = reorder(features, permutation_importance, 
#   #                                       FUN = median, 
#   #                                       decreasing = TRUE)) %>%
#   #   ungroup()
#   # 
#   # ggplot(df) +
#   #   facet_wrap(~num_features, scales = "free_x") +
#   #   geom_boxplot(aes(x = reordered_features, y = permutation_importance), lwd = 1.2, outlier.shape = outlier_shape) +
#   #   theme(
#   #     legend.position = "none",
#   #     strip.background = element_blank(),
#   #     strip.text.x = element_blank(),
#   #     plot.title = element_text(hjust = 0.5),
#   #     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
#   #     axis.text.y = element_text(size = 45),
#   #     axis.title.x = element_text(size = 20)
#   #   )
#   # 
#   # 
#   # ggplot(df) +
#   #   facet_wrap(~num_features, scales = "free_x") +
#   #   geom_boxplot(aes(x = reorder(features, permutation_importance, 
#   #                                FUN = median, decreasing=T), y = permutation_importance), 
#   #                lwd=1.2, outlier.shape = outlier_shape)+
#   #   theme(
#   #     legend.position="none",
#   #     strip.background = element_blank(),
#   #     strip.text.x = element_blank(),
#   #     plot.title = element_text(hjust = 0.5),
#   #     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#   #     axis.text.y = element_text(size = 45),
#   #     axis.title.x=element_text(size=20),
#   #     # axis.text.y=element_blank()
#   #   ) 
#   
#   
#   
#     # geom_text(aes(x = reorder_within(x=!!sym(x), 
#     #                                  by=!!sym(y),
#     #                                  within=!!sym(facet_var), 
#     #                                  median), y = y_position, 
#     #               label = paste0("n=", !!sym(n_name))), 
#     #           vjust = -0.5)
#   
#   # reorder_within <- function(x, by1, by2, within, FUN, sep = "___") {
#   #   new_x <- paste(x, within, sep = sep)
#   #   ordered_factor <- factor(new_x, levels = unique(new_x)[order(by1, decreasing = TRUE)])
#   #   for 
#   #   ties <- duplicated(by1) | duplicated(by1, fromLast = TRUE)
#   #   if (any(ties)) {
#   #     for (tie_group in unique(ordered_factor[ties])) {
#   #       indices <- which(ordered_factor == tie_group)
#   #       reordered_within_group <- order(by2[indices], decreasing = TRUE)
#   #       levels(ordered_factor)[indices] <- levels(ordered_factor)[indices[reordered_within_group]]
#   #     }
#   #   }
#   #   return(factor(ordered_factor, levels = levels(ordered_factor)))
#   # }
#   # 
#   # plot <- ggplot(df) +
#   #           geom_boxplot(aes(x = reorder_within(x=!!sym(x), by1=n_name, 
#   #                                               by2=!!sym(y),
#   #                                               within=!!sym(facet_var), 
#   #                                               FUN=median), y = !!sym(y)), 
#   #                        lwd=1.2, outlier.shape = outlier_shape) #+
#     # geom_text(aes(x = reorder_within(x=!!sym(x), by1=n_name, 
#     #                                  by2=!!sym(y),
#     #                                  within=!!sym(facet_var), 
#     #                                  median), y = y_position, 
#     #               label = paste0("n=", !!sym(n_name))), 
#     #           vjust = -0.5)
#   
# } else {
# reorder_within_option_3 <- function(y) {
#   # levels_reversed <- rev(levels(ordered_factor))
#   return(factor(y, levels = rev(c("skin-Melanocyte-BR",
#                                        "liver-Hepatoblasts-SH",
#                                        "normal_colon-Stem-GL_Co",
#                                        "stomach-Foveolar-Cell-BR",
#                                        "cerebrum-Astrocytes-Oligodendrocytes-SH",
#                                        "lung-AT2-TS",
#                                        "lung-Basal-TS",
#                                        "mammary_tissue-Basal-Epithelial-(Mammary)-BR",
#                                        "bonemarrow-B-GL_BlBm"))))
# }
#




# plot <- ggplot(df) +
#           geom_boxplot(aes(x = !!sym(x), y = !!sym(y),),
#                                                           
#            lwd=1.2, outlier.shape = outlier_shape) +#++
# facet_wrap(as.formula(paste0("~", facet_var)), nrow=1, 
#            scales = "free_x",
#            labeller = as_labeller(to)) +
#   # xlab(xlabel) +
#   ylab("") +
#   xlab(xlabel) +
#   ggtitle(title) +
#   scale_x_continuous(breaks = seq(round(min(df[[x]]), 2),
#                                   round(max(df[[x]]), 2), by = 0.05)) +
#   # scale_y_discrete(labels = custom_label_function) +
#   
#   # scale_color_manual(values = grayscale_palette) +
#   # scale_y_discrete(labels = custom_label_function) + 
#   scale_y_discrete(labels = rev(c("Foveolar",
#                                   "Melanocyte",
#                                   "Colon, Stem",
#                                   "Mammary Tissue, Basal",
#                                   "Lung, Basal",
#                                   "Lung, AT2",
#                                   "Hepatoblast",
#                                   "Astrocyte/Oligodendrocyte",
#                                   "B"))) +
#   # theme_classic() +
#   theme(
#     strip.background = element_blank(),
#     strip.text.x = element_blank(),
#     plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
#     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=25),
#     # axis.text.y = element_text(size = 25),
#     axis.title.x=element_text(size=20),
#     axis.text.y=element_blank()
#   )
# plot


# plot <- plot +
#   facet_wrap(as.formula(paste0("~", facet_var)), nrow=1,
#              scales = "free_x",
#              labeller = as_labeller(to)) +
# xlab(xlabel) +
# ylab(ylabel) +
#   xlab(xlabel) +
#   ggtitle(title) +
#   # scale_x_continuous(breaks = seq(round(min(df[[y]]), 2),
#   #                                 round(max(df[[y]]), 2), by = 0.05)) +
#   # scale_y_discrete(labels = custom_label_function) +
#   
#   # scale_color_manual(values = grayscale_palette) +
#   # scale_y_discrete(labels = custom_label_function) +
#   # theme_classic() +
#   scale_x_discrete(labels = custom_label_function) + # Adjust 'width' as needed
#   
#   theme(
#     strip.background = element_blank(),
#     strip.text.x = element_blank(),
#     plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
#     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
#     # axis.text.y = element_text(size = 20),
#     # axis.text.y=element_blank()
#   )


# reorder_within_option_3 <- function(y, by, within, sep = "___") {
#   new_y <- paste(y, within, sep = sep)
#   return(factor(new_y, levels = c("Skin, Melanocyte", 
#                                   "Liver HCC",
#                                   "Colon, Stem",
#                                   "Stomach, Foveolar",
#                                   "Frontal Cortex, Oligodendrocyte Precur")))
# }

# plot <- ggplot(df) +
#         geom_boxplot(aes(x = reorder_within(x=!!sym(x), 
#                                             by=!!sym(y), 
#                                             within=!!sym(facet_var), 
#                                             median), y = !!sym(y)), 
#                      lwd=1.2, outlier.shape = outlier_shape) +
#         geom_text(aes(x = reorder_within(x=!!sym(x), 
#                                          by=!!sym(y), 
#                                          within=!!sym(facet_var), 
#                                          median), y = y_position, label = paste0("n=", !!sym(n_name))), vjust = -0.5) +
#         facet_wrap(as.formula(paste0("~", facet_var)), nrow=1, 
#                    scales = "free_x",
#                    labeller = as_labeller(to)) +
#         ylab(ylabel) +
#         xlab("") +
#         ggtitle(title) +
#         scale_y_continuous(breaks = seq(round(min(df[[y]]), 2),
#                                         round(max(df[[y]]), 2), by = 0.05)) +
#         scale_color_gradient(low = "black", high = "white") +
#         scale_x_discrete(labels = label_wrap_gen(width = 10)) + # Adjust 'width' as needed
#         theme_classic() +
#         theme(
#             strip.background = element_blank(),
#             strip.text.x = element_blank(),
#             plot.title = element_text(hjust = 0.5),
#             axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
#         )

# plot <- ggplot(df) +
#           geom_boxplot(aes(x = reorder_within_option_2(x=!!sym(x), 
#                                               by1=n_feature,
#                                               by2=!!sym(y),
#                                               # by=!!sym(y), 
#                                               within=!!sym(facet_var), 
#                                               median), y = !!sym(y)), 
#                        lwd=1.2, outlier.shape = outlier_shape) +
#           geom_text(aes(x = reorder_within_option_2(x=!!sym(x), 
#                                                    by1=n_feature,
#                                                    by2=!!sym(y), 
#                                                    within=!!sym(facet_var), 
#                                                    median), y = y_position, label = 
#                                                    paste0("n=", !!sym(n_name))), 
#                                                    vjust = -0.5) +
#           facet_wrap(as.formula(paste0("~", facet_var)), nrow=1, 
#                      scales = "free_x",
#                      labeller = as_labeller(to)) +
#           ylab(ylabel) +
#           xlab("") +
#           ggtitle(title) +
#           scale_y_continuous(breaks = seq(round(min(df[[y]]), 2),
#                                           round(max(df[[y]]), 2), by = 0.05)) +
#           # scale_color_manual(values = grayscale_palette) + # Apply grayscale palette
#           scale_x_discrete(labels = custom_label_function) + 
#           theme_classic() +
#           theme(
#             strip.background = element_blank(),
#             strip.text.x = element_blank(),
#             plot.title = element_text(hjust = 0.5),
#             axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
#           ) +
#           coord_flip()
# 
#           # theme(plot.title = element_text(hjust = 0.5),
#           #       axis.text.x = element_blank())
#   print(plot)


# if (robustness_accumulated_feature_importance_barplot) {
#   plot_accumulated_feature_importance(df_feature_importances_all_seeds,
#                                       savepath,
#                                       cancer_type)
# }

# df_feat_imp = df_feature_importances_all_seeds %>% 
#   group_by(num_features, seed, fold_for_test_set) %>%
#   mutate(max_feature_importance=max(permutation_importance)) %>%
#   filter(permutation_importance == max_feature_importance) %>%
#   ungroup() %>%
#   group_by(num_features, features) %>%
#   mutate(n_feature = n(), y_position = max(permutation_importance)) %>%
#   filter(num_features %in% top_features_to_plot_feat_imp)
# 
# construct_boxplots(df=df_feat_imp, x="features", 
#                    y="permutation_importance", 
#                    color="features", title=cancer_type, savepath=savepath,
#                    savefile="feature_importance_top_feat_only.png", 
#                    n_name="n_feature", facet_var="num_features",
#                    ylabel="Permutation Importance")
# if (grid_analysis) {
#   savepath = c()
#   for (grid_cell_type in grid_cell_types) {
#     savepath = append(savepath, get_relevant_backwards_elim_dirs(cancer_types=paste(cancer_type, 
#                                                                    grid_cell_type, 
#                                                                    sep="_"), 
#                                                 # combined_datasets=combined_datasets,
#                                                 tissues_to_consider=tissues_to_consider,
#                                                 datasets=datasets,
#                                                 cell_number_filter=cell_number_filter,
#                                                 tss_fragment_filter=tss_fragment_filter,
#                                                 annotation=annotation,
#                                                 ML_model=ML_model,
#                                                 hundred_kb=hundred_kb,
#                                                 seed=
#                                                 accumulated_seeds=F))
#   }
#}


