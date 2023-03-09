library(ggplot2)
library(tidytext)
library(magrittr)
library(RColorBrewer)
library(dplyr)
library(optparse)
library(stringr)
source("/ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/utils.R")

parser <- OptionParser()

# parser <- add_option(parser, c("--bing_ren"), action="store_true",
#                                default=F)
# parser <- add_option(parser, c("--tsankov"), action="store_true",
#                      default=F)
# parser <- add_option(parser, c("--shendure"), action="store_true",
#                      default=F)
parser <- add_option(parser, c("--datasets"), type="character")
parser <- add_option(parser, c("--cancer_types"), type="character")
parser <- add_option(parser, c("--cell_number_filter"), type="integer")
parser <- add_option(parser, c("--pie_chart"), action="store_true", default=F)
# parser <- add_option(parser, c("--bar_plot_num_features"), type="integer")
parser <- add_option(parser, c("--tss_fragment_filter"), 
                     type="character", default="-1")
parser <- add_option(parser, c("--ML_model"), type="character")
parser <- add_option(parser, c("--num_iter_skips"), 
                     type="integer", default=5)
parser <- add_option(parser, c("--tissues_to_consider"), 
                     type="character", default="all")
parser <- add_option(parser, c("--annotation"), 
                     type="character", default="default_annotation")
parser <- add_option(parser, c("--waddell_sarc_biph"), action="store_true",
                     default=F)
parser <- add_option(parser, c("--waddell_sarc"), action="store_true",
                     default=F)
parser <- add_option(parser, c("--waddell_sarc_tsankov_sarc"), 
                     action="store_true",
                     default=F)
parser <- add_option(parser, c("--waddell_sarc_biph_tsankov_sarc_biph"), 
                     action="store_true",
                     default=F)
parser <- add_option(parser, c("--iters_dont_skip"), default="18")

args = parse_args(parser)
# args = parse_args(parser, args = c("--cancer_types=Biliary-AdenoCA",
#                                    "--all_cells", "--cell_number_filter=1",
#                                    "--bing_ren", "--tsankov", "--shendure"))

construct_backwards_elim_dir <- function(cancer_type, scATAC_source, 
                                         cell_number_filter,
                                         tss_fragment_filter, 
                                         waddell_sarc_biph,
                                         waddell_sarc,
                                         waddell_sarc_tsankov_sarc,
                                         waddell_sarc_biph_tsankov_sarc_biph,
                                         annotation,
                                         tissues_to_consider, 
                                         ML_model) {
  scATAC_source = paste("scATAC_source", scATAC_source, "cell_number_filter", 
                        cell_number_filter, sep="_")
  
  if (tss_fragment_filter != -1) {
    scATAC_source = paste(scATAC_source, "tss_fragment_filter", 
                          tss_fragment_filter, sep="_")
  }
    
  if (waddell_sarc_biph) {
      scATAC_source = paste(scATAC_source, "waddell_sarc_biph", 
                            sep="_")
  } else if (waddell_sarc) {
      scATAC_source = paste(scATAC_source, "waddell_sarc", sep="_")
  } else if (waddell_sarc_tsankov_sarc) {
      scATAC_source = paste(scATAC_source, "waddell_sarc_tsankov_sarc",
                            sep="_")
  } else if (waddell_sarc_biph_tsankov_sarc_biph) {
      scATAC_source = paste(scATAC_source, "waddell_sarc_biph_tsankov_sarc_biph", 
                            sep="_")
  }
  
  scATAC_source = paste(scATAC_source, "annotation", annotation, 
                        sep="_")
  
  dir = paste("/ahg", "regevdata", "projects", "ICA_Lung", "Mohamad", "cell_of_origin",
              "figures", "models", ML_model, cancer_type, scATAC_source,
              "backwards_elimination_results", sep="/")
  if (tissues_to_consider != "all") {
    dir = paste(dir, tissues_to_consider, sep="_")
  }

  return(dir)
}

get_relevant_backwards_elim_dirs <- function(args) {
    cancer_types = args$cancer_types
    cancer_types = unlist(strsplit(cancer_types, split = ","))
    combined_datasets = args$combined_datasets
    tissues_to_consider = paste(unlist(strsplit(args$tissues_to_consider, 
                                                split=","), 
                                       "_"))
    # bing_ren = args$bing_ren
    # shendure = args$shendure
    # tsankov = args$tsankov
    datasets = unlist(strsplit(args$datasets, split = ","))
    cell_number_filter = args$cell_number_filter
    tss_fragment_filter = unlist(strsplit(args$tss_fragment_filter, split = ","))
    waddell_sarc_biph = args$waddell_sarc_biph
    waddell_sarc = args$waddell_sarc
    annotation = args$annotation
    waddell_sarc_tsankov_sarc = args$waddell_sarc_tsankov_sarc
    waddell_sarc_biph_tsankov_sarc_biph = args$waddell_sarc_biph_tsankov_sarc_biph
    ML_model = args$ML_model
    # if (bing_ren) {
    #   scATAC_sources = paste(scATAC_sources, "bing_ren", sep="_")
    # }
    # if (shendure) {
    #   scATAC_sources = paste(scATAC_sources, "shendure", sep="_")
    # }
    # if (tsankov) {
    #   scATAC_sources = paste(scATAC_sources, "tsankov", sep="_")
    # }
    # if (bing_ren && shendure && tsankov) {
    #   scATAC_sources = "combined_datasets"
    # }
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
    
    # scATAC_sources = sub("^_", "", scATAC_sources)
    
    backward_elim_dirs = list()
    for (cancer_type in cancer_types) {
      for (tss_filter in tss_fragment_filter) {
        backward_elim_dirs = append(backward_elim_dirs,
                                    construct_backwards_elim_dir(cancer_type, 
                                                                 scATAC_sources,
                                                                 cell_number_filter,
                                                                 tss_filter,
                                                                 waddell_sarc_biph,
                                                                 waddell_sarc,
                                                                 waddell_sarc_tsankov_sarc,
                                                                 waddell_sarc_biph_tsankov_sarc_biph,
                                                                 annotation,
                                                                 tissues_to_consider,
                                                                 ML_model))
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
                  ggtitle(unlist(strsplit(dir, split ="/"))[11]) +
                  # geom_text(aes(y = ypos, label = features), color = "white", size=1) +
                  theme(plot.title = element_text(hjust = 0.5))
    ggsave(paste(dir, "pie_chart.png", sep="/"), width = 20,
           height = 6, plot)
  }
}

construct_bar_plots <- function(args) {
  dirs = get_relevant_backwards_elim_dirs(args)
  for (dir in dirs) {
    file = paste(dir, "df_for_feature_importance_plots.csv", sep="/")
    df = as_tibble(read.csv(file))
    # df = df %>% filter(num_features == args$bar_plot_num_features)
    df$num_features_f = factor(df$num_features, levels=unique(df$num_features))
    colors = get_n_colors(20, 1)
    
    from = as.character(unique(df$num_features))
    to = paste(paste(levels(df$num_features_f), "features"),
               paste("(R^2=", as.character(round(unique(df$score*100), 1)), 
                     ")", sep=""), sep=" ")
    names(to) <- from
    # df = df %>% group_by(num_features_f) %>% arrange(-importance, .by_group=T)
    print(paste("Plotting data in", dir))
    plot = ggplot(df, aes(x=reorder_within(features, -importance, within=num_features_f,
                                           sep="."), 
                          y=importance, fill=features)) +
          facet_wrap(~num_features_f, nrow=1, 
                     labeller = as_labeller(to), scales = "free") +
           geom_bar(stat="identity", width=1, color="white") +
           xlab("Cell type") +
           ylab("Percent importance (%)") +
           # theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1),
           #       aspect.ratio = 1.1/1) +
           theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1,
                                            size=15)) +
           guides(fill="none") +
           scale_fill_manual(values=colors) +
           # ggtitle(paste0(unlist(strsplit(dir, split ="/"))[3], " (R^2=",
           #               as.character(round(unique(df$score*100), 1)), ")")) +
           ggtitle(unlist(strsplit(dir, split ="/"))[11]) +
           theme(plot.title = element_text(hjust = 0.5))
    ggsave(paste(dir, "bar_plot.png", sep="/"), width = 20, height = 15, plot)
  }
}

cancer_types = args$cancer_types
cancer_types = unlist(strsplit(cancer_types, split = ","))
iters_dont_skip = args$iters_dont_skip
iters_dont_skip = as.integer(unlist(strsplit(iters_dont_skip, split = ",")))
datasets = unlist(strsplit(args$datasets, split = ","))
cell_number_filter = args$cell_number_filter
tss_fragment_filter = unlist(strsplit(args$tss_fragment_filter, split = ","))
waddell_sarc_biph = args$waddell_sarc_biph
waddell_sarc = args$waddell_sarc
num_iter_skips = args$num_iter_skips
annotation = args$annotation
waddell_sarc_tsankov_sarc = args$waddell_sarc_tsankov_sarc
waddell_sarc_biph_tsankov_sarc_biph = args$waddell_sarc_biph_tsankov_sarc_biph
tissues_to_consider = strsplit(args$tissues_to_consider,  split=",")
ML_model = args$ML_model

prep_dfs_command = paste("python3 ../../data/scripts/prep_dfs_for_feature_importance_plots.py", 
                         "--cancer_types", paste(cancer_types, collapse=" "), "--datasets", 
                         paste(datasets, collapse=" "), 
                         "--annotation", annotation, "--tissues_to_consider",
                         paste(tissues_to_consider, collapse = " "), 
                         "--cell_number_filter", cell_number_filter, 
                         "--num_iter_skips", num_iter_skips, "--iters_dont_skip",
                         iters_dont_skip, "--tss_fragment_filter", 
                         paste(tss_fragment_filter, collapse = " "),
                         "--ML_model", ML_model)

if (waddell_sarc_biph) {
  prep_dfs_command = paste(prep_dfs_command, "--waddell_sarc_biph")
} else if (waddell_sarc) {
  prep_dfs_command = paste(prep_dfs_command, "--waddell_sarc")
} else if (waddell_sarc_tsankov_sarc) {
  prep_dfs_command = paste(prep_dfs_command, "--waddell_sarc_tsankov_sarc")
} else if (waddell_sarc_biph_tsankov_sarc_biph) {
  prep_dfs_command = paste(prep_dfs_command, "--waddell_sarc_biph_tsankov_sarc_biph")
}

system(prep_dfs_command)

tissues_to_consider = paste(unlist(tissues_to_consider, "_"))

if (args$pie_chart) {
  construct_pie_charts(args)
} else {
  construct_bar_plots(args)
}

# df_temp <- df %>%
#             filter(num_features == 5) %>%
#             select(-c(num_features, num_features_f))
# 
# df_temp <- df_temp %>% 
#   arrange(desc(features)) %>%
#   mutate(prop = importance / sum(df_temp$importance) *100) %>%
#   mutate(ypos = cumsum(prop)- 0.5*prop )
# 
# ggplot(df_temp, aes(x="", y=prop, fill=features)) +
#   geom_bar(stat="identity", width=1, color="white") +
#   coord_polar("y", start=0) +
#   theme_void() + 
#   theme(legend.position="none") +
#   geom_text(aes(y = ypos, label = features), color = "white", size=2.5) +
#   scale_fill_brewer(palette="Set1")


