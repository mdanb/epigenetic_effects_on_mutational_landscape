library(ggplot2)
library(magrittr)
library(RColorBrewer)
library(dplyr)
library(optparse)
library(stringr)
source("utils.R")

parser <- OptionParser()
parser <- add_option(parser, c("--all_cells"), action="store_true", default=T)
parser <- add_option(parser, c("--combined_datasets"), action="store_true",
                               default=F)
parser <- add_option(parser, c("--bing_ren"), action="store_true",
                               default=F)
parser <- add_option(parser, c("--tsankov"), action="store_true",
                     default=F)
parser <- add_option(parser, c("--shendure"), action="store_true",
                     default=F)
parser <- add_option(parser, c("--cancer_types"), type="character")
parser <- add_option(parser, c("--cell_number_filter"), type="integer")
parser <- add_option(parser, c("--pie_chart"), action="store_true", default=F)
parser <- add_option(parser, c("--bar_plot_num_features"), type="integer")
parser <- add_option(parser, c("--tss_fragment_filter"), 
                     type="integer", default=-1)
group.add_argument('--meso_waddell_and_biphasic', action="store_true",
                   default=False)
group.add_argument('--meso_waddell_only', action="store_true", default=False)
group.add_argument('--meso_waddell_and_broad_only', action="store_true", default=False)
group.add_argument('--meso_waddell_biph_786_846', action="store_true", default=False)

args = parse_args(parser)
args = parse_args(parser, args = c("--cancer_types=Lung-AdenoCA,Lung-SCC",
                                   "--all_cells", "--cell_number_filter=1",
                                   "--bar_plot_num_features=20", "--tsankov"))

construct_backwards_elim_dir <- function(cancer_type, scATAC_source, 
                                         cell_number_filter,
                                         tss_fragment_filter, 
                                         meso_waddell_and_biphasic,
                                         meso_waddell_only,
                                         meso_waddell_and_broad_only,
                                         meso_waddell_biph_786_846) {
  scATAC_source = paste("scATAC_source", scATAC_source, "cell_number_filter", 
                        cell_number_filter, sep="_")
  
  if (tss_fragment_filter != -1) {
    scATAC_source = paste(scATAC_source, "tss_fragment_filter", 
                          tss_fragment_filter, sep="_")
  }
    
  if (meso_waddell_and_biphasic) {
      scATAC_source = paste(scATAC_source, "meso_waddell_and_biphasic", 
                            sep="_")
  } else if (meso_waddell_only) {
      scATAC_source = paste(scATAC_source, "meso_waddell_only", sep="_")
  } else if (meso_waddell_and_broad_only) {
      scATAC_source = paste(scATAC_source, "meso_waddell_and_broad_only",
                            sep="_")
  } else if (meso_waddell_biph_786_846) {
      scATAC_source = paste(scATAC_source, "meso_waddell_biph_786_846", 
                            sep="_")
  }
  
  dir = paste("figures", "models", cancer_type, scATAC_source,
              "backwards_elimination_results", sep="/")
  return(dir)
}

get_relevant_backwards_elim_dirs <- function(args) {
    cancer_types = args$cancer_types
    cancer_types = unlist(strsplit(cancer_types, split = ","))
    all_cells = args$all_cells
    combined_datasets = args$combined_datasets
    bing_ren = args$bing_ren
    tsankov = args$tsankov
    cell_number_filter = args$cell_number_filter
    tss_filtered = args$tss_filtered
    tss_filtered_num_fragment_filter = args$tss_filtered_num_fragment_filter
    meso = args$meso
    
    backward_elim_dirs = list()
    for (cancer_type in cancer_types) {
      if (all_cells) {
        if (bing_ren) {
          backward_elim_dirs = append(backward_elim_dirs,
                        construct_backwards_elim_dir(cancer_type, "bing_ren", 
                                                     cell_number_filter,
                                                     tss_filtered,
                                                     tss_filtered_num_fragment_filter,
                                                     str_to_title(meso)))
        }
        if (tsankov) {
          backward_elim_dirs = append(backward_elim_dirs,
                                      construct_backwards_elim_dir(cancer_type, "tsankov", 
                                                                   cell_number_filter,
                                                                   tss_filtered,
                                                                   tss_filtered_num_fragment_filter,
                                                                   str_to_title(meso)))
        }
        # if (shendure) {
        #   backward_elim_dirs.append(construct_backwards_elim_dir(cancer_type, 
        #                                                          "shendure"))
        # }
        if (combined_datasets) {
          backward_elim_dirs = append(backward_elim_dirs,
                 construct_backwards_elim_dir(cancer_type, "combined_datasets",
                                              cell_number_filter,
                                              tss_filtered,
                                              tss_filtered_num_fragment_filter,
                                              str_to_title(meso)))
        }
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
                  ggtitle(unlist(strsplit(dir, split ="/"))[3]) +
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
    df = df %>% filter(num_features == args$bar_plot_num_features)
    df$num_features_f = factor(df$num_features, levels=unique(df$num_features))
    colors = get_n_colors(20, 1)
    plot = ggplot(df, aes(x=reorder(features, -importance), 
                          y=importance, fill=features)) +
           geom_bar(stat="identity", width=1, color="white") +
           xlab("Cell type") +
           ylab("Percent importance (%)") +
           # theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1),
           #       aspect.ratio = 1.1/1) +
           theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1)) +
           guides(fill="none") +
           scale_fill_manual(values=colors) +
           ggtitle(paste0(unlist(strsplit(dir, split ="/"))[3], " (R^2=",
                         as.character(round(unique(df$score*100), 1)), ")")) +
           theme(plot.title = element_text(hjust = 0.5))
    ggsave(paste(dir, "bar_plot.png", sep="/"), width = 6, height = 8, plot)
  }
}
# labeller = as_labeller(c("20" = "A",
#                          "15" = "B",
#                          "10" = "C",
#                          "5" = "D",
#                          "2" = "E"))) +
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


