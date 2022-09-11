library(ggplot2)
library(magrittr)
library(RColorBrewer)
library(dplyr)
library(optparse)
source("utils.R")

parser <- OptionParser()
parser <- add_option(parser, c("--all_cells"), action="store_true", default=T)
parser <- add_option(parser, c("--combined_datasets"), action="store_true",
                               default=F)
parser <- add_option(parser, c("--bing_ren"), action="store_true",
                               default=F)
parser <- add_option(parser, c("--cancer_types"), type="character")

args = parse_args(parser, args = c("--cancer_types=Breast-AdenoCA",
                                   "--all_cells", "--bing_ren"))

construct_backwards_elim_dir <- function(cancer_type, scATAC_source) {
  scATAC_source = paste("scATAC_source", scATAC_source, sep="_")
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
    
    backward_elim_dirs = list()
    for (cancer_type in cancer_types) {
      if (all_cells) {
        if (bing_ren) {
          backward_elim_dirs = append(backward_elim_dirs,
                        construct_backwards_elim_dir(cancer_type, "bing_ren"))
        }
        # if (shendure) {
        #   backward_elim_dirs.append(construct_backwards_elim_dir(cancer_type, 
        #                                                          "shendure"))
        # }
        if (combined_datasets) {
          backward_elim_dirs = append(backward_elim_dirs,
                 construct_backwards_elim_dir(cancer_type, "combined_datasets"))
        }
      }
    return(unlist(backward_elim_dirs))
    }
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
    pie_chart_file = paste(dir, "df_for_pie_charts.csv", sep="/")
    df = read.csv(pie_chart_file)
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
  }
  ggsave(paste(dir, "pie_chart.png", sep="/"), width = 20,
         height = 6, plot)
}

# labeller = as_labeller(c("20" = "A",
#                          "15" = "B",
#                          "10" = "C",
#                          "5" = "D",
#                          "2" = "E"))) +
construct_pie_charts(args)

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


