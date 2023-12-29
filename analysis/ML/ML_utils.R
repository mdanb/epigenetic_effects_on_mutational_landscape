

get_relevant_backwards_elim_dirs <- function(cancer_types, 
                                             # combined_datasets,
                                             tissues_to_consider,
                                             datasets,
                                             cell_number_filter,
                                             tss_fragment_filter,
                                             annotation,
                                             ML_model,
                                             hundred_kb,
                                             fold_for_test_set = "-1",
                                             seed = "-1",
                                             accumulated_seeds=F,
                                             per_donor=F,
                                             cell_types_keep=NULL) {
  if (accumulated_seeds) {
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
                                                               seed,
                                                               hundred_kb,
                                                               fold_for_test_set,
                                                               # test=F,
                                                               per_donor,
                                                               cell_types_keep=cell_types_keep))
    }
  }
  return(unlist(backward_elim_dirs))
}

construct_dir <- function(scATAC_source,
                          cell_number_filter,
                          tss_fragment_filter,
                          annotation,
                          seed,
                          fold_for_test_set,
                          ML_model,
                          cancer_type,
                          hundred_kb,
                          tissues_to_consider,
                          expanded_hundred_kb=F,
                          grid_analysis=F,
                          grid_cell_type=NULL,
                          cell_types_keep=NULL) {
  dir = paste("seed", seed, sep="_")
  
  if (!is.null(cell_types_keep)) {
    dir = paste("ctk", cell_types_keep, dir, sep="_")
  }
  
  dir = paste("annotation", annotation, dir, sep="_")
  if (tss_fragment_filter != -1) {
    dir = paste(dir, "tss_fragment_filter", tss_fragment_filter, sep="_")
  }
  
  if (tissues_to_consider != "all") {
    dir = paste("tissues_to_consider", tissues_to_consider, dir, sep="_")
  }
  
  if (!grid_analysis) {
    dir = paste("cell_number_filter", cell_number_filter, dir, sep="_")
    dir = paste("scATAC_source", scATAC_source, dir, sep="_")
  }

  if (fold_for_test_set != "-1") {
    dir = paste(dir, "fold_for_test_set", fold_for_test_set, sep="_")
  }
  
  if (hundred_kb) {
    dir = paste("interval_ranges_100kb", dir, sep="_")
  } else if (expanded_hundred_kb) {
    dir = paste("expanded_interval_ranges_100kb", dir, sep="_")
  }
  
  return(paste("models", ML_model, cancer_type, dir, sep="/"))
}

construct_backwards_elim_dir <- function(cancer_type, 
                                         scATAC_source, 
                                         cell_number_filter,
                                         tss_fragment_filter, 
                                         annotation,
                                         tissues_to_consider,
                                         ML_model,
                                         seed,
                                         hundred_kb,
                                         fold_for_test_set="-1",
                                         # test=F,
                                         per_donor=F,
                                         grid_analysis=F,
                                         grid_cell_type=NULL,
                                         cell_types_keep=NULL) {
  dir = construct_dir(scATAC_source=scATAC_source,
                      cell_number_filter=cell_number_filter,
                      tss_fragment_filter=tss_fragment_filter,
                      annotation=annotation,
                      seed=seed,
                      fold_for_test_set=fold_for_test_set,
                      ML_model=ML_model,
                      cancer_type=cancer_type,
                      hundred_kb=hundred_kb,
                      tissues_to_consider=tissues_to_consider,
                      grid_analysis=grid_analysis,
                      grid_cell_type=grid_cell_type,
                      cell_types_keep=cell_types_keep)
  
  if (per_donor) {
    part_one = unlist(strsplit(dir, split="/"))[1:2]
    pattern = paste0(cancer_type, "_")
    part_two = unlist(strsplit(dir, split="/"))[4]
    dir = paste(part_one[1:2], collapse="/")
    if (basename(getwd()) == "analysis") {
      dir = paste("ML", dir, sep="/")
    }
    dir = list.files(dir, pattern=pattern)
    dir = paste(paste(part_one, collapse="/"), dir, part_two, sep="/")
  }
  
  dir = paste(dir, "backwards_elimination_results", sep="/")
  
  # if (tissues_to_consider != "all") {
  #   dir = paste(dir, tissues_to_consider, sep="_")
  # }
  
  # if (!test) {
  #   dir = paste("../../figures", dir, sep="/")
  # }
  
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


# construct_pie_charts <- function(args) {
#   dirs = get_relevant_backwards_elim_dirs(args)
#   for (dir in dirs) {
#     file = paste(dir, "df_for_feature_importance_plots.csv", sep="/")
#     df = read.csv(file)
#     df$num_features_f = factor(df$num_features, levels=unique(df$num_features))
#     colors = get_n_colors(20, 1)
#     from = as.character(unique(df$num_features))
#     to = paste(paste(levels(df$num_features_f), "features"),
#                paste("(R^2=", as.character(round(unique(df$score*100), 1)), 
#                      ")", sep=""), sep=" ")
#     
#     names(to) <- from
#     title = unlist(strsplit(dir, split ="/"))
#     title = title[length(title)]
#     print("title: ")
#     print(title)
#     plot = ggplot(df, aes(x="", y=importance, fill=features)) +
#       facet_wrap(~num_features_f, nrow=1, 
#                  labeller = as_labeller(to)) +
#       geom_bar(stat="identity", width=1, color="white") +
#       coord_polar("y", start=0) +
#       xlab("") +
#       ylab("") +
#       theme(axis.text = element_blank(),
#             axis.ticks = element_blank(),
#             panel.grid  = element_blank()) +
#       scale_fill_manual(values=colors) +
#       ggtitle() +
#       theme(plot.title = element_text(hjust = 0.5))
#     ggsave(paste(dir, "pie_chart.png", sep="/"), width = 20,
#            height = 6, plot)
#   }
# }