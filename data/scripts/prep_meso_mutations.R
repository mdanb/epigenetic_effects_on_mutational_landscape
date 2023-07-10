library(tidyverse)
library(readxl)
library(GenomicRanges)

load('../mutation_data/hg19.1Mb.ranges.Polak.Nature2015.RData')
rownames = names(interval.ranges)
biphasic_mesomics = read.csv("../mutation_data/Mesothelioma_MMB_MutCountperBin.txt",
                             sep="\t")
sarcomatoid_mesomics = read.csv("../mutation_data/Mesothelioma_MMS_MutCountperBin.txt",
                             sep="\t")
epithelioid_mesomics = read.csv("../mutation_data/Mesothelioma_MME_MutCountperBin.txt",
                                sep="\t")

all_mesomics = cbind(biphasic_mesomics, 
                     sarcomatoid_mesomics[, 4:ncol(sarcomatoid_mesomics)], 
                     epithelioid_mesomics[, 4:ncol(epithelioid_mesomics)])
                     
aggregate_mutations <- function(mut_df, name) {
  agg = rowSums(mut_df[, 4:ncol(mut_df)])
  agg = data.frame(cancer_type = agg)
  colnames(agg) = name
  rownames(agg) = rownames
  return(agg)
}

get_top_bottom_n_perc_ids <- function(df, n, direction, rank_name) {
  if (direction == "top") {
    patients = df[df[[rank_name]] >= quantile(df[[rank_name]], 1 - n), 
                     "Donor"]
  }
  else {
    patients = df[df[[rank_name]] <= quantile(df[[rank_name]], n), 
                  "Donor"]
  }
  print("Sanity Check")
  # sanity check
  print("Top/Bottom")
  print(df[df[["Donor"]] %in% patients, 
                  rank_name])
  patients = gsub("_T.*", "", patients)
  return(patients) 
}

get_top_bottom_n_perc <- function(df, n, rank_name, all_mesomics, name) {
  top_n_perc_ids = get_top_bottom_n_perc_ids(df, n, "top", rank_name)
  bottom_n_perc_ids = get_top_bottom_n_perc_ids(df, n, "bottom", rank_name)
  mutation_data_top_n_perc = all_mesomics[, top_n_perc_ids]
  top_n_perc_agg = aggregate_mutations(mutation_data_top_n_perc, 
                                       paste(name, "top", n * 100, "perc", 
                                             sep = "_"))
  mutation_data_bottom_n_perc = all_mesomics[, bottom_n_perc_ids]
  bottom_n_perc_agg = aggregate_mutations(mutation_data_bottom_n_perc, 
                                          paste(name, "bottom", n * 100, "perc", 
                                          sep = "_"))
  return(list(top_n_perc_agg, bottom_n_perc_agg))
}


agg_epithelioid = aggregate_mutations(epithelioid_mesomics, 
                                      "epithelioid_mesomics")
agg_sarcomatoid = aggregate_mutations(sarcomatoid_mesomics, 
                                      "sarcomatoid_mesomics")
agg_biphasic = aggregate_mutations(biphasic_mesomics, "biphasic_mesomics")

waddell = read.csv("../mutation_data/mesothelioma_WGS_Waddell.csv")

ranked_mesomics = read.table(
                  "../mutation_data/MESOMICS_Sarcomatoid_ModuleScore.txt",
                  header=1)
ranked_mesomics = ranked_mesomics[gsub("_T.*", "", ranked_mesomics[["Donor"]]) 
                                  %in% colnames(all_mesomics), ]

### Blum ###
blum_agg_data_10_perc = get_top_bottom_n_perc(ranked_mesomics, 0.1, 
                                      "Sarcomatoid.Blum.rank",
                                      all_mesomics,
                                      "blum")
blum_top_10_perc = blum_agg_data_10_perc[[1]]
blum_bottom_10_perc = blum_agg_data_10_perc[[2]]

blum_agg_data_15_perc = get_top_bottom_n_perc(ranked_mesomics, 0.15, 
                                              "Sarcomatoid.Blum.rank",
                                              all_mesomics,
                                              "blum")
blum_top_15_perc = blum_agg_data_15_perc[[1]]
blum_bottom_15_perc = blum_agg_data_15_perc[[2]]

blum_agg_data_20_perc = get_top_bottom_n_perc(ranked_mesomics, 0.20, 
                                              "Sarcomatoid.Blum.rank",
                                              all_mesomics,
                                              "blum")
blum_top_20_perc = blum_agg_data_20_perc[[1]]
blum_bottom_20_perc = blum_agg_data_20_perc[[2]]
### NMF ###
nmf_agg_data_10_perc = get_top_bottom_n_perc(ranked_mesomics, 0.1, 
                                      "Sarcomatoid.nmf.old.rank",
                                      all_mesomics,
                                      "nmf")
nmf_top_10_perc = nmf_agg_data_10_perc[[1]]
nmf_bottom_10_perc = nmf_agg_data_10_perc[[2]]

nmf_agg_data_15_perc = get_top_bottom_n_perc(ranked_mesomics, 0.15, 
                                             "Sarcomatoid.nmf.old.rank",
                                             all_mesomics,
                                             "nmf")
nmf_top_15_perc = nmf_agg_data_15_perc[[1]]
nmf_bottom_15_perc = nmf_agg_data_15_perc[[2]]

nmf_agg_data_20_perc = get_top_bottom_n_perc(ranked_mesomics, 0.20, 
                                             "Sarcomatoid.nmf.old.rank",
                                             all_mesomics,
                                             "nmf")
nmf_top_20_perc = nmf_agg_data_20_perc[[1]]
nmf_bottom_20_perc = nmf_agg_data_20_perc[[2]]

### Blum 50 ###
blum_50_agg_data_10_perc = get_top_bottom_n_perc(ranked_mesomics, 0.1, 
                                         "Sarcomatoid.Blum.50.rank",
                                         all_mesomics,
                                         "blum_50")
blum_50_top_10_perc = blum_50_agg_data_10_perc[[1]]
blum_50_bottom_10_perc = blum_50_agg_data_10_perc[[2]]

blum_50_agg_data_15_perc = get_top_bottom_n_perc(ranked_mesomics, 0.15, 
                                                 "Sarcomatoid.Blum.50.rank",
                                                 all_mesomics,
                                                 "blum_50")
blum_50_top_15_perc = blum_50_agg_data_15_perc[[1]]
blum_50_bottom_15_perc = blum_50_agg_data_15_perc[[2]]

blum_50_agg_data_20_perc = get_top_bottom_n_perc(ranked_mesomics, 0.20, 
                                                 "Sarcomatoid.Blum.50.rank",
                                                 all_mesomics,
                                                 "blum_50")
blum_50_top_20_perc = blum_50_agg_data_20_perc[[1]]
blum_50_bottom_20_perc = blum_50_agg_data_20_perc[[2]]



ranked_mesomics_agg = cbind(blum_top_10_perc, blum_bottom_10_perc,
                            nmf_top_10_perc, nmf_bottom_10_perc, 
                            blum_50_top_10_perc, blum_50_bottom_10_perc,
                            blum_top_15_perc, blum_bottom_15_perc,
                            nmf_top_15_perc, nmf_bottom_15_perc,
                            blum_50_top_15_perc, blum_50_bottom_15_perc,
                            blum_top_20_perc, blum_bottom_20_perc,
                            nmf_top_20_perc, nmf_bottom_20_perc,
                            blum_50_top_20_perc, blum_50_bottom_20_perc)

unranked_agg = cbind(agg_epithelioid, agg_sarcomatoid, agg_biphasic)
meso_df = cbind(ranked_mesomics_agg, unranked_agg)
write.csv(meso_df, "../processed_data/mesothelioma.csv")
# nmf_top_10_perc_ids = get_top_bottom_n_perc(ranked_mesomics, 0.1, 
#                                          "top", "Sarcomatoid.nmf.old.rank")
# nmf_bottom_10_perc_ids = get_top_bottom_n_perc(ranked_mesomics, 0.1, 
#                                             "bottom", 
#                                            "Sarcomatoid.nmf.old.rank")
# nmf_mutation_data_top_10_perc = all_mesomics[, nmf_top_10_perc_ids]
# nmf_top_10_perc = aggregate_mutations(nmf_mutation_data_top_10_perc, 
#                                        "nmf_top_10_perc")
# nmf_mutation_data_bottom_10_perc = all_mesomics[, nmf_bottom_10_perc_ids]
# nmf_bottom_10_perc = aggregate_mutations(nmf_mutation_data_bottom_10_perc, 
#                                           "nmf_bottom_10_perc")


# 
# blum_50_top_10_perc_ids = get_top_bottom_n_perc(ranked_mesomics, 0.1, 
#                                         "top", "Sarcomatoid.Blum.50.rank")
# blum_50_bottom_10_perc_ids = get_top_bottom_n_perc(ranked_mesomics, 0.1, 
#                                            "bottom", 
#                                            "Sarcomatoid.Blum.50.rank")
# blum_50_mutation_data_top_10_perc = all_mesomics[, blum_50_top_10_perc_ids]
# blum_50_top_10_perc = aggregate_mutations(blum_50_mutation_data_top_10_perc, 
#                                        "blum_50_top_10_perc")
# blum_50_mutation_data_bottom_10_perc = all_mesomics[, blum_50_bottom_10_perc_ids]
# blum_50_bottom_10_perc = aggregate_mutations(blum_50_mutation_data_bottom_10_perc, 
#                                           "blum_50_bottom_10_perc")


# mut_df = mut_df[, 1:3]

# mut_df["agg"] = agg
# 
# mut_overlaps = as.data.frame(countOverlaps(interval.ranges, mut_granges))



