library(GenomicRanges)
library(optparse)
library(tibble)
library(dplyr)
library(magrittr)
library(parallel)

option_list <- list( 
  make_option(c("--cancer_cohorts", type="str"))
)

cancer_cohorts = parse_args(OptionParser(option_list=option_list))
cancer_cohorts = unlist(strsplit(cancer_cohorts$cancer_cohorts, split = ","))

count_overlaps_and_save <- function(df, interval_ranges, cohort) {
  # filename = unlist(strsplit(filepath, "/"))
  # filename = paste("binned_mutations", filename[length(filename)], sep="_")
  # mutations = read.csv(filepath)
  # if (!file.exists(paste(per_patient_mutations_dir, 
    #                      filename, sep="/"))) {
    # if (!grepl("chr", mutations["chromosome"])) {
    #   seqnames = paste(rep("chr", dim(mutations["chromosome"])[1]), 
    #                    mutations[["chromosome"]], sep="")
    # }
    # else {
    #   seqnames = mutations[["chromosome"]]
    # }
    # if ("chrMT" %in% seqnames) {
    #   idx = match("chrMT", seqnames)
    #   mutations = mutations[-idx, ]
    #   seqnames = seqnames[-idx]
    # }
    if (!grepl("chr", df["Chromosome"])) {
      seqnames = paste(rep("chr", dim(df["Chromosome"])[1]),
                       df[["Chromosome"]], sep="")
    } else {
      seqnames = df[["Chromosome"]]
    }
    irange_patient_mutation = IRanges(start=
                                        as.vector(unlist(df["Start_position"])), 
                                      end=
                                        as.vector(unlist(df["End_position"])))
    grange_patient_mutation = GRanges(seqnames=seqnames,
                                      ranges=irange_patient_mutation)
    count_overlaps = as.data.frame(countOverlaps(interval_ranges, 
                                                 grange_patient_mutation))
    regions_to_use = read.csv("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/mutation_data/regions_to_use.csv")["x"]
    count_overlaps = count_overlaps[rownames(count_overlaps) %in% regions_to_use[["x"]], ]
    # colnames(count_overlaps) <- "num_mutations"
    filename = paste0("binned_mutations_", unique(df["Donor_ID"]), ".csv")
    write.csv(count_overlaps, paste("..", "processed_data",
                                             "per_patient_mutations",
                                             cohort,
                                             filename,
                                             sep="/"))
}

count_overlaps_and_save_helper <- function(cancer_cohort) {
  root = "/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/mutation_data/icgc/per_cohort_mutations/"
  cohort_file = paste0(cancer_cohort, "_SNV_without_SEX.txt")
  cohort_df = as_tibble(read.csv(paste(root, cohort_file, sep="/"), sep="\t"))
  header = colnames(read.csv(paste(root, "maf_header.txt", sep="/"), sep="\t"))
  colnames(cohort_df) = header
  # filespath = paste("processed_data/per_patient_mutations", cancer_cohort, sep="/")
  # files = grep("^((?!binned).)*$", 
  #              list.files(filespath, full.names = TRUE),
  #              perl = TRUE, 
  #              value = TRUE)
  dir_donors_with_WGS = "/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/mutation_data/icgc/MutPerBarcodeGroup"
  donors_with_WGS = read.csv(paste(dir_donors_with_WGS, paste0(cancer_cohort, "_SNV.txt"),
                                   sep="/"), sep="\t")
  cohort_df = cohort_df %>% filter(Donor_ID %in% donors_with_WGS[["DonorID"]])
  cohort_df = cohort_df %>%
                group_by(Donor_ID) %>%
                group_split()
  dir.create(paste("../processed_data/per_patient_mutations", cancer_cohort, sep="/"))
  load('/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/mutation_data/hg19.1Mb.ranges.Polak.Nature2015.RData')
  lapply(cohort_df, count_overlaps_and_save, interval.ranges, cancer_cohort)
  # lapply(files, count_overlaps_and_save, filespath)
}
# Lung_AdenoCA = list.files("processed_data/per_patient_mutations/Lung-AdenoCA", 
#                           full.names = TRUE)
# Skin_Melanoma = list.files("processed_data/per_patient_mutations/Skin-Melanoma", 
#                            full.names = TRUE)

lapply(cancer_cohorts, count_overlaps_and_save_helper)

# for (cancer_type in cancer_types) {
#   filespath = paste("processed_data/per_patient_mutations", cancer_type, sep="/")
#   files = list.files(filespath, full.names = TRUE)
#   lapply(files, count_overlaps_and_save, filespath)
# }

# lapply(Lung_AdenoCA, count_overlaps_and_save, 
#        "processed_data/per_patient_mutations/Lung-AdenoCA")
# lapply(Skin_Melanoma, count_overlaps_and_save, 
#        "processed_data/per_patient_mutations/Skin-Melanoma")
