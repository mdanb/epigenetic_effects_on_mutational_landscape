library(GenomicRanges)
library(optparse)

load('raw_dir/mutation_data/hg19.1Mb.ranges.Polak.Nature2015.RData')

option_list <- list( 
  make_option(c("--cancer_types", type="str"))
)

cancer_types = parse_args(OptionParser(option_list=option_list))
cancer_types = unlist(strsplit(cancer_types$cancer_types, split = ","))

count_overlaps_and_save <- function(filepath, per_patient_mutations_dir) {
  filename = unlist(strsplit(filepath, "/"))
  filename = paste("binned_mutations", filename[length(filename)], sep="_")
  mutations = read.csv(filepath)
  if (!file.exists(paste(per_patient_mutations_dir, 
                         filename, sep="/"))) {
    if (!grepl("chr", mutations["chromosome"])) {
      seqnames = paste(rep("chr", dim(mutations["chromosome"])[1]), 
                       mutations[["chromosome"]], sep="")
    }
    else {
      seqnames = mutations[["chromosome"]]
    }
    if ("chrMT" %in% seqnames) {
      idx = match("chrMT", seqnames)
      mutations = mutations[-idx, ]
      seqnames = seqnames[-idx]
    }
    irange_patient_mutation = IRanges(start=
                                        as.vector(unlist(mutations["chromosome_start"])), 
                                      end=
                                        as.vector(unlist(mutations["chromosome_end"])))
    grange_patient_mutation = GRanges(seqnames=seqnames,
                                      ranges=irange_patient_mutation)
    count_overlaps = as.data.frame(countOverlaps(interval.ranges, 
                                                 grange_patient_mutation))
    colnames(count_overlaps) <- "num_mutations"
    write.csv(as.data.frame(count_overlaps), paste(per_patient_mutations_dir, 
                                                   filename, sep="/"))
  }
}

parallel_count_overlaps_and_save <- function(cancer_type) {
  filespath = paste("processed_data/per_patient_mutations", cancer_type, sep="/")
  files = grep("^((?!binned).)*$", list.files(filespath, full.names = TRUE),
               perl = TRUE, value = TRUE)
  lapply(files, count_overlaps_and_save, filespath)
}
# Lung_AdenoCA = list.files("processed_data/per_patient_mutations/Lung-AdenoCA", 
#                           full.names = TRUE)
# Skin_Melanoma = list.files("processed_data/per_patient_mutations/Skin-Melanoma", 
#                            full.names = TRUE)

lapply(cancer_types, parallel_count_overlaps_and_save)

# for (cancer_type in cancer_types) {
#   filespath = paste("processed_data/per_patient_mutations", cancer_type, sep="/")
#   files = list.files(filespath, full.names = TRUE)
#   lapply(files, count_overlaps_and_save, filespath)
# }

# lapply(Lung_AdenoCA, count_overlaps_and_save, 
#        "processed_data/per_patient_mutations/Lung-AdenoCA")
# lapply(Skin_Melanoma, count_overlaps_and_save, 
#        "processed_data/per_patient_mutations/Skin-Melanoma")
