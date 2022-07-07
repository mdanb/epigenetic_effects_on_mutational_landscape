library(GenomicRanges)

Lung_AdenoCA = list.files("processed_data/per_patient_mutations/Lung-AdenoCA", 
                           full.names = TRUE)
Skin_Melanoma = list.files("processed_data/per_patient_mutations/Skin-Melanoma", 
                          full.names = TRUE)
load('raw_dir/mutation_data/hg19.1Mb.ranges.Polak.Nature2015.RData')

count_overlaps_and_save <- function(filepath, per_patient_mutations_dir) {
  mutations = read.csv(filepath)
  if (!grepl("chr", mutations["chromosome"])) {
    seqnames = paste(rep("chr", dim(mutations["chromosome"])[1]), 
                     mutations[["chromosome"]], sep="")
  }
  else {
    seqnames = mutations[["chromosome"]]
  }
  if ("chrMT" %in% mutations[["chromosome"]]) {
    idx = match("chrMT", mutations[["chromosome"]])
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
  filename = unlist(strsplit(filepath, "/"))
  filename = paste("binned_mutations", filename[length(filename)], sep="_")
  write.csv(as.data.frame(count_overlaps), paste(per_patient_mutations_dir, 
                                                 filename, sep="/"))
}

lapply(Lung_AdenoCA, count_overlaps_and_save, 
       "processed_data/per_patient_mutations/Lung-AdenoCA")
lapply(Skin_Melanoma, count_overlaps_and_save, 
       "processed_data/per_patient_mutations/Skin-Melanoma")
