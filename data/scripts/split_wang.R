library(rtracklayer)
library(magrittr)
library(dplyr)

data = import("../bed_files/wang_adult_lung/lung_snATAC.200bp_reads.bed.gz", format="bed")

data$sample = unlist(lapply(strsplit(data$name, "_"), "[", 1))
export(data, 
       "../bed_files/wang_adult_lung/lung_snATAC_with_sample.bed.gz",
       format="bed")

data = as_tibble(data)

data <- data %>%
        group_by(sample) %>%
        group_split()

for (i in seq_along(data)) {
  sample_name = data[[i]][[1, "sample"]]
  export(data[[i]], paste("../bed_files/wang_adult_lung", 
                          paste0(sample_name, ".bed.gz"), sep="/"), 
                          format="bed")
}