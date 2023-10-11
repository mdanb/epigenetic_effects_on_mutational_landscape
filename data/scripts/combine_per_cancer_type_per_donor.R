library(optparse)

option_list <- list( 
  make_option("--cohorts", type="character"),
  make_option("--cancer_type", type="character")
)

args = parse_args(OptionParser(option_list=option_list))
cancer_cohorts = unlist(strsplit(args$cohorts, split = ","))
cancer_type = args$cancer_type
chr_keep = read.csv("../processed_data/chr_keep.csv", row.names = 1)
mutations = data.frame(row.names=chr_keep[["chr"]]) 

for (cohort in cancer_cohorts) {
  for (file in list.files(paste("../processed_data/per_patient_mutations",
                                cohort,
                                sep="/"), 
                          full.names=T)) {
    print(paste("adding", file))
    fname = basename(file)
    patient_id = unlist(lapply(strsplit(fname, split="_"), "[", 3))
    patient_id = unlist(lapply(strsplit(patient_id, split="\\."), "[", 1))
    patient_mutations = read.csv(file, row.names = 1)
    colnames(patient_mutations) = patient_id
    mutations = cbind(mutations, patient_mutations)
  }
}
fname = paste(cancer_type, "per_donor.csv", sep="_")
write.csv(mutations, 
          paste("../processed_data/per_patient_mutations", fname, sep="/"))

