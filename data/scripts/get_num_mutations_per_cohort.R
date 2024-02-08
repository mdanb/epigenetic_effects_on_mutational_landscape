library(dplyr)
library(tibble)

chr_keep = read.csv("../processed_data/chr_keep.csv", row.names = 1)[["chr"]]

df = read.table("../mutation_data/PCAWG_Cancergroup_MutationCountperMB.txt",
                header=T)
df = df[, 5:ncol(df)]
df = apply(df, 2, sum)
names(df) = gsub("\\.", "-", names(df))
df = tibble(CancerType=names(df), mut_counts=df)
sample_stats = read.table("../processed_data/Sample_Statistics_of_PCAWG_Really_Used.tsv",
                          header=T)
sample_stats["CancerType"] = unlist(lapply(strsplit(
                                    sample_stats[["CancerType"]], "_"), "[", 1))
sample_stats = left_join(sample_stats, df)
sample_stats = sample_stats %>% filter(!(CancerType %in% c("Adenocarcinoma",
                                                           "Bone-Cart",
                                                           "Bone-Epith",
                                                           "Breast-DCIS", 
                                                           "Breast-LobularCa",
                                                           "Breast", 
                                                           "CNS", 
                                                           "Carcinoma",
                                                           "Cervix-AdenoCA",
                                                           "Digestive",
                                                           "Female",
                                                           "Glioma",
                                                           "Hematopoietic",
                                                           "Kidney",
                                                           "Lung",
                                                           "Lymph-NOS",
                                                           "Lymph",
                                                           "Myeloid-MDS",
                                                           "Myeloid",
                                                           "Pancan",
                                                           "Sarcoma",
                                                           "Squamous",
                                                           "Kidney-RCC")))

meso = read.csv("../processed_data/mesothelioma.csv")[["waddell_combined"]]
meso = sum(meso, na.rm=T)
sample_stats = sample_stats %>% add_row(CancerType="mesothelioma", 
                                        NofDonor=58,
                                        mut_counts=meso)

sclc = read.csv("../processed_data/sclc_count_overlaps.csv", row.names=1)
sclc = sclc[chr_keep, ]
sclc = sum(sclc)

sample_stats = sample_stats %>% add_row(CancerType="SCLC", 
                                        NofDonor=109,
                                        mut_counts=sclc)

df = read.csv("../processed_data/histologically_subtyped_mutations.csv", 
               row.names=1)
pRCC = sum(df["kidney_papillary"])
ccRCC = sum(df["kidney_clear_cell"])
kidney = read.csv("../processed_data/mutations_with_subtypes/kidney_all.csv")
table(kidney[["subtype"]])
sample_stats = sample_stats %>% add_row(CancerType="PRCC", 
                                        NofDonor=32,
                                        mut_counts=pRCC)
sample_stats = sample_stats %>% add_row(CancerType="ccRCC", 
                                        NofDonor=111,
                                        mut_counts=ccRCC)

msi_high = read.csv("../processed_data/msi_high.csv", row.names=1)
msi_high = sum(msi_high[chr_keep, ])
sample_stats = sample_stats %>% add_row(CancerType="msi_high", 
                                        NofDonor=7,
                                        mut_counts=msi_high)

mm = read.csv("../processed_data/multiple_myeloma.csv", row.names=1)
mm = sum(mm)
sample_stats = sample_stats %>% add_row(CancerType="multiple_myeloma", 
                                        NofDonor=NA,
                                        mut_counts=mm)

sample_stats["avg_mut_per_sample"] = sample_stats["mut_counts"] / sample_stats["NofDonor"]
sample_stats = sample_stats %>% arrange(desc(avg_mut_per_sample))
write.csv(sample_stats, "../processed_data/mutation_and_sample_counts.csv")