library(dplyr)
library(magrittr)
library(tibble)

filter_by_WGS_availability <- function(specimen_file, WGS_patients_file) {
  df = as_tibble(read.csv(specimen_file, sep="\t"))
  WGS_donor_df = read.csv(WGS_patients_file, sep="\t")
  df = df %>% 
          filter(icgc_donor_id %in% WGS_donor_df[["DonorID"]])
  return(df)
} 

custom_map_histologies <- function(df, from, to) {
  histology_map = mapvalues(df["tumour_histological_type"], 
                                         from=from, 
                                         to=to)
  
  df["tumour_histological_type"] = histology_map
  return(df)
}

get_combined_WGS_cancer_df <- function(cohorts) {
  root = "/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/mutation_data/icgc"
  cohort_dfs = list()
  for (cohort in cohorts) {
    curr_df = filter_by_WGS_availability(
                              paste(root, paste("specimen", cohort, "tsv", "gz", 
                                                sep="."), sep="/"),
                              paste(root, "MutPerBarcodeGroup", paste(cohort, "SNV.txt", sep="_"), sep="/")
                              )
    cohort_dfs = append(cohort_dfs, list(curr_df))
  }
  df = bind_rows(cohort_dfs)
  return(df)
}

                               # Histologies #
                      ##################################
# PANCREAS
# 8500/3: DUCT CARCINOMA / Invasive carcinoma of no special type
# 8000/3: NEOPLASM / Neoplasm, malignant
# 8140/3: ADENOCARCINOMA, NOS / Adenocarcinoma, NOS

# BRAIN
# 9421/1: FIBRILLARY ASTROCYTOMA / Pilocytic astrocytoma
# 9470/3: MEDULLOBLASTOMA, NOS / Medulloblastoma, NOS
# 9471/3: MEDULLOBLASTOMA, NOS / Desmoplastic medulloblastoma
# 9474/3: MEDULLOBLASTOMA, NOS / Large cell medulloblastoma

                      ##################################

### PancAdenoCA ###
Panc_AdenoCA_WGS = get_combined_WGS_cancer_df(c("PACA-AU", "PACA-CA"))
Panc_AdenoCA_WGS = custom_map_histologies(Panc_AdenoCA_WGS,
                                          c("", 
                                            "Pancreatic Ductal Adenocarcinoma",
                                            "Intraductal Papillary Mucinous Neoplasm with invasion",
                                            "Mucinous Non-cystic carcinoma",
                                            "Acinar Cell Carcinoma",
                                            "Undifferentiated (anaplastic) carcinoma",
                                            "Undifferentiated carcinoma with osteoclast like giant cells",
                                            "8500/3",
                                            "Mar-00",
                                            ""),
                                          c("Adenocarcinoma", 
                                            "Pancreatic ductal carcinoma",
                                            "Invasive carcinoma arising in IPMN",
                                            "Mucinous adenocarcinoma",
                                            "Acinar cell carcinoma",
                                            "",
                                            "",
                                            "Pancreatic ductal carcinoma",
                                            "Adenocarcinoma",
                                            ""))

Panc_AdenoCA_grouped = Panc_AdenoCA %>% 
                            group_by(icgc_donor_id) %>% 
                            add_count()

# No such cases 
temp = Panc_AdenoCA_grouped %>%
          filter(n == 1, tumour_histological_type == "")

# One patient, two samples, both same histological label
temp = Panc_AdenoCA_grouped %>%                  
  filter(tumour_histological_type == "Undifferentiated (anaplastic) carcinoma")

# One patient, one sample
temp = Panc_AdenoCA_grouped %>%                  
  filter(tumour_histological_type == "Undifferentiated carcinoma with osteoclast like giant cells")

# No such cases
temp = Panc_AdenoCA_grouped %>%
  filter(n == 1, tumour_histological_type == "Mar-00")

# one patient, one sample
temp = Panc_AdenoCA_grouped %>%
  filter(tumour_histological_type == "Mar-60")

### Medulloblastoma ###
PBCA_DE_WGS = get_combined_WGS_cancer_df(c("PBCA-DE"))
PBCA_DE_WGS_grouped = PBCA_DE_WGS %>% 
                          group_by(icgc_donor_id) %>% 
                          add_count()

# No such cases
temp = PBCA_DE_WGS_grouped %>%
              filter(n == 1, tumour_histological_type == "N/A")
PBCA_DE_WGS = PBCA_DE_WGS %>%

