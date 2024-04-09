library(readxl)
library(dplyr)

bingren = readRDS("../processed_data/count_overlap_data/combined_count_overlaps/finalized_annotation/Bingren_combined_count_overlaps_metadata.rds")
bingren["dataset"] = "D1"
bingren = bingren[bingren["num_cells"] >= 100, ]
bingren = bingren[, colnames(bingren)]

shendure = readRDS("../processed_data/count_overlap_data/combined_count_overlaps/finalized_annotation/Shendure_combined_count_overlaps_metadata.rds")
shendure["dataset"] = "D2"
shendure = shendure[shendure["num_cells"] >= 100, ]
shendure = shendure[, colnames(bingren)]

gl_colon = readRDS("../processed_data/count_overlap_data/combined_count_overlaps/finalized_annotation/Greenleaf_colon_combined_count_overlaps_metadata.rds")
gl_colon["dataset"] = "D3"
gl_colon = gl_colon[gl_colon["num_cells"] >= 100, ]
gl_colon = gl_colon[, colnames(bingren)]

gl_colon_w_polyp =  readRDS("../processed_data/count_overlap_data/combined_count_overlaps/Greenleaf_colon_remove_cancer_merge_normal_unaffected/Greenleaf_colon_combined_count_overlaps_metadata.rds")
gl_colon_w_polyp["dataset"] = "D3"
gl_colon_w_polyp = gl_colon_w_polyp[gl_colon_w_polyp["num_cells"] >= 100, ]
gl_colon_w_polyp = gl_colon_w_polyp[, colnames(bingren)]

gl_blood = readRDS("../processed_data/count_overlap_data/combined_count_overlaps/finalized_annotation/Greenleaf_pbmc_bm_combined_count_overlaps_metadata.rds")
gl_blood["dataset"] = "D4"
gl_blood = gl_blood[gl_blood["num_cells"] >= 100, ]
gl_blood = gl_blood[, colnames(bingren)]

gl_blood_other = readRDS("../processed_data/count_overlap_data/combined_count_overlaps/new_intermediate_blood_bm_annotation/Greenleaf_pbmc_bm_combined_count_overlaps_metadata.rds")
gl_blood_other["dataset"] = "D4"
gl_blood_other = gl_blood_other[gl_blood_other["num_cells"] >= 100, ]
gl_blood_other = gl_blood_other[, colnames(bingren)]

tsankov = readRDS("../processed_data/count_overlap_data/combined_count_overlaps/finalized_annotation/Tsankov_combined_count_overlaps_metadata.rds")
tsankov["dataset"] = "D5"
tsankov = tsankov[tsankov["num_cells"] >= 100 | tsankov["cell_type"] == "Neuroendocrine", ]
tsankov = tsankov[, colnames(bingren)]

yang_kidney = readRDS("../processed_data/count_overlap_data/combined_count_overlaps/finalized_annotation/Yang_kidney_combined_count_overlaps_metadata.rds")
yang_kidney["dataset"] = "D6"
yang_kidney = yang_kidney[yang_kidney["num_cells"] >= 100, ]
yang_kidney = yang_kidney[, colnames(bingren)]

gl_brain = readRDS("../processed_data/count_overlap_data/combined_count_overlaps/finalized_annotation/Greenleaf_brain_combined_count_overlaps_metadata.rds")
gl_brain["dataset"] = "D7"
gl_brain = gl_brain[gl_brain["num_cells"] >= 100, ]
gl_brain = gl_brain[, colnames(bingren)]

bingren_brain = readRDS("../processed_data/count_overlap_data/combined_count_overlaps/finalized_annotation/Bingren_adult_brain_combined_count_overlaps_metadata.rds")
# colnames(bingren_brain)[3] = "tissue_name"
bingren_brain["dataset"] = "D8"
bingren_brain = bingren_brain[bingren_brain["num_cells"] >= 100, ]
bingren_brain = bingren_brain[, colnames(bingren)]

rawlins_fetal_lung = readRDS("../processed_data/count_overlap_data/combined_count_overlaps/finalized_annotation/Rawlins_fetal_lung_combined_count_overlaps_metadata.rds")
rawlins_fetal_lung["dataset"] = "D9"
rawlins_fetal_lung = rawlins_fetal_lung[rawlins_fetal_lung["num_cells"] >= 100, ]
rawlins_fetal_lung = rawlins_fetal_lung[, colnames(bingren)]


with_polyp_and_other_blood = rbind(bingren, shendure, gl_colon_w_polyp, gl_blood_other, 
                                   tsankov, yang_kidney, gl_brain, 
                                   bingren_brain, rawlins_fetal_lung)

with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "ENS glia"
                           & with_polyp_and_other_blood["dataset"] == "D2",
                           "Abbreviation Meaning"] = "Enteric nervous system glia"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "ENS neurons"
                           & with_polyp_and_other_blood["dataset"] == "D2",
                           "Abbreviation Meaning"] = "Enteric nervous system neurons"

with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "TA1"
                           & with_polyp_and_other_blood["dataset"] == "D3",
                           "Abbreviation Meaning"] = "Transit Amplifying 1 Cells"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "TA2"
                           & with_polyp_and_other_blood["dataset"] == "D3",
                           "Abbreviation Meaning"] = "Transit Amplifying 2 Cells"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "Secretory TA"
                           & with_polyp_and_other_blood["dataset"] == "D3",
                           "Abbreviation Meaning"] = "Secretory Transit Amplifying Cells"

with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "NK"
                           & with_polyp_and_other_blood["dataset"] == "D4",
                           "Abbreviation Meaning"] = "Natural Killer Cells"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "cDC"
                           & with_polyp_and_other_blood["dataset"] == "D4",
                           "Abbreviation Meaning"] = "Classical Dendritic Cells"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "CMP.LMPP"
                           & with_polyp_and_other_blood["dataset"] == "D4",
                           "Abbreviation Meaning"] = "Common Myeloid Progenitor / lymphoid-primed multipotent progenitors"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "Early.Baso"
                           & with_polyp_and_other_blood["dataset"] == "D4",
                           "Abbreviation Meaning"] = "Early Basophils"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "Early.Eryth"
                           & with_polyp_and_other_blood["dataset"] == "D4",
                           "Abbreviation Meaning"] = "Early Erythrocytes"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "GMP"
                           & with_polyp_and_other_blood["dataset"] == "D4",
                           "Abbreviation Meaning"] = "Granulocyte-Monocyte progenitors"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "GMP.Neut"
                           & with_polyp_and_other_blood["dataset"] == "D4",
                           "Abbreviation Meaning"] = "Granulocyte-Monocyte progenitors / Neutrophils"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "HSC"
                           & with_polyp_and_other_blood["dataset"] == "D4",
                           "Abbreviation Meaning"] = "Hematopoietic Stem Cells"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "Late.Eryth"
                           & with_polyp_and_other_blood["dataset"] == "D4",
                           "Abbreviation Meaning"] = "Late Erythrocytes"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "Pre.B"
                           & with_polyp_and_other_blood["dataset"] == "D4",
                           "Abbreviation Meaning"] = "Pre B-cell Progenitor"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "pDC"
                           & with_polyp_and_other_blood["dataset"] == "D4",
                           "Abbreviation Meaning"] = "Plasmacytoid Dendritic Cells"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "CLP"
                           & with_polyp_and_other_blood["dataset"] == "D4",
                           "Abbreviation Meaning"] = "Common Lymphoid Progenitor"

with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "AT1"
                           & with_polyp_and_other_blood["dataset"] == "D5",
                           "Abbreviation Meaning"] = "Alveolar type I"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "AT2"
                           & with_polyp_and_other_blood["dataset"] == "D5",
                           "Abbreviation Meaning"] = "Alveolar type II"

with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "ENDO"
                           & with_polyp_and_other_blood["dataset"] == "D6",
                           "Abbreviation Meaning"] = "Endothelial Cells"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "LEUK"
                           & with_polyp_and_other_blood["dataset"] == "D6",
                           "Abbreviation Meaning"] = "Leukocytes"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "LOH"
                           & with_polyp_and_other_blood["dataset"] == "D6",
                           "Abbreviation Meaning"] = "Loop of Henle Cells"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "MES"
                           & with_polyp_and_other_blood["dataset"] == "D6",
                           "Abbreviation Meaning"] = "Mesangial Cells"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "PEC"
                           & with_polyp_and_other_blood["dataset"] == "D6",
                           "Abbreviation Meaning"] = "Glomerular Parietal Epithelial Cells"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "PODO"
                           & with_polyp_and_other_blood["dataset"] == "D6",
                           "Abbreviation Meaning"] = "Podocytes"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "PTPL"
                           & with_polyp_and_other_blood["dataset"] == "D6",
                           "Abbreviation Meaning"] = "Proximal Tubule Progenitor-Like"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "PT"
                           & with_polyp_and_other_blood["dataset"] == "D6",
                           "Abbreviation Meaning"] = "Proximal Tubule Cells"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "ICA"
                           & with_polyp_and_other_blood["dataset"] == "D6",
                           "Abbreviation Meaning"] = "Collecting duct, Intercalated Cells Type A"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "DCT"
                           & with_polyp_and_other_blood["dataset"] == "D6",
                           "Abbreviation Meaning"] = "Distal Convoluted Tubule Cells"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "CD_PC"
                           & with_polyp_and_other_blood["dataset"] == "D6",
                           "Abbreviation Meaning"] = "Collecting Duct Principal Cells"

with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "GluN"
                           & with_polyp_and_other_blood["dataset"] == "D7",
                           "Abbreviation Meaning"] = "Cortical Glutamatergic Neurons"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "IN"
                           & with_polyp_and_other_blood["dataset"] == "D7",
                           "Abbreviation Meaning"] = "GABAergic neurons"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "nIPC"
                           & with_polyp_and_other_blood["dataset"] == "D7",
                           "Abbreviation Meaning"] = "Neuronal Intermedate Progenitor Cell"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "late RG"
                           & with_polyp_and_other_blood["dataset"] == "D7",
                           "Abbreviation Meaning"] = "Late Radial Glia"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "early RG"
                           & with_polyp_and_other_blood["dataset"] == "D7",
                           "Abbreviation Meaning"] = "Early Radial Glia"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "mGPC"
                           & with_polyp_and_other_blood["dataset"] == "D7",
                           "Abbreviation Meaning"] = "Multipotent Glial Progenitor Cells"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "OPC/Oligo"
                           & with_polyp_and_other_blood["dataset"] == "D7",
                           "Abbreviation Meaning"] = "Oligodendrocyte Progenitor Cell/Oligodendrocyte"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "Peric"
                           & with_polyp_and_other_blood["dataset"] == "D7",
                           "Abbreviation Meaning"] = "Pericyte"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "MG"
                           & with_polyp_and_other_blood["dataset"] == "D7",
                           "Abbreviation Meaning"] = "Microglia"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "EC"
                           & with_polyp_and_other_blood["dataset"] == "D7",
                           "Abbreviation Meaning"] = "Endothelial Cell"

# adult_brain_metadata = read_excel("../metadata/adult_brain_metadata.xlsx")
# adult_brain_metadata = unique(adult_brain_metadata[, c("Cell type", 
#                                                        "Cell type name")])
# colnames(adult_brain_metadata) = c("cell_type", "Abbreviation Meaning")
# joined_data <- left_join(with_polyp_and_other_blood, adult_brain_metadata, by = "cell_type")
# joined_data <- joined_data %>%
#   mutate(Abbreviation_Meaning = coalesce(`Abbreviation Meaning.x`, `Abbreviation Meaning.y`)) %>%
#   select(-`Abbreviation Meaning.x`, -`Abbreviation Meaning.y`)

with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "ACBGM"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "Bergmann glia"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "AMY"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "Glutamatergic neurons from amygdala"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "ASCNT"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "Non-telencephalon astrocytes"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "ASCT"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "Telencephalon astrocytes"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "BFEXA"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "Inhibitory neurons from basal forebrain and extended amygdala"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "BNGA"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "GABAergic neurons from basal nuclei"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "CBGRC"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "Granule neurons from cerebellum"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "CHO"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "Cholinergic neurons"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "CNGA"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "GABAergic neurons from cerebral nuclei"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "COP"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "Committed oligodendrocytes cells"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "CT"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "L6 corticothalamic (CT) projection neurons"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "D12NAC"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "D1-/D2- medium spiny neurons from nucleus Accumbens"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "D1CaB"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "D1- medium spiny neurons from body of the Caudate"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "D1Pu"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "D1- medium spiny neurons from putamen"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "D2CaB"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "D2- medium spiny neurons from body of the Caudate"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "D2Pu"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "D2- medium spiny neurons from putamen"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "EC"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "Endothelial cells"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "ERC"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "Glutamatergic neurons from entorhinal cortex"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "ET"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "Extratelencephalic projecting neurons"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "FOXP2"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "FOXP2+ GABAergic neurons from cerebral nuclei"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "ICGA"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "Dopaminergic neurons from Inferior colliculus and nearby nuclei"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "ITL23"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "Intratelencephalic projecting neurons, cortical layer 2/3"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "ITL34"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "Intratelencephalic projecting neurons, cortical layer 3/4 like"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "ITL4"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "Intratelencephalic projecting neurons, cortical layer 4"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "ITL45"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "Intratelencephalic projecting neurons, cortical layer 4/5 like"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "ITL5"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "Intratelencephalic projecting neurons, cortical layer 5"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "ITL6_1"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "Intratelencephalic projecting neurons, cortical layer 6 - subclass 1"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "ITL6_2"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "Intratelencephalic projecting neurons, cortical layer 6 - subclass 2"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "ITV1C"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "Intratelencephalic projecting neurons from primary visual cortex"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "L6B"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "Intratelencephalic projecting neurons, cortical layer 6B"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "LAMP5"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "LAMP5+ GABAergic neurons"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "LAMP5_LHX6"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "LAMP5+ GABAergic neurons with LHX6+"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "MDGA"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "Dopaminergic neurons from mediodorsal nucleus "
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "MGC"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "Microglia"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "MSN"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "Medium spiny neurons"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "NP"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "Near-projecting neurons"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "OGC"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "Oligodendrocytes"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "OPC"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "Oligodendrocytes precursor cells"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "PIR"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "Glutamatergic neurons from piriform cortex"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "PKJ"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "Purkinje cells"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "PVALB"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "PVALB+ GABAergic neurons"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "PV_ChCs"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "PVALB+ chandelier cells"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "SEPGA"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "Dopaminergic neurons from septal nuclei"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "SIGA"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "Dopaminergic neurons from Inferior colliculus and nearby nuclei"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "SMC"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "Vascular smooth muscle cells"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "SNCG"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "SNCG+ GABAergic neurons"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "SST"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "SST+ GABAergic neurons"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "SST_CHODL"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "SST+ GABAergic neurons with CHODL+"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "SUB"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "Glutamatergic neurons from subicular cortex"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "THMGA"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "GABAergic neurons from thalamus"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "TP"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "Glutamatergic neurons from temporal pole (TP)"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "VIP"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "VIP+ GABAergic neurons"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "CTXMIX"
                           & with_polyp_and_other_blood["dataset"] == "D8",
                           "Abbreviation Meaning"] = "Dopaminergic neurons from midbrain mixed with cell from cerebral nuclei"


with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "AT1/AT2"
                           & with_polyp_and_other_blood["dataset"] == "D9",
                           "Abbreviation Meaning"] = "Alveolar type I/type II"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "Adventitial fibro"
                           & with_polyp_and_other_blood["dataset"] == "D9",
                           "Abbreviation Meaning"] = "Adventitial fibroblast"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "Airway SMC"
                           & with_polyp_and_other_blood["dataset"] == "D9",
                           "Abbreviation Meaning"] = "Airway Smooth Muscle Cells"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "Airway fibro"
                           & with_polyp_and_other_blood["dataset"] == "D9",
                           "Abbreviation Meaning"] = "Airway fibroblasts"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "Alveolar fibro"
                           & with_polyp_and_other_blood["dataset"] == "D9",
                           "Abbreviation Meaning"] = "Alveolar fibroblasts"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "Arterial endo"
                           & with_polyp_and_other_blood["dataset"] == "D9",
                           "Abbreviation Meaning"] = "Arterial endothelial cells"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "CD16+ NK"
                           & with_polyp_and_other_blood["dataset"] == "D9",
                           "Abbreviation Meaning"] = "CD16+ Natural Killer Cells"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "CD56bright NK"
                           & with_polyp_and_other_blood["dataset"] == "D9",
                           "Abbreviation Meaning"] = "CD56bright Natural Killer Cells"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "CX3CR1+ Mac"
                           & with_polyp_and_other_blood["dataset"] == "D9",
                           "Abbreviation Meaning"] = "CX3CR1+ Macrophages"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "DC1"
                           & with_polyp_and_other_blood["dataset"] == "D9",
                           "Abbreviation Meaning"] = "Dendritic cells 1"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "DC2"
                           & with_polyp_and_other_blood["dataset"] == "D9",
                           "Abbreviation Meaning"] = "Dendritic cells 2"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "Def. ery."
                           & with_polyp_and_other_blood["dataset"] == "D9",
                           "Abbreviation Meaning"] = "Definitive Erythrocytes"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "Early fibro"
                           & with_polyp_and_other_blood["dataset"] == "D9",
                           "Abbreviation Meaning"] = "Early fibroblasts"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "GHRL+ NE"
                           & with_polyp_and_other_blood["dataset"] == "D9",
                           "Abbreviation Meaning"] = "GHRL+ Neuroendocrine"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "ILC2"
                           & with_polyp_and_other_blood["dataset"] == "D9",
                           "Abbreviation Meaning"] = "Type 2 Innate Lymphoid Cells"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "ILC3/Th17"
                           & with_polyp_and_other_blood["dataset"] == "D9",
                           "Abbreviation Meaning"] = "IL-17 producing innate lymphoid cells 3"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "ILCP"
                           & with_polyp_and_other_blood["dataset"] == "D9",
                           "Abbreviation Meaning"] = "Innate lymphoid cells (progenitors)"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "Interm fibro"
                           & with_polyp_and_other_blood["dataset"] == "D9",
                           "Abbreviation Meaning"] = "Intermediate fibroblasts"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "Interm lymphatic endo"
                           & with_polyp_and_other_blood["dataset"] == "D9",
                           "Abbreviation Meaning"] = "Intermediate lymphatic endothelial cells"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "Lymphatic endo"
                           & with_polyp_and_other_blood["dataset"] == "D9",
                           "Abbreviation Meaning"] = "Lymphatic endothelial cells"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "Megk"
                           & with_polyp_and_other_blood["dataset"] == "D9",
                           "Abbreviation Meaning"] = "Megakaryocytes"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "Mid/late meso"
                           & with_polyp_and_other_blood["dataset"] == "D9",
                           "Abbreviation Meaning"] = "Mid/late mesothelial cells"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "Myofibro"
                           & with_polyp_and_other_blood["dataset"] == "D9",
                           "Abbreviation Meaning"] = "Myofibroblasts"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "NKT"
                           & with_polyp_and_other_blood["dataset"] == "D9",
                           "Abbreviation Meaning"] = "Natural killer T cells"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "Pri. ery."
                           & with_polyp_and_other_blood["dataset"] == "D9",
                           "Abbreviation Meaning"] = "Primary Erythrocytes"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "Pulmonary NE"
                           & with_polyp_and_other_blood["dataset"] == "D9",
                           "Abbreviation Meaning"] = "Pulmonary Neuroendocrine"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "SCG3+ lymphatic endo"
                           & with_polyp_and_other_blood["dataset"] == "D9",
                           "Abbreviation Meaning"] = "SCG3+ lymphatic endothelial cells"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "SPP1+ Mac"
                           & with_polyp_and_other_blood["dataset"] == "D9",
                           "Abbreviation Meaning"] = "SPP1+ Macrophages"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "Secretory 1/2"
                           & with_polyp_and_other_blood["dataset"] == "D9",
                           "Abbreviation Meaning"] = "Secretory cells (subtypes 1 and 2)"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "Treg"
                           & with_polyp_and_other_blood["dataset"] == "D9",
                           "Abbreviation Meaning"] = "Regulatory T cells"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "Vascular SMC"
                           & with_polyp_and_other_blood["dataset"] == "D9",
                           "Abbreviation Meaning"] = "Vascular smooth muscle cells"
with_polyp_and_other_blood[with_polyp_and_other_blood["cell_type"] == "Venous endo"
                           & with_polyp_and_other_blood["dataset"] == "D9",
                           "Abbreviation Meaning"] = "Venous endothelial cells"

write.csv(with_polyp_and_other_blood, "extended_data_table_2.csv")
write.csv(gl_blood_other, "extended_data_table_3.csv")