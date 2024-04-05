

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
write.csv(with_polyp_and_other_blood, "extended_data_table_2.csv")
write.csv(gl_blood_other, "extended_data_table_3.csv")