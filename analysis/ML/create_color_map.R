library(dplyr)
library(magrittr)

load_data <- function(dataset, cell_filter=50) {
  file = paste(dataset, "combined_count_overlaps_metadata.rds", sep="_")
  fp = paste("../../data/processed_data/count_overlap_data/combined_count_overlaps/finalized_annotation", 
             file, sep="/")
  df = readRDS(fp)
  df["num_cells"] = as.numeric(df[["num_cells"]])
  df = df %>% filter(num_cells >= 50)
  return(df)
}

bingren = load_data("Bingren")
bingren["cell_type_origin"] = paste(bingren[["cell_type"]], "BR")
colnames(bingren)[3] = "tissue_name"
greenleaf_brain = load_data("Greenleaf_brain")
greenleaf_brain["cell_type_origin"] = paste(greenleaf_brain[["cell_type"]], "GL_Br")
shendure = load_data("Shendure")
shendure["cell_type_origin"] = paste(shendure[["cell_type"]], "SH")
greenleaf_colon = load_data("Greenleaf_colon")
greenleaf_colon["cell_type_origin"] = paste(greenleaf_colon[["cell_type"]], "GL_Co")
tsankov = load_data("Tsankov")
tsankov["cell_type_origin"] = paste(tsankov[["cell_type"]], "TS")
yang_kidney = load_data("Yang_kidney")
yang_kidney["cell_type_origin"] = paste(yang_kidney[["cell_type"]], "Y_K")

combined = rbind(bingren, greenleaf_brain, shendure, greenleaf_colon, 
                 tsankov, yang_kidney)

n_colors = length(unique(combined[["tissue_name"]]))
hues = seq(0, 360, length.out = n_colors + 1)[1:n_colors]

hsv_to_rgb <- function(h, s = 1, v = 1) {
  h = h / 360
  rgb = hsv(h, s, v)
  return(rgb)
}

colors <- sapply(hues, hsv_to_rgb)
barplot(rep(1, n_colors), col = colors, border = NA, space = 0, axes = FALSE)

generate_shades <- function(h, n_shades = 5) {
  values <- seq(1, 0.5, length.out = n_shades)
  shades <- sapply(values, function(v) {
    grDevices::hsv(h / 360, 1, v)
  })
  return(shades)
}

# Generate 5 shades for each of the 38 unique hues
all_shades <- lapply(hues, generate_shades)

# Flatten the list to a vector for plotting
shades_vector <- do.call(c, all_shades)

# Display the shades
barplot(rep(1, n_colors * 5), col = shades_vector, border = NA, space = 0, axes = FALSE)
