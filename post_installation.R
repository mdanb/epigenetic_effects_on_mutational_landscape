# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install ArchR from GitHub
devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())

