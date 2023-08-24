
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools", version="2.3.2")
}

library(devtools)

# Install ArchR from GitHub
devtools::install_github("GreenleafLab/ArchR", ref="v1.0.1", 
                         repos = BiocManager::repositories())

