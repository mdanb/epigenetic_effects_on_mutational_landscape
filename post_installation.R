# Set default CRAN mirror
local({
  r <- getOption("repos")
  r["CRAN"] <- "https://cloud.r-project.org/"
  options(repos = r)
})

if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools", version="2.3.2")
}

library(devtools)
devtools::install_github("GreenleafLab/ArchR", ref="v1.0.1", 
                         repos = BiocManager::repositories(), dependencies=T)
remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE, dependencies=T)
