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

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", version="1.30.22")

if (!requireNamespace("ArchR", quietly = TRUE)) {
	devtools::install_github("GreenleafLab/ArchR", ref="v1.0.1", 
                         	 repos = BiocManager::repositories(), dependencies=T)
}

if (!requireNamespace("seurat5", quietly = TRUE)) {
	remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE, dependencies=T)
}

if (!requireNamespace("exomeCopy", quietly = TRUE)) {
        BiocManager::install("exomeCopy", version="3.17")
}

