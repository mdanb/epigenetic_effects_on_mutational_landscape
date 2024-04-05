# Set default CRAN mirror
local({
  r <- getOption("repos")
  r["CRAN"] <- "https://cloud.r-project.org/"
  options(repos = r)
})

if (!requireNamespace("tidyverse", quietly = TRUE)) {
    install.packages('tidyverse', version="2.0.0")
}

if (!requireNamespace("data.table", quietly = TRUE)) {
    install.packages("data.table", version="1.14.8")
}

if (!requireNamespace("optparse", quietly = TRUE)) {
    install.packages("optparse", version="1.7.3")
}

if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    install.packages("RColorBrewer", version="1.1.3")
}

if (!requireNamespace("IRkernel", quietly = TRUE)) {
    install.packages("IRkernel", version="1.3.2")
}

if (!requireNamespace("viridis", quietly = TRUE)) {
    install.packages("viridis", version="0.6.3")
}

if (!requireNamespace("tidytext", quietly = TRUE)) {
    install.packages("tidytext", version="0.4.1")
}

if (!requireNamespace("usethis", quietly = TRUE)) {
    install.packages("usethis", version="2.2.2")
}

if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools", version="2.3.2")
}

if (!requireNamespace("hash", quietly = TRUE)) {
  install.packages("hash", version="2.2.6")
} 

if (!requireNamespace("hash", quietly = TRUE)) {
  install.packages("paletteer", version="1.6.0")
}

if (!requireNamespace("svglite", quietly = TRUE)) {
  install.packages("svglite", version="2.1.2")
}

if (!requireNamespace("this.path", quietly = TRUE)) {
  install.packages("this.path", version="2.4.0")
}


library(devtools)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", version="1.30.22")

if (!requireNamespace("ArchR", quietly = TRUE)) {
	devtools::install_github("GreenleafLab/ArchR", ref="v1.0.1", 
                         	 repos = BiocManager::repositories(), dependencies=T)
}


if (!requireNamespace("seurat5", quietly = TRUE)) {
	install.packages('Seurat', version="4.3.0")
}

if (!requireNamespace("DoubletFinder", quietly = TRUE)) {
        remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', quiet=T, dependencies=T)
}

if (!requireNamespace("exomeCopy", quietly = TRUE)) {
        BiocManager::install("exomeCopy", version="3.18")
}

if (!requireNamespace("harmony", quietly = TRUE)) {
	devtools::install_github("immunogenomics/harmony@63ebd73")
}

if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
        BiocManager::install("ComplexHeatmap", version="3.18")
}

if (!requireNamespace("edgeR", quietly = TRUE)) {
        BiocManager::install("edgeR", version="3.18")
}

if (!requireNamespace("preprocessCore", quietly = TRUE)) {
        BiocManager::install("preprocessCore", version="3.18")
}

