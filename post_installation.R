# Set default CRAN mirror
local({
  r <- getOption("repos")
  r["CRAN"] <- "https://cloud.r-project.org/"
  options(repos = r)
})

if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools") 
}

library(devtools)


install_version("lattice", version="0.22.6")
# Set the URL for the specific version you need
url <- "https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.6-4.tar.gz"

# Set the destination file path
destfile <- tempfile()

# Download the package tarball
download.file(url, destfile, mode = "wb")

# Install the package from the downloaded tarball
install.packages(destfile, repos = NULL, type = "source")

# Remove the temporary file if you want to clean up
unlink(destfile)

# Set the URL for the specific version you need
url <- "https://cran.r-project.org/src/contrib/Archive/MASS/MASS_7.3-50.tar.gz"

# Set the destination file path
destfile <- tempfile()

# Download the package tarball
download.file(url, destfile, mode = "wb")

# Install the package from the downloaded tarball
install.packages(destfile, repos = NULL, type = "source")

# Remove the temporary file if you want to clean up
unlink(destfile)


#if (!requireNamespace("devtools", quietly = TRUE)) {
#  install.packages("devtools")
#}

#library(devtools)

if (!requireNamespace("tidyverse", quietly = TRUE)) {
    install_version('tidyverse', version="2.0.0", dependencies=T)
}
remove.packages("ggplot2")
install_version('ggplot2', version="3.4.0", dependencies=T)


if (!requireNamespace("data.table", quietly = TRUE)) {
    install_version("data.table", version="1.14.8")
}

if (!requireNamespace("optparse", quietly = TRUE)) {
    install_version("optparse", version="1.7.3")
}

if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    install_version("RColorBrewer", version="1.1.3")
}

if (!requireNamespace("IRkernel", quietly = TRUE)) {
    install_version("IRkernel", version="1.3.2")
}

if (!requireNamespace("viridis", quietly = TRUE)) {
    install_version("viridis", version="0.6.3")
}

if (!requireNamespace("tidytext", quietly = TRUE)) {
    install_version("tidytext", version="0.4.1")
}

if (!requireNamespace("usethis", quietly = TRUE)) {
    install_version("usethis", version="2.2.2")
}

#if (!requireNamespace("devtools", quietly = TRUE)) {
#  install_version("devtools", version="2.4.5")
#}

if (!requireNamespace("hash", quietly = TRUE)) {
  install_version("hash", version="2.2.6.3")
} 

if (!requireNamespace("paletteer", quietly = TRUE)) {
  install_version("paletteer", version="1.6.0")
}

if (!requireNamespace("svglite", quietly = TRUE)) {
  install_version("svglite", version="2.1.2")
}

if (!requireNamespace("this.path", quietly = TRUE)) {
  install_version("this.path", version="2.4.0")
}


if (!requireNamespace("BiocManager", quietly = TRUE))
    install_version("BiocManager", version="1.30.22")
BiocManager::install(version = "3.14")

if (!requireNamespace("ArchR", quietly = TRUE)) {
	devtools::install_github("GreenleafLab/ArchR", ref="v1.0.1", 
                         	 repos = BiocManager::repositories(), dependencies=T)
}

if (!requireNamespace("multtest", quietly = TRUE)) {
        BiocManager::install("multtest", version="3.14")
}

if (!requireNamespace("seurat", quietly = TRUE)) {
	install_version('Seurat', version="4.4.0", dependencies=T)
}

if (!requireNamespace("DoubletFinder", quietly = TRUE)) {
        remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', quiet=T, dependencies=T)
}

if (!requireNamespace("exomeCopy", quietly = TRUE)) {
        BiocManager::install("exomeCopy", version="3.14")
}

if (!requireNamespace("harmony", quietly = TRUE)) {
	devtools::install_github("immunogenomics/harmony@63ebd73")
}

if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
        BiocManager::install("ComplexHeatmap", version="3.14")
}

if (!requireNamespace("edgeR", quietly = TRUE)) {
        BiocManager::install("edgeR", version="3.14")
}

if (!requireNamespace("preprocessCore", quietly = TRUE)) {
        BiocManager::install("preprocessCore", version="3.14")
}

if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)) {
        BiocManager::install("BSgenome.Hsapiens.UCSC.hg19", version="3.14")
}

if (!requireNamespace("plyranges", quietly = TRUE)) {
        BiocManager::install("plyranges", version="3.14")
}

if (!requireNamespace("colorRamps", quietly = TRUE)) {
        BiocManager::install("colorRamps", version="2.3.4")
}

