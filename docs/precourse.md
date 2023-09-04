
## R and RStudio

### Previous knowledge / Competencies

This tutorial is designed for individuals interested in spatial transcriptomics data analysis. Participants must possess a basic understanding of R and RNA-seq gene expression method and data format. A basic understanding of single-cell RNA sequencing data analysis with R is a plus. Participants should bring a laptop with Wi-Fi, the latest versions of R, Rstudio, and relevant R packages installed (see below).

### Technical

To be able to perform the exercises, participants should install the latest the version of [R](https://cran.r-project.org/)
and the free version of [RStudio](https://www.rstudio.com/products/rstudio/download/).

Install the necessary packages using:

```r
# Packages to install

install.packages("here")
# install.packages("Seurat") # do not install version v. 4.9.9.9050
install.packages("https://cran.r-project.org/src/contrib/Seurat_4.3.0.1.tar.gz", type = "source", repos = NULL)
install.packages("https://cran.r-project.org/src/contrib/SeuratObject_4.1.3.tar.gz", type="source", repos = NULL)
install.packages("BiocManager")
BiocManager::install("BayesSpace")
install.packages("ggplot2")
install.packages("patchwork")
install.packages("dplyr")
install.packages("tibble")
BiocManager::install("scater")
install.packages("hrbrthemes")
install.packages("gridExtra")
install.packages("RColorBrewer")
install.packages("devtools")
install.packages("hdf5r")

install.packages("https://cran.r-project.org/src/contrib/Seurat_4.3.0.1.tar.gz", type = "source", repos = NULL)
install.packages("https://cran.r-project.org/src/contrib/SeuratObject_4.1.3.tar.gz", type="source", repos = NULL)

# If using RCTD:
# devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)

```

