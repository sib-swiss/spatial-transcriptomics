
In this section, you will find the R code that we will use during the course. We will discuss the code and output during practical sessions.

## Source of data

We will work with spatial transcriptomics data of human breast tumor, that was profiled using 3 10x technologies: Chromium scRNAseq on FFPE-fixed single cells, Visium on FFPE slide, and Xenium  [Janesick et al, bioRxiv](https://www.biorxiv.org/content/10.1101/2022.10.06.510405v1). 

Before starting the exercises, set the working directory to where you have downloaded and unzipped the [data folder]() with the files for the exercises and load the necessary packages. Make sure you download also the R script with the [custom functions]() that will be sourced.
On the materials page, you may download an Rmd file that contains the code presented here.

!!! warning
    When using setwd(), change the path within quotes to where you have saved the data for the exercises


```r
# Change the path here:
setwd("/path/to/whereDataFolderIsSaved/")

library(here)
library(Seurat)
library(BayesSpace)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tibble)
library(scater)
library(grid)
library(hrbrthemes)
library(gridExtra)
library(RColorBrewer)

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
source("utils_new.R")

```

## Import pre-processed Chromium data

We will generate a Seurat object containing the Chromium scRNAseq data. The annotation of the cell types was performed by 10X using marker genes and a cell atlas reference. For the sake of time, we will directly load the saved Chromium data.

```r
# chrom_path <- "path/to/where/filtered_feature_bc_matrix.h5/isSaved"

if(FALSE) {
Counts <- Seurat::Read10X_h5(paste0(chrom_path, "filtered_feature_bc_matrix.h5")) # # 18082 30365 barcodes

MT <- read.csv(paste0(chrom_path, "FFPE_cell_annotations.csv")) # 27472 barcodes
MT <- column_to_rownames(MT, "Barcode")

overlap_cell_barcode <- intersect(rownames(MT), colnames(Counts)) # 27460
length(overlap_cell_barcode)

Counts_ <- Counts[, colnames(Counts) %in% overlap_cell_barcode] # 18082 27460

MT_ <- MT %>%
  filter(rownames(MT) %in% overlap_cell_barcode) # 27460 1

chrom <- CreateSeuratObject(counts = Counts_,
                          meta.data = MT_)

saveRDS(chrom, "./data/chrom_raw.rds")
} else {
chrom <- readRDS("data/chrom_raw.rds")
}

``` 
We generate some QC plots of the scRNAseq data.

```r
chrom <- readRDS("data/chrom_raw.rds")
chrom <- preprocess_SEU(chrom)

chrom_libsize_drop <- chrom$libsize_drop
chrom_mito_drop <- chrom$mito_drop
chrom_lowgenecount_drop <- unlist(chrom[["RNA"]][["lowgenecount_drop"]])

if(any(chrom_libsize_drop)){
  p1 <- plot_Hist_Low_Lib_Sizes(chrom)
} else{
  print("No Chromium cells with low library size detected by scuttle::isOutlier()")
}

if(any(chrom_mito_drop)){
  p2 <- plot_Hist_High_Mito_Props(chrom) # TODO: check plot label
}

if(any(chrom_lowgenecount_drop)){
  p3 <- plot_Hist_Low_Abun_Genes(chrom)
}

```
We retain the cells that pass QC, and save the Chromium reference for later use with cell-type deconvolution per spot.

```r
# Subset
chrom <- chrom[!chrom_lowgenecount_drop, !chrom_libsize_drop & !chrom_mito_drop]

dim(chrom) # 14026 25817

# saveRDS(chrom, "data/chrom_qcd.rds")
```

# Visium analysis

We import a Visium Seurat object, and visualize the per spot total count in 2D coordinates.

```r
vis_path <- "path/to/Visium_FFPEfiles"

?Load10X_spatial

vis <- Load10X_Spatial(vis_path)
Idents(vis) <- ""
SpatialFeaturePlot(vis, features = "nCount_Spatial") + theme(legend.position = "right")

# saveRDS(vis, "data/vis_raw.rds")
```





## Exercise  (the last one :sun_with_face: :scientist_tone3:) - 











