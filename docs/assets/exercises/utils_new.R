preprocess_SEU <- function(seu){
  seu$percent.mt <- PercentageFeatureSet(seu, pattern = "^MT-")
  gene_means <- as.numeric(unlist(rowMeans(GetAssayData(seu, "counts"), na.rm = TRUE)))
  
  ## Per cell
  libsize_drop <- isOutlier(
    metric = as.numeric(seu@meta.data[, grepl("nCount", colnames(seu@meta.data))]), 
    type = "lower",
    log = TRUE) 
  
  mito_drop <- isOutlier(
    metric = seu$percent.mt,
    type = "higher") 
  
  seu$libsize_drop <- libsize_drop
  seu$mito_drop <- mito_drop
  
  ## Per gene
  lowgenecount_drop <- log(gene_means) < -5 | gene_means <= 0
  
  seu[[names(seu@assays)]][["lowgenecount_drop"]] <- lowgenecount_drop
  seu[[names(seu@assays)]][["means"]] <- gene_means
  
  return(seu)
}


plotQC_bar <- function(plot_df, xvar = "logtotal", yvar = "libsize_drop"){
  p <- plot_df %>%
    ggplot(aes(x = get(xvar), fill = get(yvar))) +
    geom_histogram(color = "#e9ecef", alpha = 0.6, position = 'identity', bins = 30) +
    scale_fill_manual(values = c("#404080", "#69b3a2")) +
    theme_ipsum() +
    labs(fill = "") +
    theme(plot.title = element_text(hjust = 0.5, face = "plain", size = 12))
  
  p
}


plot_Hist_High_Mito_Props <- function(seu){
  plot_df <- data.frame(mito_percent = log(seu$percent.mt),
                        mito_drop = factor(seu$mito_drop))
  
  plot_df$mito_drop <- relevel(plot_df$mito_drop, "TRUE")
  
  p <- plotQC_bar(plot_df, xvar = "mito_percent", yvar = "mito_drop")
    
  p + xlab("Cell level mitochondria percent") + ylab("Frequency") + 
    ggtitle('Cells With High Mitochondria Proportion')
  
  p
}


plot_Hist_Low_Abun_Genes <- function(seu){
  plot_df <- data.frame(mean_genecount = log(unlist(seu[[names(seu@assays)]][["means"]])),
                        lowgenecount_drop = factor(unlist(seu[[names(seu@assays)]][["lowgenecount_drop"]])))
  
  plot_df$lowgenecount_drop <- relevel(plot_df$lowgenecount_drop, "TRUE")
  
  p <- plotQC_bar(plot_df, xvar = "mean_genecount", yvar = "lowgenecount_drop")
  
  p <- p + xlab("Log mean count across all cells") + ylab("Frequency") + 
    ggtitle('Low Abundance Genes')
  
  p
}


plot_Hist_Low_Lib_Sizes <- function(seu){
  plot_df <- data.frame(logtotal=log(as.numeric(seu@meta.data[, grepl("nCount", colnames(seu@meta.data))])),
                        libsize_drop=factor(seu$libsize_drop))
  
  plot_df$libsize_drop <- relevel(plot_df$libsize_drop, "TRUE")
  
  p <- plotQC_bar(plot_df, xvar = "logtotal", yvar = "libsize_drop")
  
  p <- p + xlab("Cell level log total count") + ylab("Frequency") + 
    ggtitle('Low Library Size Cells')
  
  p
}

# Spatial plot for continuous features
SpatialFeaturePlot_cont <- function(xe, color_by = "ABCC11"){
  CD <- xe@meta.data
  feat_mat <- GetAssayData(object = xe, slot = "counts")
  CD <- cbind(CD, t(feat_mat))
  
  ggplot(CD,
         aes(x = array_row, y = array_col,
             color = get(color_by))) +
    geom_point(size = 0.3) +
    theme_void() + 
    theme(legend.position="right") +
    scale_colour_gradientn(name = color_by, colors = myPalette(100))
}

# Spatial plot for binary flags 
plotQC_seu <- function(xe, flag = "libsize_drop"){
  if(!("array_row" %in% colnames(xe@meta.data))){
    xe[[c("array_row", "array_col")]] <- GetTissueCoordinates(xe)
    xe[["cell_id"]] <- colnames(xe)
  }
  df <- data.frame(xe[[c("cell_id", "array_col", "array_row", "nCount_Spatial", "nFeature_Spatial", flag)]])
  ggplot(df, aes(x = array_row, y = array_col, color = get(flag))) + 
    geom_point(size = 0.3) + 
    coord_fixed() + 
    scale_color_manual(name = flag, values = c("gray85", "red")) + 
    ggtitle("QC dots") + 
    theme_bw() +
    theme(panel.grid = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
}

# Spatial plot for categorical value with self-defined palette
SpatialFeaturePlot_cate <- function(xe, color_by = "seurat_clusters", palette = NULL){
  CD <- xe@meta.data
  
  ggplot(CD,
         aes(x = array_row, y = array_col,
             color = get(color_by))) +
    geom_point(size = 0.3) +
    theme_void() + 
    theme(legend.position="right") +
    scale_color_manual(name = color_by, values = palette) +
    guides(colour = guide_legend(override.aes = list(size = 3)))
}


extractCol <- function(pumap){
  group <- as.numeric(ggplot_build(pumap)$data[[1]]$group)
  col <- ggplot_build(pumap)$data[[1]]$colour
  col <- forcats::fct_reorder(col, group)
  cols_Red <- levels(col)
  
  return(cols_Red)
}
