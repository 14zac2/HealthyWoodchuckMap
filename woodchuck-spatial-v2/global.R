library(shiny)
library(Seurat)
library(ggplot2)

create_metadata_plot <- function(sobj, col) {
  if (col %in% c("seurat_clusters","orig.ident","SCT_snn_res.0.8","SCT_snn_res.0.2")) {
    plot1 <- SpatialDimPlot(sobj, alpha = 0.7, group.by = col, images = 'slice1') + theme(legend.position="top")
    plot2 <- SpatialDimPlot(sobj, alpha = 0.7, group.by = col, images = 'slice1.1') + theme(legend.position="top")
  } else if (col %in% colnames(sobj@meta.data)) {
    plot1 <- SpatialFeaturePlot(sobj, features = col, images = 'slice1')
    plot2 <- SpatialFeaturePlot(sobj, features = col, images = 'slice1.1')
  } else {
    plot1 <- ggplot() +
      theme_void() +
      geom_text(aes(x = 0.5, y = 0.5, label = "Metadata doesn't exist."), size = 5, color = "gray73",
                fontface = "bold") +
      theme(plot.margin = unit(c(0,0,0,0), "cm"))
    plot2 <- plot1
  }
  return(plot1 + plot2)
}

create_feature_plot <- function(sobj, gene) {
  if (gene %in% rownames(sobj)) {
    plot3 <- SpatialFeaturePlot(sobj, features = gene, images = 'slice1')
    plot4 <- SpatialFeaturePlot(sobj, features = gene, images = 'slice1.1')
  } else {
    plot3 <- ggplot() +
      theme_void() +
      geom_text(aes(x = 0.5, y = 0.5, label = "Gene doesn't exist."), size = 5, color = "gray73",
                fontface = "bold") +
      theme(plot.margin = unit(c(0,0,0,0), "cm"))
    plot4 <- plot3
  }
  return(plot3 + plot4)
}  

create_violin_plot <- function(sobj, gene) {
  if (gene %in% rownames(sobj)) {
    plot <- VlnPlot(sobj, features = gene, group.by = "seurat_clusters")
  } else {
    plot <- ggplot() +
      theme_void() +
      geom_text(aes(x = 0.5, y = 0.5, label = "Gene doesn't exist."), size = 5, color = "gray73",
                fontface = "bold") +
      theme(plot.margin = unit(c(0,0,0,0), "cm"))
  }
  return(plot)
}


