library(shiny)
library(shinydashboard)
library(shinyjs)
library(tools)
library(Seurat)
library(ggplot2)

load_seurat_obj <- function(path) {
  errors <- c()
  # Check file extension
  if (!tolower(tools::file_ext(path)) == 'rds') {
    errors <- c(errors, "Invalid rds file.")
    return(errors)
  }
  # Try to read in file
  tryCatch(
    {
      obj <- readRDS(path)
    },
    error = function(e) {
      errors <- c(errors, "Invalid rds file.")
      return(errors)
    }
  )
  # Validate obj is a Seurat object
  if (!inherits(obj, "Seurat")) {
    errors <- c(errors, "File is not a Seurat object.")
    return(errors)
  }
  return(obj)
}

create_metadata_plot <- function(obj, col) {
  if (col %in% c("seurat_clusters")) {
    plot1 <- SpatialDimPlot(sobj, alpha = 0.7, group.by = "seurat_clusters", images = 'slice1') & NoLegend()
    plot2 <- SpatialDimPlot(sobj, alpha = 0.7, group.by = "seurat_clusters", images = 'slice1.1') +
      labs(fill = "Cluster")
  } else if (col %in% colnames(obj@meta.data)) {
    plot1 <- SpatialFeaturePlot(sobj, features = col, images = 'slice1') & NoLegend()
    plot2 <- SpatialFeaturePlot(sobj, features = col, images = 'slice1.1') & NoLegend()
  } else {
    plot1 <- ggplot() +
      theme_void() +
      geom_text(aes(x = 0.5, y = 0.5, label = "Metadata doesn't exist."), size = 20, color = "gray73",
                fontface = "bold") +
      theme(plot.margin = unit(c(0,0,0,0), "cm"))
    plot2 <- plot1
  }
  return(plot1 + plot2)
}

create_feature_plot <- function(obj, gene) {
  if (gene %in% rownames(obj)) {
    plot3 <- SpatialFeaturePlot(sobj, features = gene, images = 'slice1') & NoLegend()
    plot4 <- SpatialFeaturePlot(sobj, features = gene, images = 'slice1.1') & NoLegend()
  } else {
    plot3 <- ggplot() +
      theme_void() +
      geom_text(aes(x = 0.5, y = 0.5, label = "Gene doesn't exist."), size = 20, color = "gray73",
                fontface = "bold") +
      theme(plot.margin = unit(c(0,0,0,0), "cm"))
    plot4 <- plot3
  }
  return(plot3 + plot4)
}  
  
  
  
