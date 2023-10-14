# app.R

source('global.R')

library(shiny)
library(Seurat)
library(ggplot2)
library(markdown)
library(rsconnect)

# Increase maximum size of object to load
options(shiny.maxRequestSize = 260 * 1024^2)

# Load spatial data
sobj <- readRDS("merged_spatial_sobj.RDS")

ui <- fluidPage(
  titlePanel("Interactive visualization of woodchuck spatial transcriptomics data"),
  
  br(),
  
  fluidRow(
    column(width = 12,
           uiOutput("text_info_md"),
           br()
           )
  ),
  
  fluidRow(
    column(width = 2,
           # Dropdown menu for metadata
           selectInput("metadata_col", "Select metadata to plot:", choices = colnames(sobj@meta.data)),
           # Create button to download metadata plot
           downloadButton("download_md", "Download metadata plot")
           ),
    column(width = 8,
           # Create metadata plot
           plotOutput("mdPlot")
           )
  ),
  
  # Add a break between the rows
  br(),
  
  # Add explanatory text between the two plots
  fluidRow(
    column(width = 12,
           uiOutput("text_info_genes"),
           br()
           )
  ),
  
  fluidRow(
    column(width = 2,
           # Dropdown menu of genes
           selectInput("gene_spatial", "Select gene to plot:", choices = rownames(sobj)),
           # Create button to download gene plot
           downloadButton("download_genePlot", "Download feature plot")
           ),
    column(width = 8,
           # Create feature plot
           plotOutput("genePlot")
           )
  ),
  
  br(),
  
  fluidRow(
    column(width = 12,
           uiOutput("text_info_vln"),
           br()
           )
  ),
  
  fluidRow(
    column(width = 2,
           # Dropdown menu of genes
           selectInput("gene_vln", "Select gene to plot:", choices = rownames(sobj)),
           # Create button to download gene plot
           downloadButton("download_vlnPlot", "Download violin plot")
           ),
    column(width = 8,
           # Create violin plot
           plotOutput("vlnPlot")
           )
  )
)

server <- function(input, output, session) {
  
  # Set resolution for spatial object
  Idents(sobj) <- "SCT_snn_res.0.2"
  
  # Create metadata plot
  output$mdPlot <- renderPlot({
    create_metadata_plot(sobj, input$metadata_col)
  })
  
  # Create spatial feature plots
  output$genePlot <- renderPlot({
    create_feature_plot(sobj, input$gene_spatial)
  })
  
  # Create violin plots
  output$vlnPlot <- renderPlot({
    create_violin_plot(sobj, input$gene_vln)
  })
  
  # Explanatory text for metadata
  output$text_info_md <- renderUI({
    HTML("<p>Below you can explore different metadata for the spatial transcriptomics data. If you scroll down to Periportal_UCell or Pericentral UCell, this gives you an idea of where the periportal and pericentral regions are. These align to clusters 0 and 1 respectively if you navigate to 'seurat_clusters'.</p>")
  })
  
  # Explanatory text for genes
  output$text_info_genes <- renderUI({
    HTML("<p>The following plots show the expression of specific genes selected in the dropdown menu. You can select a gene by scrolling through the dropdown menu or by typing in part of a gene name you are interested in until it appears. Note that this genome contains similar gene symbols separated by semi-colons if they are a one-to-many ortholog or -1, -2, -3, etc. if they are a many-to-one ortholog.</p>")
  })
  
  # Explanatory text for violin plot
  output$text_info_vln <- renderUI({
    HTML("<p>Here you can look at how the expression of each gene is related to the clustering of the spatial transcriptomics spots. Cluster 0 aligns with the periportal region, and cluster 1 aligns with the pericentral region.</p>")
  })
  
  # Enable downloading of spatial metadata plots
  output$download_md <- downloadHandler(
    filename = function() {
      paste0(input$metadata_col, '_spatialPlot.png')
    },
    content = function(file) {
      plot <- create_metadata_plot(sobj, input$metadata_col)
      ggsave(filename = file, width = 10, height = 5, type = "cario")
    }
  )
  
  # Enable downloading of feature plots
  output$download_genePlot <- downloadHandler(
    filename = function() {
      paste0(input$gene_spatial, '_spatialPlot.png')
    },
    content = function(file) {
      plot <- create_feature_plot(sobj, input$gene_spatial)
      ggsave(filename = file, width = 10, height = 5, type = "cario")
    }
  )
  
  # Enable downloading of violin plots
  output$download_vlnPlot <- downloadHandler(
    filename = function() {
      paste0(input$gene_vln, '_violinPlot.png')
    },
    content = function(file) {
      plot <- create_violin_plot(sobj, input$gene_vln)
      ggsave(filename = file, width = 10, height = 5, type = "cario")
    }
  )
  
  # Stop when the browswer tab is closed
  session$onSessionEnded(stopApp)
  
}

shinyApp(ui, server)

