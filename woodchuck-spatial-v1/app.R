# app.R

source('global.R')

library(shiny)
library(shinydashboard)
library(shinyjs)
library(shinybusy)

#ui <- fluidPage(
#  fileInput(inputId = "spatial_data",
#            label = "Spatial Transcriptomics Data"),
#  tableOutput("files")
#)

ui <- dashboardPage(
  dashboardHeader(title = "Woodchuck liver"),
  dashboardSidebar(
    # Stop side bar from being able to minimize
    tags$head(
      tags$style(HTML(".skin-blue .main-header .sidebar-toggle {display: none;}"))
    ),
    sidebarMenu(id = 'tab',
                useShinyjs(), # Allows anything to be disabled
                menuItem("Home Page", tabName = "home", icon = icon("list")),
                menuItem("Gene Expression", tabName = "input", icon = icon("edit")),
                conditionalPanel(condition = "input.tab == 'input'",
                                 fileInput("file", "Upload File", multiple = FALSE, accept = ".rds")),
                actionButton("reset", "Reset", icon = icon("undo"), style = "color: #fff; background-color: #dc3545; width: 87.25%"),
                actionButton("run", "Run", icon = icon("play"), style = "color: #fff; background-color: #28a745; width: 87.25%")
                )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "input",
              tabsetPanel(id = "main_tabs",
                          tabPanel("Instructions",
                                   # Use text from markdown file
                                   includeMarkdown('instructions.md')
                                   )
                          )),
      tabItem(tabName = "home",
              tags$h1(HTML("Woodchuck liver spatial transcriptomics analysis")))
    )
  )
)

server <- function(input, output, session) {
  options(shiny.maxRequestSize = 260 * 1024^2)
  shinyjs::disable("run") # Disable run by default
  observe({
    # Enable run if there is an input file
    if (is.null(input$file) != TRUE){
      shinyjs::enable("run")
      print(input$file)
    } else {
      shinyjs::disable("run")
    }
  })
  
  # When click reset, remove file and disable run button
  observeEvent(input$reset, {
    shinyjs::reset("file")
    shinyjs::disable("run")
    # Remove specific tabs
    removeTab('main_tabs', 'Metadata')
    removeTab('main_tabs', 'Gene Expression')
  })
  
  # When click run, disable run button and enable spinner widget while object is loading
  observeEvent(input$run, {
    shinyjs::disable("run")
    show_modal_spinner(text = "Preparing plots...")
    obj <- load_seurat_obj(input$file$datapath)
    # If there is an error, print out what the error is
    if (is.vector(obj)) {
      showModeal(modalDialog(
        title = "Error with file",
        HTML("<h5>There is an error with the file you uploaded. See below for more details.</h5><br>",
             paste(unlist(obj), collapse = "<br><br>"))
      ))
      # Enable run widget again
      shinyjs::enable("run")
    } else {
      output$mdPlot <- renderPlot({
        create_metadata_plot(obj, input$metadata_col)
      })
      output$genePlot <- renderPlot({
        create_feature_plot(obj, input$gene)
      })
      
      # Enable downloading of spatial metadata plots
      output$download_md <- downloadHandler(
        filename = function() {
          paste0(input$metadata_col, '_spatialPlot.png')
        },
        content = function(file) {
          plot <- create_metadata_umap(obj, input$metadata_col)
          ggsave(filename = file, width = 10, height = 5, type = "cario")
        }
      )
      
      # Enable downloading of feature plots
      output$download_md <- downloadHandler(
        filename = function() {
          paste0(input$gene, '_spatialPlot.png')
        },
        content = function(file) {
          plot <- create_metadata_umap(obj, input$gene)
          ggsave(filename = file, width = 10, height = 5, type = "cario")
        }
      )
      
      insertTab(
        inputId = "main_tabs",
        tabPanel(
          "Metadata",
          fluidRow(
            column(
              width = 8, # Set plot width
              plotOutput(outputId = "mdPlot"), # Needs to match ID from output
              downloadButton("download_md", "Download metadata plot")
            ),
            column(
              width = 4, # Width of dropdown menu
              selectizeInput("metadata_col",
                              "Metadata Column",
                              colnames(obj@meta.data)
                              )
            )
          )
          )
        )
      insertTab(
        inputId = "main_tabs",
        tabPanel(
          "Gene Expression",
          fluidRow(
            column(
              width = 8, # Set plot width
              plotOutput(outputId = "genePlot"),
              downloadButton("download_genePlot", "Download feature plot")
            ),
            column(
              width = 4, # Width of dropdown menu
              selectizeInput("gene",
                              "Genes",
                              rownames(obj)
              )
            )
          )
        )
      )
      remove_modal_spiner()
      shinyjs::enable("run")
    }
  })
}

shinyApp(ui, server)
