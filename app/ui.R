library(shiny)
library(shinyFiles)
library(shinyWidgets)
library(shinythemes)
library(shinysky)
library(shinyjs)
library(dqshiny)
library(dplyr)
library(bslib)

fpath <- '../'


# Define UI
ui <- fluidPage(theme = shinytheme("paper"),
                shinyjs::useShinyjs(),
                navbarPage(
                  "PLO(SC)²: PLOts and SCripts for SCrna-seq",
                  
                  tabPanel(
                    "Select File",
                    sidebarPanel(
                      selectInput('selectfile','Select File',choice = c("-", list.files(fpath, pattern = ".Rds"), selected="-")),
                      actionButton("load_selected_file", "Load File"),
                      actionButton("reload_files", "Reload List"),
                      textOutput("curpath"),
                      selectInput("u_groupby1", "Group-by Attribute", choices = NULL)
                    ), #sidebarPanel
                    mainPanel("Main Panel",
                              plotOutput("dimplot", width="auto", height = "800px"),
                              downloadButton("download_dimplot", "Download DimPlot"),
                              dataTableOutput("ftxtout")
                    ) # mainPanel
                  ), #tabPanel
                  tabPanel("Enhanced Violin Plot",
                           sidebarPanel(
                             p("Main Options"),
                             selectInput("v_plotmode1", "Plot Mode", choices = c("Violin-Boxplot", "Violin-Boxplot-Comparison", "Side-by-side Violin-Boxplot")),
                             selectInput("v_groupby1", "Group-by Attribute", choices = NULL),
                             selectInput("v_splitby1", "Split-by Attribute", choices = NULL),
                             autocomplete_input("v_feature1", "Feature", options = c(), value = ""),
                             hr(),
                             p("Details"),
                             numericInput("v_ymin1", "ymin", 0),
                             numericInput("v_ymax1", "yMax", 6)
                           ), #sidebarPanel
                           
                           mainPanel(tags$br(),tags$br(),
                                     h4("Enhanced Violin Plot"),
                                     plotOutput("vplot", width="100%", height = "800px"),
                                     downloadButton("download_vplot", "Download Enhanced Violin Plot"),
                                     style = "font-size:70%"
                           ) # mainPanel,
                           
                  ), # Navbar 1, tabPanel
                  tabPanel("Enhanced DotPlot/Heatmap",
                           sidebarPanel(
                             
                             selectInput("v_plotmode2", "Plot Mode", choices = c("Enhanced DotPlot", "Enhanced Heatmap")),
                             selectInput("v_scalingmode2", "Scaling Mode", choices=c("ALL", "GROUP", "FEATURE", "GLOBAL"), selected = "ALL"),
                             selectInput("v_groupby2", "Group-by Attribute", choices = NULL),
                             selectInput("v_splitby2", "Split-by Attribute", choices = NULL),
                             selectizeInput(
                               inputId = 'v_features2',
                               label = 'Search',
                               choices = NULL,
                               selected = NULL,
                               multiple = TRUE, # allow for multiple inputs
                               options = list(create = FALSE) # if TRUE, allows newly created inputs
                             ),
                             textOutput("t_dotplot_warning")
                             
                             
                           ), #sidebarPanel
                           
                           mainPanel(tags$br(),tags$br(),
                                     h4("Enhanced DotPlot/Heatmap"),
                                     plotOutput("enhanced_plot", width="100%", height = "800px"),
                                     downloadButton("download_enhanced_plot", "Download Enhanced DotPlot/Heatmap"),
                                     style = "font-size:70%"
                           ) # mainPanel
                           
                  ), # Navbar 1, tabPanel
                  tabPanel("2D UMAP",
                           sidebarPanel(
                             selectInput("du_groupby", "Group-by Attribute", choices = NULL),
                             selectInput("du_scalingmode", "Scaling Mode", choices=c("All cells", "Downsample"), selected = "All cells"),
                             selectInput("du_xsplit", "X-Split-by Attribute", choices = NULL),
                             selectInput("du_ysplit", "Y-Split-by Attribute", choices = NULL)
                             
                           ), #sidebarPanel
                           
                           mainPanel(tags$br(),tags$br(),
                                     h4("2D UMAP"),
                                     plotOutput("du_umap", width="100%", height = "800px"),
                                     downloadButton("download_2d_umap", "Download 2D UMAP"),
                                     style = "font-size:70%"
                           ) # mainPanel
                           
                  ), # Navbar 1, tabPanel
                  tabPanel("Process from Scratch",
                           sidebarPanel(
                             selectInput('scratch_folder','Select Folder', choice = c("-", grep(x=list.dirs(fpath, recursive = FALSE), pattern = "\\.\\/\\.", value = TRUE, invert = TRUE)), selected = "-"),
                             actionButton("scratch_reload_folders", "Reload Folders"),
                             
                             
                             selectInput("scratch_annotation_samples", "Sample Annotation", choices = NULL),
                             
                             selectInput("scratch_organism", "Organism", choices=c("Human", "Mouse"), selected = "Human"),
                             
                             #selectInput("scratch_annotation_hto", "HTO Annotation", choices = NULL),
                             
                             selectInput("scratch_integration_method", "Integration Method", choices=c("rpca", "cca"), selected = "rpca"),
                             numericInput("scratch_integration_features", "Number of Integration Features", 3000),
                             
                             
                             numericInput("scratch_minfeatures", "nFeature_RNA min", 100, min=100),
                             numericInput("scratch_maxfeatures", "nFeature_RNA max", 6000, min=0),
                             numericInput("scratch_minncount", "nCount_RNA min", 500, min=0),
                             numericInput("scratch_maxpercentmt", "PercentMT max", 7, min=0),
                             
                             numericInput("scratch_umap_dims", "PCs for UMAP", 7, min=2),
                             numericInput("scratch_cluster_resolution", "Cluster Resolution", 0.5, min=0),
                             
                             actionButton("scratch_go_button", "Process!"),tags$br(),
                             downloadButton("scratch_downloadobj", "Download Seurat Object"),
                             
                           ), #sidebarPanel
                           
                           mainPanel(tags$br(),tags$br(),
                                     h4("Processed Data UMAP"),
                                     plotOutput("scratch_umap", width="100%", height = "800px"),
                                     style = "font-size:70%"
                           ) # mainPanel
                           
                  ),
                  tabPanel("Differential Gene Expression",
                           sidebarPanel(
                             selectInput("de_groupby", "Group-by Attribute", choices = NULL),
                             selectInput("de_method", "DE Test Method", choices=c("t", "wilcox", "bimod", "negbinom", "poisson", "LR"), selected = "t"),
                             actionButton("de_process", "Process!"),tags$br(),
                             
                             
                           ), #sidebarPanel
                           
                           mainPanel(tags$br(),tags$br(),
                                     h4("Differential Gene Expression"),
                                     dataTableOutput("de_out"),
                                     shinyjs::hidden(downloadButton("de_download", "Download DEG Results")),
                                     
                                     style = "font-size:70%"
                           ) # mainPanel
                           
                  ),
                  tabPanel("Gene Set Enrichment",
                           sidebarPanel(
                             selectInput('gse_folder','Select Folder', choice = grep(x=list.dirs(fpath, recursive = FALSE), pattern = "\\.\\/\\.", value = TRUE, invert = TRUE)),
                             selectInput('gse_de_file','DE Results', choice = c("-")),
                             checkboxInput("gse_plot_volcanos", "Add Volcano Plots", value = FALSE, width = NULL),
                             
                             selectInput("gse_organism", "Organism", choices=c("human", "mouse"), selected = "Human"),
                             numericInput("gse_min_foldchange", "Min FoldChange", 0.25, min=0, max=1),
                             
                             actionButton("gse_process", "Process!"),tags$br(),
                             
                           ), #sidebarPanel
                           
                           mainPanel(tags$br(),tags$br(),
                                     h4("Gene Set Enrichment"),
                                     dataTableOutput("gse_inputfile"),
                                     style = "font-size:70%"
                           ) # mainPanel
                           
                  ),
                  tabPanel("Settings",
                           sidebarPanel(
                             selectInput("plot_output", "Plot Download Format", choices=c("png", "pdf"), selected = "png"),
                             
                           ), #sidebarPanel
                           
                           mainPanel() # mainPanel
                           
                  ),
                ) # navbarPage
) # fluidPage



