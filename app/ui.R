library(shiny)
library(shinyFiles)
library(shinyWidgets)
library(shinythemes)
library(shinysky)
library(dqshiny)
library(dplyr)

fpath <- '../'

# Define UI
ui <- fluidPage(theme = shinytheme("paper"),
                navbarPage(
                  "PLO(SC)²: PLOts and SCripts for SCrna-seq",
                  tabPanel(
                    "Select File",
                    sidebarPanel(
                      selectInput('selectfile','Select File',choice = list.files(fpath, pattern = ".Rds")),
                      actionButton("load_selected_file", "Load File")
                    ), #sidebarPanel
                    mainPanel("Main Panel",
                              plotOutput("dimplot", width="auto", height = "800px"),
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
                                     h4("Data Selection"),
                                     plotOutput("vplot", width="100%", height = "800px"),
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
                                     h4("Data Selection"),
                                     plotOutput("enhanced_plot", width="100%", height = "800px"),
                                     style = "font-size:70%"
                           ) # mainPanel
                           
                  ), # Navbar 1, tabPanel
                  
                  
                  
                ) # navbarPage
) # fluidPage



