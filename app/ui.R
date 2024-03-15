library(shiny)
library(shinyFiles)
library(shinyWidgets)
library(shinythemes)
library(shinysky)
library(dqshiny)
library(dplyr)

fpath <- '../'

# Define UI
ui <- fluidPage(theme = shinytheme("spacelab"),
                navbarPage(
                  "PLO(SC)²: PLOts and SCripts for SCrna-seq",
                  tabPanel(
                    "Select File",
                    sidebarPanel(
                      selectInput('selectfile','Select File',choice = list.files(fpath, pattern = ".rds"))
                    ), #sidebarPanel
                    mainPanel("Main Panel",
                              plotOutput("dimplot"),
                              dataTableOutput("ftxtout")
                              ) # mainPanel
                  ), #tabPanel
                  tabPanel("Violin Plot",
                           sidebarPanel(
      
                             selectInput("v_plotmode1", "Plot Mode", choices = c("Violin-Boxplot", "Violin-Boxplot-Comparison", "Side-by-side Violin-Boxplot")),
                             selectInput("v_groupby1", "Group-by Attribute", choices = NULL),
                             selectInput("v_splitby1", "Split-by Attribute", choices = NULL),
                             autocomplete_input("v_feature1", "Feature", options = c(), value = "")
                             

                           ), #sidebarPanel
                           
                           mainPanel(tags$br(),tags$br(),
                                     h4("Data Selection"),
                                     plotOutput("vplot", height = 400),
                                     style = "font-size:70%"
                           ) # mainPanel
                           
                  ), # Navbar 1, tabPanel
                  tabPanel("DotPlot Plot",
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
                             )
                                                          
                             
                           ), #sidebarPanel
                           
                           mainPanel(tags$br(),tags$br(),
                                     h4("Data Selection"),
                                     plotOutput("enhanced_plot", height = 400),
                                     style = "font-size:70%"
                           ) # mainPanel
                           
                  ), # Navbar 1, tabPanel
                  
                  
                  
                ) # navbarPage
) # fluidPage



