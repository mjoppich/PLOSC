
library(shiny)
library(shinyFiles)
library(shinyWidgets)
library(shinythemes)
library(shinysky)
library(dqshiny)
library(dplyr)
library(Seurat)

#devtools::install_github("AnalytixWare/ShinySky")
#remotes::install_github("daqana/dqshiny")

#setwd("D:/git/PLOSC")

fpath <- '../'

# remotes::install_local(".", force=TRUE); devtools::reload(".");
# devtools::load_all("./")

# Define server function
server <- function(input, output, session) {
  
  makeViolinPlot = function(input, output)
  {
    print(input$v_plotmode1)
    print(input$v_groupby1)
    print(input$v_splitby1)
    print(input$v_feature1)
    
    plotMode = input$v_plotmode1
    groupby = input$v_groupby1
    splitby = input$v_splitby1
    feature = input$v_feature1
    
    if (!(groupby %in% colnames(seuratObj()@meta.data)))
    {
      groupby = NULL
    }
    if (!(splitby %in% colnames(seuratObj()@meta.data)))
    {
      splitby = NULL
    }
    
    print(paste(feature, groupby, splitby))
    
    if (feature %in% rownames(seuratObj()))
    {
      
      if (plotMode == "Violin-Boxplot")
      {
        output$vplot = renderPlot({
          PLOSC::VlnBoxPlot(seuratObj(), gene = feature, group.by = groupby, split.by=splitby)
        }, res = 96)
        
      } else if (plotMode == "Violin-Boxplot-Comparison")
      {
        output$vplot = renderPlot({
          PLOSC::comparativeVioBoxPlot(seuratObj(), feature = feature, group.by = groupby, split.by = splitby)
        }, res = 96)
        
      } else if (plotMode == "Side-by-side Violin-Boxplot")
      {
        output$vplot = renderPlot({
          PLOSC::SplitVlnBoxPlot(seuratObj(), gene = feature, group.by = groupby, split.by=splitby)
        }, res = 96)
        
      }
      

    } else {
      output$vplot = NULL
    }
  }
  
  makeEnhancedPlots = function(input, output)
  {
    print("makeEnhancedPlots")
    print(input$v_plotmode2)
    print(input$v_groupby2)
    print(input$v_splitby2)
    print(input$v_features2)
    print(input$v_scalingmode2)
    
    plotMode = input$v_plotmode2
    scalingMode = input$v_scalingmode2
    
    groupby = input$v_groupby2
    splitby = input$v_splitby2
    features = input$v_features2
    
    if (!(groupby %in% colnames(seuratObj()@meta.data)))
    {
      groupby = NULL
    }
    if (!(splitby %in% colnames(seuratObj()@meta.data)))
    {
      splitby = NULL
    }
    
    print(paste(paste(features, collapse=","), groupby, splitby))
    obj = seuratObj()
    
    p=NULL
    
    
    if (length(features) < 2)
    {
      print("Enhanced Plots, no features")
      output$enhanced_plot = p
      return();
    }
    
    if (is.null(groupby))
    {
      print("Enhanced Plots, no group")
      output$enhanced_plot = p
      return();
    }

    if (plotMode == "Enhanced DotPlot")
    {
      
      plotElems = list()
      
      if (is.null(splitby))
      {
        plotElems[["all"]] = list(cells=colnames(obj),label="All")
      } else {
        
        metaData = obj@meta.data[ , splitby]
        names(metaData) = rownames(obj@meta.data)
        splitElems = unique(metaData)
        
        for (splitElem in splitElems)
        {
          cellnames = names(metaData[ metaData == splitElem])
          plotElems[[splitElem]] = list(cells=cellnames, label=splitElem)
        }
        
      }
      
      p=PLOSC::enhancedDotPlot(obj, plotElems,featureGenes = features, group.by = groupby, scale.by = scalingMode)
      
      
    } else if (plotMode == "Enhanced Heatmap") {
      
      
      if (!is.null(groupby))
      {
        allGroups = obj@meta.data[ , groupby]
        allGroups = unique(as.character(allGroups))
        
        plotGOIs = list()
        plotGOIs[["all"]] = list(clusters=allGroups, genes=features)
        
        
        if (is.null(splitby))
        {
          
          if (scalingMode %in% c("GLOBAL", "ALL"))
          {
            p=PLOSC::enhancedHeatMap(obj.in = obj, plot_gois = plotGOIs, group.by = groupby, scale.by = scalingMode, include_all_clusters=TRUE)
          } else {
            showNotification(paste("Scaling must be GLOBAL or ALL for this plot."), duration = NULL)
            p=NULL
          }
        } else {
          p=PLOSC::makeComplexExprHeatmapSplit(obj.in = obj, plot_gois = plotGOIs, split.by = splitby, group.by = groupby, title = "Enhanced Heat Map", scale.by = scalingMode,include_all_clusters=TRUE)
        }
      } else {
        p=NULL
      }
      
    }
    
    output$enhanced_plot = renderPlot({p})
    
  }
  
  
  
  output$fileselected<-renderText({
    paste0('You have selected: ', input$selectfile)
  })
  
  
  seuratObj <- eventReactive(input$selectfile, {
    fullpath <- file.path(fpath,input$selectfile)
    print(fullpath)
    
    obj=readRDS(fullpath)
    
    if (length(obj@reductions) == 0)
    {
      print("No Reductions in Seurat Object, preprocessing")
      obj = PLOSC::preprocessIntegrated(obj, useAssay = "RNA", with.hto = FALSE)
    }
    
    return(obj)
  })
  
  observeEvent(seuratObj(), {
    df <- seuratObj()@meta.data
    vars <- colnames(df)
    features <- rownames(seuratObj())
    # Update select input immediately after clicking on the action button.
    updateSelectInput(session, "v_groupby1","v_groupby1", choices = c("-", vars), selected="-")
    updateSelectInput(session, "v_splitby1","v_splitby1", choices = c("-", vars), selected="-")
    
    updateSelectInput(session, "v_groupby2","v_groupby2", choices = c(vars), selected=vars[1])
    updateSelectInput(session, "v_splitby2","v_splitby2", choices = c("-", vars), selected="-")
    
    updateSelectizeInput(session, "v_features2", choices=features)
  })
  
  

  
  #
  ##
  ### Main Page
  ##
  #
  output$ftxtout <- renderDataTable({
    head(seuratObj()@meta.data)
  }, options =list(pageLength = 5))
  
  output$dimplot <- renderPlot({
    DimPlot(seuratObj())
  }, res = 96)
  
  output$txtout <- renderDataTable({
    f <- seuratObj()@meta.data %>% subset(select = input$columns) 
    f$var1 <- f[[input$v_attribute1]]
    f$var2 <- f[[input$v_attribute2]]
    ff <- f %>% dplyr::filter(var1 == input$v_filter1 & var2 == input$v_filter2) 
    fff <- ff %>% subset(select=-c(var1,var2))
    head(fff)
  }, options =list(pageLength = 5)
  )
  
  #
  ##
  ### Violin Plot Page
  ##
  #
  observeEvent(input$v_splitby1, {
    makeViolinPlot(input, output)
  })
  
  observeEvent(input$v_groupby1, {
    makeViolinPlot(input, output)
  })
  
  observeEvent(input$v_plotmode1, {
    makeViolinPlot(input, output)
  })
  
  observeEvent(input$v_feature1, {
    
    if (input$v_feature1 %in% rownames(seuratObj()))
    {
      makeViolinPlot(input, output)
    } else {
      update_autocomplete_input(session, "v_feature1", options = rownames(seuratObj()))  
    }
    
  })
  
  #
  ##
  ### Enhanced Plots Page
  ##
  #
  observeEvent(input$v_splitby2, {
    makeEnhancedPlots(input, output)
  })
  
  observeEvent(input$v_groupby2, {
    makeEnhancedPlots(input, output)
  })
  
  observeEvent(input$v_plotmode2, {

    makeEnhancedPlots(input, output)
  })
  
  observeEvent(input$v_scalingmode2, {
    makeEnhancedPlots(input, output)
  })
  
  observeEvent(input$v_features2, {
    makeEnhancedPlots(input, output)
  })
  
  
  
  
  
  
  
  
  
} # server