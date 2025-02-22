#
library(shiny)
library(shinyFiles)
library(shinyWidgets)
library(shinythemes)
library(shinysky)
library(dqshiny)
library(dplyr)
library(Seurat)
library(stringr)
library(DT)

#devtools::install_github("AnalytixWare/ShinySky")
#remotes::install_github("daqana/dqshiny")

#setwd("D:/git/PLOSC")

fpath <- '../'

# remotes::install_local(".", force=TRUE, threads=4); devtools::reload(".");
# devtools::load_all("./")
# runApp("./app")


#library(executablePackeR)
#pack(
#  app_name = "PLOSC",
#  electron_settings = list(
#    c("product_name_template", "PLO(SC)2"),
#    c("app_description_template", "PLO(SC)2"),
#    c("author_name_template", "Markus Joppich"),
#    c("author_email_template", "joppich@bio.ifi.lmu.de"),
#    c("repository_url_template", "https://github.com/mjoppich/PLOSC")
#  ),
#  option = list(is_dev=TRUE)
#)

downloadHandler2 <- function(filename, content, contentType=NULL, outputArgs=list()) {
  print("downloadHandler2")
  print(content)
  if (!is.null(content))
  {
    renderFunc <- function(shinysession, name, ...) {
      shinysession$registerDownload(name, filename, contentType, content)
    }
    shiny::snapshotExclude(
      shiny::markRenderFunction(downloadButton, renderFunc, outputArgs, cacheHint = FALSE)
    )
  }
  
}

downloadablePlots = list()

# Define server function
server <- function(input, output, session) {
  
  makeViolinPlot = function(input, output)
  {
    
    if (is.null(seuratObj()))
    {
      print("Seurat Object NULL")
      return(NULL)
    }
    
    print(input$v_plotmode1)
    print(input$v_groupby1)
    print(input$v_splitby1)
    print(input$v_feature1)
    print(input$v_ymin1)
    print(input$v_ymax1)
    
    plotMode = input$v_plotmode1
    groupby = input$v_groupby1
    splitby = input$v_splitby1
    feature = input$v_feature1
    ymin = input$v_ymin1
    ymax = input$v_ymax1
    
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
        p=PLOSC::VlnBoxPlot(seuratObj(), gene = feature, group.by = groupby, split.by=splitby)
        
      } else if (plotMode == "Violin-Boxplot-Comparison")
      {
        p=PLOSC::comparativeVioBoxPlot(seuratObj(), feature = feature, group.by = groupby, split.by = splitby)
        
      } else if (plotMode == "Side-by-side Violin-Boxplot")
      {
        p=PLOSC::SplitVlnBoxPlot(seuratObj(), gene = feature, group.by = groupby, split.by=splitby)
      }
      
      
      p = p + ggplot2::ylim(c(ymin, ymax))
      
      output$vplot = renderPlot({
        p
      }, res = 96)
      return(p)
    } else {
      output$vplot = NULL
    }
    return(NULL)
  }
  
  makeEnhancedPlots = function(input, output)
  {
    
    if (is.null(seuratObj()))
    {
      print("Seurat Object NULL")
      return()
    }
    
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
      return(p);
    }
    
    if (is.null(groupby))
    {
      print("Enhanced Plots, no group")
      output$enhanced_plot = p
      return(p);
    }
    output$t_dotplot_warning = renderText("Warning: -")
    if (scalingMode == "GLOBAL")
    {
      unavailFeatures = setdiff(features, head(rownames(GetAssayData(ifnb, layer = "scale"))) )
      if (length(unavailFeatures) > 0)
      {
        output$t_dotplot_warning = renderText(paste("Warning! Features not available in scaled layer", paste(unavailFeatures, collapse=", ")))
      }
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
    return(p)
  }
  
  makeEnhancedUMAP = function(input, output)
  {
    
    if (is.null(seuratObj()))
    {
      print("Seurat Object NULL")
      return(NULL)
    }
    
    print("makeEnhancedUMAP")
    print(input$du_groupby)
    print(input$du_scalingmode)
    print(input$du_xsplit)
    print(input$du_ysplit)
    
    scalingMode = input$du_scalingmode
    
    groupby = input$du_groupby
    xsplit = input$du_xsplit
    ysplit = input$du_ysplit
    
    if (!(groupby %in% colnames(seuratObj()@meta.data)))
    {
      groupby = NULL
    }
    if (!(xsplit %in% colnames(seuratObj()@meta.data)))
    {
      xsplit = NULL
    }
    if (!(ysplit %in% colnames(seuratObj()@meta.data)))
    {
      ysplit = NULL
    }
    
    print(paste(scalingMode, groupby, xsplit, ysplit))
    obj = seuratObj()
    
    p=NULL
    
    if (is.null(groupby) || is.null(xsplit) || is.null(ysplit))
    {
      print("Enhanced Plots, no group")
      output$du_umap = p
      return(NULL);
    }
    
    downsample=FALSE
    if (scalingMode == "Downsample")
    {
      downsample=TRUE
    }
    
    p=PLOSC::makeUMAPPlot(obj, dim1 = ysplit, dim2=xsplit, group.by=groupby, downsample = downsample)
    
    downloadablePlots[["2DUMAP"]] = p
    output$du_umap = renderPlot({p})
    return(p)
    
  }
  
  output$fileselected<-renderText({
    paste0('You have selected: ', input$selectfile)
  })
  
  seuratObj = reactiveVal(NULL)
  scratchObj = reactiveVal(NULL)
  
  
  observeEvent(input$load_selected_file, {
    fullpath <- file.path(fpath,input$selectfile)
    print(fullpath)
    
    obj=readRDS(fullpath)
    
    if (length(obj@reductions) == 0)
    {
      print("No Reductions in Seurat Object, preprocessing")
      obj = PLOSC::preprocessIntegrated(obj, useAssay = "RNA", with.hto = FALSE)
    }
    
    seuratObj(obj)
  })
  
  
  
  
  
  observeEvent(seuratObj(), {
    
    obj = seuratObj()
    if (is.null(obj))
    {
      print("Seurat Object NULL")
      return()
    }
    
    df <- obj@meta.data
    vars <- colnames(df)
    features <- rownames(seuratObj())
    
    updateSelectInput(session, "u_groupby1", choices = c("-", vars), selected="-")
    
    
    # Update select input immediately after clicking on the action button.
    updateSelectInput(session, "v_groupby1", choices = c("-", vars), selected="-")
    updateSelectInput(session, "v_splitby1", choices = c("-", vars), selected="-")
    
    updateSelectInput(session, "v_groupby2", choices = c(vars), selected=vars[1])
    updateSelectInput(session, "v_splitby2", choices = c("-", vars), selected="-")
    
    updateSelectInput(session, "du_groupby", choices = c(vars), selected="-")
    updateSelectInput(session, "du_xsplit", choices = c(vars), selected="-")
    updateSelectInput(session, "du_ysplit", choices = c(vars), selected="-")
    
    update_autocomplete_input(session, "v_feature1", options = rownames(obj))
    updateSelectizeInput(session, "v_features2", choices=features)
    
    
    updateSelectInput(session, "de_groupby", choices = c(vars), selected=vars[1])
    
    output$ftxtout <- renderDT({
      head(obj@meta.data, n=1000)
    }, options =list(pageLength = 5,
                     dom = 'Bfrtip',
                     buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
    )
    )
    
    
  }, ignoreInit = TRUE)
  
  
  
  
  #
  ##
  ### Main Page
  ##
  #
  
  
  
  observeEvent(input$u_groupby1, {
    
    if (is.null(seuratObj()))
    {
      print("Seurat Object NULL")
      return()
    }
    
    groupby=input$u_groupby1
    if (!(groupby %in% colnames(seuratObj()@meta.data)))
    {
      groupby = NULL
    }
    
    p=Seurat::DimPlot(seuratObj(), group.by=groupby)
    download_dimplot(p)
    output$dimplot <- renderPlot({p})
    
  })
  
  download_dimplot = reactiveVal(NULL)
  output$download_dimplot <- downloadHandler2(
    filename = function() {
      "plot.png"
    },
    content = function(file) {
      
      p=download_dimplot()
      if (!is.null(p))
      {
        print("output width")
        print(session$clientData$output_dimplot_width)
        print("output height")
        print(session$clientData$output_dimplot_height)
        
        png(file, width=session$clientData$output_dimplot_width,
            height=session$clientData$output_dimplot_height,
            units="px", type="cairo-png")
        print(p)
        dev.off()
      } else {
        return(NULL)
      }
    }
  )
  
  
  #
  ##
  ### Violin Plot Page
  ##
  #
  observeEvent(input$v_splitby1, {
    p=makeViolinPlot(input, output)
    download_violinplot(p)
  })
  
  observeEvent(input$v_groupby1, {
    p=makeViolinPlot(input, output)
    download_violinplot(p) 
  })
  
  observeEvent(input$v_plotmode1, {
    p=makeViolinPlot(input, output)
    download_violinplot(p)
  })
  observeEvent(input$v_ymin1, {
    p=makeViolinPlot(input, output)
    download_violinplot(p)
  })
  observeEvent(input$v_ymax1, {
    p=makeViolinPlot(input, output)
    download_violinplot(p)
  })
  observeEvent(input$v_feature1, {
    p=makeViolinPlot(input, output)
    download_violinplot(p)
  })
  
  
  download_violinplot = reactiveVal(NULL)
  output$download_vplot <- downloadHandler2(
    filename = function() {
      "plot.png"
    },
    content = function(file) {
      
      p=download_violinplot()
      if (!is.null(p))
      {
        print("output width")
        print(session$clientData$output_dimplot_width)
        print("output height")
        print(session$clientData$output_dimplot_height)
        
        png(file, width=session$clientData$output_dimplot_width,
            height=session$clientData$output_dimplot_height,
            units="px", type="cairo-png")
        print(p)
        dev.off()
      } else {
        return(NULL)
      }
    }
  )
  
  #
  ##
  ### Enhanced Plots Page
  ##
  #
  observeEvent(input$v_splitby2, {
    p=makeEnhancedPlots(input, output)
    download_enhancedplot(p)
  })
  
  observeEvent(input$v_groupby2, {
    p=makeEnhancedPlots(input, output)
    download_enhancedplot(p)
  })
  
  observeEvent(input$v_plotmode2, {
    p=makeEnhancedPlots(input, output)
    download_enhancedplot(p)
  })
  
  observeEvent(input$v_scalingmode2, {
    p=makeEnhancedPlots(input, output)
    download_enhancedplot(p)
  })
  
  observeEvent(input$v_features2, {
    p=makeEnhancedPlots(input, output)
    download_enhancedplot(p)
  })
  
  download_enhancedplot = reactiveVal(NULL)
  output$download_enhanced_plot <- downloadHandler2(
    filename = function() {
      "plot.png"
    },
    content = function(file) {
      
      p=download_enhancedplot()
      if (!is.null(p))
      {
        print("output width")
        print(session$clientData$output_dimplot_width)
        print("output height")
        print(session$clientData$output_dimplot_height)
        
        png(file, width=session$clientData$output_dimplot_width,
            height=session$clientData$output_dimplot_height,
            units="px", type="cairo-png")
        print(p)
        dev.off()
      } else {
        return(NULL)
      }
    }
  )
  #
  ##
  ### 2D UMAP
  ##
  #
  observeEvent(input$du_scalingmode, {
    p=makeEnhancedUMAP(input, output)
    download_2dumap(p)
  })
  
  observeEvent(input$du_xsplit, {
    p=makeEnhancedUMAP(input, output)
    download_2dumap(p)
  })
  
  observeEvent(input$du_ysplit, {
    p=makeEnhancedUMAP(input, output)
    download_2dumap(p)
  })
  
  observeEvent(input$du_groupby, {
    p=makeEnhancedUMAP(input, output)
    download_2dumap(p)
  })
  
  
  download_2dumap = reactiveVal(NULL)
  output$download_2d_umap <- downloadHandler2(
    filename = function() {
      "plot.png"
    },
    content = function(file) {
      
      p=download_enhancedplot()
      if (!is.null(p))
      {
        print("output width")
        print(session$clientData$output_dimplot_width)
        print("output height")
        print(session$clientData$output_dimplot_height)
        
        png(file, width=session$clientData$output_dimplot_width,
            height=session$clientData$output_dimplot_height,
            units="px", type="cairo-png")
        print(p)
        dev.off()
      } else {
        return(NULL)
      }
    }
  )
  
  output$downloadData <- downloadHandler(
    filename = function() {
      "obj.integrated.Rds"
    },
    content = function(file) {
      saveRDS(scratchObj(), file)
    }
  )
  
  #
  ##
  ### SCRATCH
  ##
  #
  observeEvent(input$scratch_folder, {
    print(input$scratch_folder)
    tsvFiles = list.files(input$scratch_folder,pattern="*.tsv")
    print(tsvFiles)
    updateSelectInput(session, "scratch_annotation_samples", choices = c("-", tsvFiles), selected="-")
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      "obj.integrated.Rds"
    },
    content = function(file) {
      saveRDS(scratchObj(), file)
    }
  )
  
  observeEvent(input$scratch_go_button, {
    
    print("Here!")
    
    organism = input$scratch_organism
    
    if (organism == "Human")
    {
      print("Human")
      patternList.use = PLOSC::patternList_human()
      cc.genes.use = Seurat::cc.genes
    } else {
      print("Mouse")
      patternList.use = PLOSC::patternList_mouse()
      
      mouse.cc.genes = list()
      mouse.cc.genes[["s.genes"]] = PLOSC::convertMouseGeneList(Seurat::cc.genes$s.genes)
      mouse.cc.genes[["g2m.genes"]] = PLOSC::convertMouseGeneList(Seurat::cc.genes$g2m.genes)
      
      cc.genes.use = mouse.cc.genes
    }
    
    
    inputFolder = input$scratch_folder
    integrationFolder = paste(inputFolder, "integration", sep="/")
    
    
    # read files
    files <- Sys.glob(paste(inputFolder, "/*.h5", sep=""))
    
    print(files)
    
    sample_element = str_count(files[1], "/")+1
    inputMatrices = PLOSC::readH5Files(files, sample_element=sample_element, sample_processor=function(x){return(substr(x, 1, str_locate(x, "\\.")[1]-1  ))})
    
    # convert matrices to seurat objects
    objlist.raw = PLOSC::toObjList(inputMatrices, patternList.use, input$scratch_integration_features)
    
    # filter seurat objects for basic qc
    objlist = PLOSC::scatterAndFilter(objlist.raw,
                                      nfeature_rna.lower= input$scratch_minfeatures,
                                      nfeature_rna.upper=input$scratch_maxfeatures,
                                      ncount_rna.lower=input$scratch_minncount,
                                      percent_mt.upper=input$scratch_maxpercentmt,
                                      plot=FALSE
    )
    
    print("cells per experiment")
    print(mapply(sum, lapply(objlist, function(x) {dim(x)[2]})))
    print("total cells")
    print(sum(mapply(sum, lapply(objlist, function(x) {dim(x)[2]}))))
    
    # cleanup RAM
    objlist.raw = NULL
    inputMatrices = NULL
    
    # Prepare Integration
    print("Prepare Integration")
    finalList = PLOSC::prepareIntegration(objlist, cc.use.genes = cc.genes.use,
                                          nfeatures.variable = input$scratch_integration_features,
                                          nfeatures.scale=input$scratch_integration_features,
                                          run.parallel=FALSE)
    
    
    # Perform Integration
    print("Perform Integration")
    integratedList_sample = PLOSC::performIntegration(finalList$data, integrationFolder,
                                                      features.integration = finalList$features,
                                                      gex.method.normalization="LogNormalize",
                                                      gex.method.integration="rpca",
                                                      add.do=FALSE, run.parallel=FALSE)
    
    # Preprocess Integrated
    print("Preprocess Integrated")
    obj.integrated = PLOSC::preprocessIntegrated(integratedList_sample$integrated, "integratedgex", integrationFolder,
                                                 resolution=input$scratch_cluster_resolution, num.pcs=input$scratch_umap_dims,
                                                 dim.reduction="igpca", with.hto=FALSE)
    obj.integrated@reductions$ig.umap = NULL
    
    DefaultAssay(obj.integrated) = "RNA"
    
    dir.create(integrationFolder, showWarnings = FALSE, recursive = TRUE)
    saveRDS(obj.integrated, file.path(integrationFolder, "obj_integrated.Rds"))
    
    metaPath <- file.path(inputFolder,input$scratch_annotation_samples)
    mdf = read.csv(metaPath, header=TRUE, sep="\t")
    print(head(mdf))
    
    if ("library" %in% colnames(mdf))
    {
      metadata = dplyr::left_join(x = obj.integrated@meta.data, y = mdf, by=dplyr::join_by(library))
      print(head(metadata))
      
      obj.integrated@meta.data = metadata
      saveRDS(obj.integrated, file.path(integrationFolder, "obj_integrated.Rds"))
    }
    
    scratchObj(obj.integrated)
    seuratObj(obj.integrated)
  })
  
  
  
  observeEvent(scratchObj(), {
    
    obj = scratchObj()
    if (!is.null(obj) && !is.na(obj)) 
    {
      p = Seurat::DimPlot(obj, group.by = "idents")
      output$scratch_umap = renderPlot({p})      
    }
    
  })
  
  
  #
  ## DE Analysis
  #
  de_markers = reactiveVal(NULL)
  
  observeEvent(input$de_process, {
    
    obj = seuratObj()
    if (!is.null(obj) && !is.na(obj)) 
    {
      markers = PLOSC::makeDEResults(obj, group.by=input$de_groupby, assay="RNA", test=input$de_method)
      de_markers(markers)
      
      output$de_out <- renderDT({
        de_markers()
      }, options =list(pageLength = 5,
                       dom = 'Bfrtip',
                       buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
      )
      )
      
      showNotification("Differential Gene Expression analysis finished. Results available.")
      shinyjs::show("de_download")
      
    }
    
  })
  output$de_download <- downloadHandler(
    filename = function() {
      "de_markers.tsv"
    },
    content = function(file) {
      write.table(x=de_markers(), file = file, sep="\t", row.names = FALSE)
    }
  )
  
  #
  ## GSE Analysis
  #
  
  observeEvent(input$gse_folder, {
    print(input$gse_folder)
    tsvFiles = list.files(input$gse_folder,pattern="*.tsv")
    print(tsvFiles)
    updateSelectInput(session, "gse_de_file", choices = c("-", tsvFiles), selected="-")
  })
  
  observeEvent(input$gse_de_file, {
    
    if (input$gse_de_file != "-")
    {
      gseFolder = input$gse_folder
      defile = file.path(gseFolder, input$gse_de_file)
      deDF = read.table(defile, sep = "\t", header = TRUE)
      gse_de_tsv_file(deDF)
      
      output$gse_inputfile <- renderDT({
        gse_de_tsv_file()
      }, options =list(pageLength = 5,
                       dom = 'Bfrtip',
                       buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
      )
      )
    }
  })
  
  gse_markers = reactiveVal(NULL)
  gse_de_tsv_file = reactiveVal(NULL)
  
  
  observeEvent(input$gse_process, {
    
    obj = seuratObj()
    if (is.null(obj))
    {
      showNotification("Please load Seurat Object first")
      return()
    }
    
    gseFolder = input$gse_folder
    
    deDF = gse_de_tsv_file()
    
    allClusterIDs = unique(deDF$clusterID)
    deGroups = list()
    
    print(colnames(deDF))
    
    for (cid in allClusterIDs)
    {
      deGroups[[cid]] = deDF[deDF$clusterID == cid,]
    }
    
    if (input$gse_plot_volcanos)
    {
      PLOSC::makeVolcanos(deGroups, "DE Comparison", file.path(gseFolder, "volcano/volcano"), FCcutoff=input$gse_min_foldchange)  
    }
    
    print("Do Enrichment Analysis")
    gseResults = PLOSC::performEnrichtmentAnalysis(seuratObj(), deGroups, input$gse_organism, file.path(gseFolder, "enrichment/gse_results.rds"))
    
    print("Plot Enrichment Analysis")
    PLOSC::makeEnrichmentPlots( gseResults, outfolder = gseFolder )
    print("Enrichment Analysis Done")
    
    gse_markers(gseResults)
    
    showNotification("Gene Set Enrichment processing finished. Results available.")
    
  })
  
  
  
} # server