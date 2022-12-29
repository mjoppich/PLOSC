

#
##
### Global Helpers
##
#

patternList.human = list()
patternList.human[["MT"]] = "^MT-"
patternList.human[["RPL"]] = "^RPL"
patternList.human[["RPS"]] = "^RPS"


patternList.mouse = list()
patternList.mouse[["MT"]] = "^mt-"
patternList.mouse[["RPL"]] = "^Rpl"
patternList.mouse[["RPS"]] = "^Rps"



#
##
### Seurat Read-In
##
#

#' Reads in a list of mtx files and returns two named lists of gene expression counts and Antibody Capture Counts
#'
#'
#' @param files list of mtx paths
#'
#' @return list of gex and ab (gene expression, antibody capture) lists of per cell counts
#'
#'
#' @export
readMtxFiles = function(files)
{

  allfiles.raw = list()
  allABs.raw = list()

  for (file in files)
  {
    samplename = stringr::str_split(dirname(file), "/")[[1]][3]
    foldername = dirname(file)
    
    print(paste(samplename, foldername))
    
    h5file = Seurat::Read10X(foldername,unique.features = TRUE)

    if (is.null(names(h5file)))
    {
      print(paste("WITHOUT AB", samplename))
      allfiles.raw[[samplename]] = h5file
    } else {
        print(paste("WITH AB", samplename))
      allfiles.raw[[samplename]] = h5file$`Gene Expression`
      allABs.raw[[samplename]] = h5file$`Antibody Capture`
    }

    print(paste(samplename, ncol(allfiles.raw[[samplename]]), "genes x cells"))
  }

  return(list(gex=allfiles.raw, ab=allABs.raw))
}

#' Turns a feature-cell-matrix into a Seurat objects with called mt/rp/rps/rpl-content
#'
#'
#' @param matrix feature-cell-matrix of the sample
#' @param proj project name for the seurat object
#' @param pl patternlist for mt-content and rp/rps/rpl-content
#'
#' @return Seurat object
#'
#' @export
makeSeuratObj = function(matrix, proj, pl)
{
    obj = Seurat::CreateSeuratObject(matrix, project=proj)
    print("Renaming Cells")
    obj <- Seurat::RenameCells(obj, add.cell.id=proj)
    
    print(paste("Seurat obj project", obj@project.name))
    
    mtPattern = pl[["MT"]]
    rplPattern = pl[["RPL"]]
    rpsPattern = pl[["RPS"]]
    rpPattern = paste(c(pl[["RPL"]], pl[["RPS"]]), sep="", collapse="|")

    
    selGenes = rownames(obj)[grepl(rownames(obj), pattern=mtPattern)]
    print(paste("Got a total of mt-Genes:", length(selGenes), paste(head(selGenes), collapse=", ")))
    
    selGenes = rownames(obj)[grepl(rownames(obj), pattern=rplPattern)]
    print(paste("Got a total of Rpl-Genes:", length(selGenes), paste(head(selGenes), collapse=", ")))
    
    selGenes = rownames(obj)[grepl(rownames(obj), pattern=rpsPattern)]
    print(paste("Got a total of Rps-Genes:", length(selGenes), paste(head(selGenes), collapse=", ")))
    
    selGenes = rownames(obj)[grepl(rownames(obj), pattern=rpPattern)]
    print(paste("Got a total of Rp-Genes:", length(selGenes), paste(head(selGenes), collapse =", ")))
    
    obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = mtPattern)
    obj[["percent.rpl"]] <- Seurat::PercentageFeatureSet(obj, pattern = rplPattern)
    obj[["percent.rps"]] <- Seurat::PercentageFeatureSet(obj, pattern = rpsPattern)
    obj[["percent.rp"]] <- Seurat::PercentageFeatureSet(obj, pattern = rpPattern)
    
    return(obj)
}


#' Transforms a list of gene expression matrices into Seurat objects, normalizes data and calculate variable features
#'
#'
#' @param inputMatrices list of count matrices to transform into Seurat object
#' @param patternlist 
#' @param variable.features list of mtx paths
#'
#' @return list of gex and ab (gene expression, antibody capture) lists of per cell counts
#'
#'
#' @export
toObjList = function(inputMatrices, patternlist, variable.features=3000)
{

objlist = list()

for (x in names(inputMatrices$gex))
{

    matrix = inputMatrices$gex[[x]]
    
    filteredObj = makeSeuratObj(matrix, x, patternlist)   
    
    filteredObj <- Seurat::NormalizeData(filteredObj, verbose = FALSE)
    filteredObj <- Seurat::FindVariableFeatures(filteredObj, nfeatures=variable.features, verbose = FALSE)
    
    objlist[[x]] = filteredObj

    print(x)
    print(filteredObj)
        
}

return(objlist)

}



#' Takes a list of seurat objects and creates QC plots and filters all objects according to nFeature_RNA, nCount_RNA and percent.mt
#'
#'
#' @param objlist list of Seurat objects
#' @param nfeature_rna.lower lower bound of nFeature_RNA
#' @param nfeature_rna.upper upper bound of nFeature_RNA
#' @param ncount_rna.lower lower bound of nCount_RNA
#' @param percent_mt.upper upper bound of percent.mt
#'
#' @return filtered objlist
#'
#'
#' @export
scatterAndFilter = function(objlist, nfeature_rna.lower=100, nfeature_rna.upper=6000, ncount_rna.lower=500, percent_mt.upper=7)
{

  for (name in names(objlist))
  {
    print(name)

    plot1 <- Seurat::FeatureScatter(objlist[[name]], feature1 = "nCount_RNA", feature2 = "percent.mt")
    plot2 <- Seurat::FeatureScatter(objlist[[name]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    save_plot(plot1 + plot2, paste(name, "scatter_ncount_mt", sep="_"), fig.width=10, fig.height=6)

    plot1 <- Seurat::FeatureScatter(objlist[[name]], feature1 = "nCount_RNA", feature2 = "percent.rp")
    plot2 <- Seurat::FeatureScatter(objlist[[name]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    save_plot(plot1 + plot2, paste(name, "scatter_ncount_rp", sep="_"), fig.width=10, fig.height=6)
  }


  objlist.new <- lapply(X = objlist, FUN = function(obj) {
    # mt content: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6072887/
    
    obj <- subset(obj, subset = nFeature_RNA > nfeature_rna.lower & nFeature_RNA < nfeature_rna.upper & nCount_RNA > ncount_rna.lower)
    obj <- subset(obj, subset = percent.mt < percent_mt.upper)
    print(obj)
    
    return(obj)
  })


  for (name in names(objlist.new))
  {
    p=Seurat::VlnPlot(objlist.new[[name]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
    save_plot(p, paste(name, "filtered_violins_qc", sep="_"), fig.width=10, fig.height=6)
    
    plot1 <- Seurat::FeatureScatter(objlist.new[[name]], feature1 = "nCount_RNA", feature2 = "percent.mt")
    plot2 <- Seurat::FeatureScatter(objlist.new[[name]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    save_plot(plot1 + plot2, paste(name, "filtered_scatter_ncount_mt", sep="_"), fig.width=10, fig.height=6)
    
    plot1 <- Seurat::FeatureScatter(objlist.new[[name]], feature1 = "nCount_RNA", feature2 = "percent.rp")
    plot2 <- Seurat::FeatureScatter(objlist.new[[name]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    save_plot(plot1 + plot2, paste(name, "filtered_scatter_ncount_rp", sep="_"), fig.width=10, fig.height=6)
  }

  return(objlist.new)

}



#
##
### Hashtag Oligos
##
#


#' Takes a list of antibody capture feature-cell-matrices, the list of Seurat objects and a list of relevant hashtag oligo IDs for each sample
#'
#'
#' @param inputMatrices list of antibody capture feature-cell-matrices
#' @param objlist list of Seurat objects
#' @param relevantHTOs relevant hashtag oligo IDs for each sample
#'
#' @return list of Seurat objects with HTO assay
#'
#'
#' @export
processHTO = function(inputMatrices, objlist, relevantHTOs)
{

htoObjList = list()
for (name in intersect(names(inputMatrices$ab), names(relevantHTOs)))
{

    relHTOs = relevantHTOs[[name]]
    print("Relevant HTOs:")
    print(relHTOs)

    htoMatrix = inputMatrices$ab[[name]]

    cellsGEX = length(colnames(objlist[[name]]))
    cellsAB = length(colnames(htoMatrix))

    gexnames = substring(colnames(objlist[[name]]), str_length(name)+2)
    print(head(gexnames))
    cellsJoint = intersect(gexnames, colnames(htoMatrix))
    cellsABsub = htoMatrix[, cellsJoint]
    cellsABsub = cellsABsub[relHTOs,]

    colnames(cellsABsub) = paste(objlist[[name]]@project.name, colnames(cellsABsub), sep="_")
    print(head(colnames(cellsABsub)))
    print(paste(cellsGEX, cellsAB, length(cellsJoint), ncol(cellsABsub)))

    # Normalize RNA data with log normalization
    xobj <- Seurat::NormalizeData(objlist[[name]])
    # Find and scale variable features
    xobj <- Seurat::FindVariableFeatures(xobj, selection.method = "mean.var.plot")
    xobj <- Seurat::ScaleData(xobj, features = VariableFeatures(xobj))

    xobj[["HTO"]] <- Seurat::CreateAssayObject(counts = cellsABsub)
    xobj <- Seurat::NormalizeData(xobj, assay = "HTO", normalization.method = "CLR")
    xobj <- Seurat::HTODemux(xobj, assay = "HTO", positive.quantile = 0.99)

    print(table(xobj$HTO_classification.global))
    print(table(xobj$HTO_classification))

    htoObjList[[name]] = xobj

}

return(htoObjList)

}


#' Plot quality control plots for the Seurat objects with HTO assay.
#'
#'
#' @param objlist list of Seurat objects with HTO assay
#'
#' @return list of Seurat objects with HTO assay
#'
#'
#' @export
qcFilterHTO = function(htoObjList)
{

  for (name in names(htoObjList))
  {
    p=Seurat::VlnPlot(htoObjList[[name]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0, group.by = "HTO_classification.global")
    save_plot(p, paste(name, "hto_violins_qc", sep="_"), fig.width=10, fig.height=6)


    allHTOFeatures = rownames(htoObjList[[name]][["HTO"]])
    figHeight = round(length(allHTOFeatures)*0.5 * 6)
    r=Seurat::RidgePlot(htoObjList[[name]], assay = "HTO", features = allHTOFeatures, ncol = 2)
    save_plot(r, paste(name, "hto_ridge_qc", sep="_"), fig.width=10, fig.height=figHeight)

    allHTOFeatures = rownames(htoObjList[[name]][["HTO"]])
    figHeight = round(length(unique(htoObjList[[name]]$HTO_classification))*0.5 * 4)
    r=Seurat::RidgePlot(htoObjList[[name]], assay = "HTO", features = allHTOFeatures, group.by="HTO_classification", ncol = 2)
    save_plot(r, paste(name, "hto_ridge_detail_qc", sep="_"), fig.width=10, fig.height=figHeight)

  
    print(name)
    print(table(htoObjList[[name]]$HTO_classification))

  }

  htoFiltered <- lapply(X = htoObjList, FUN = function(obj) {
    
    print(obj)
    print(table(obj$HTO_classification))
    obj = subset(obj, subset=HTO_classification.global == "Singlet")
    print(obj)
    print(table(obj$HTO_classification)) 
    return(obj)
  })

  return(htoFiltered)

}



#' Merges the list with Seurat objects and the list with HTO Seurat objects
#'
#'
#' @param objlist list of Seurat objects
#' @param htolist list of Seurat objects with HTO assay
#'
#' @return list of Seurat objects
#'
#' @export
mergeRNAHTO = function(objlist, htolist)
{

  returnList = list()

  for (name in names(objlist))
  {
    if (name %in% names(htolist))
    {
      returnList[[name]] = htolist[[name]]
    } else {
      returnList[[name]] = objlist[[name]]
    }
  }
  return(returnList)
}


#' Splits all Seurat objects in the input list by a specific meta tag; useful for sample integration
#'
#'
#' @param objlist list of Seurat objects
#' @param group.by meta-column on which values the Seurat object should be split
#'
#' @return list of Seurat objects
#'
#' @export
splitObjListByGroup = function(objlist, group.by)
{
  finalList = list()
  for (objname in names(objlist))
  {

    if (group.by %in% colnames(objlist[[objname]]@meta.data))
    {
        for (abTag in unique(objlist[[objname]][[group.by]][[group.by]]))
        {
            sampleName = paste(objname, abTag, sep="_")
            print(sampleName)

            selDF = objlist[[objname]][[group.by]]
            selCells = rownames(selDF)[selDF[[group.by]] == abTag]

            finalList[[sampleName]] = subset(objlist[[objname]], cells = selCells)
            finalList[[sampleName]]$orig_project = objname
            finalList[[sampleName]]$library = paste(objname, abTag, sep="_")
            
            print(sampleName)
            print(finalList[[sampleName]])
        }

    } else {
        finalList[[objname]] = objlist[[objname]]

        print(objname)
        print(finalList[[objname]])
    }

  }

  return(finalList)
}




#
##
###
#### Seurat object integration
###
##
#


#' This function prepares an integration on the Seurat objects by running NormalizeData, FindVariableFeatures per objects, SelectIntegrationFeatures on all, and then scoring cell cycle genes, running ScaleData and PCA again on all objects.
#'
#'
#' @param finalList list of Seurat objects
#' @param cc.use.genes list of s.genes and g2m.genes
#' @param nfeatures.variable number of variable features to detect
#' @param nfeatures.scale number of features to scale
#' @param scale.regress which variables to regress out during ScaleData
#' @param normalize whether to normalize the data
#' @param findvariable whether to run FindVariableFeatures
#' @param run.parallel whether scale data should run in parallel. if FALSE, ScaleData and RunPCA will be run in a sequential context.
#'
#' @return list of data=list of Seurat objects and features=features determined by SelectIntegrationFeatures
#'
#'
#' @export
prepareIntegration = function(finalList, cc.use.genes, nfeatures.variable = 3000, nfeatures.scale=3000, scale.regress=c('percent.rp', 'percent.mt', "nCount_RNA","S.Score", "G2M.Score"), normalize=TRUE, findvariable=TRUE, run.parallel=TRUE)
{
 

    print("cells per experiment")
    print(mapply(sum, lapply(finalList, function(x) {dim(x)[2]})))
    print("total cells")
    print(sum(mapply(sum, lapply(finalList, function(x) {dim(x)[2]}))))


    objlist = list()
    for (objname in names(finalList))
    {
        x = finalList[[objname]]

        if (! "orig_project" %in% colnames(x@meta.data))
        {
          x$orig_project = objname
        }

        Seurat::Project(x) = objname
        print(paste("Seurat obj project", x@project.name))

        Seurat::DefaultAssay(x) = "RNA"

        if (normalize)
        {
          x <- Seurat::NormalizeData(x, verbose = FALSE)
        }
        
        if (findvariable)
        {
          x <- Seurat::FindVariableFeatures(x, nfeatures=nfeatures.variable, verbose = FALSE)
        }
      
        x$library = objname

        objlist[[objname]] = x
    }

    print("SelectIntegrationFeatures")
    features <- Seurat::SelectIntegrationFeatures(object.list = objlist, nfeatures = nfeatures.scale)

    objlist <- lapply(X = objlist, FUN = function(x) {

        print(paste("Seurat obj project", x@project.name))
        print(x)

        s.genes <- cc.use.genes$s.genes
        g2m.genes <- cc.use.genes$g2m.genes

        s.genes = intersect(s.genes, rownames(x))
        g2m.genes = intersect(g2m.genes, rownames(x))

        print(paste("CellCycle", length(s.genes), length(g2m.genes)))

        x <- Seurat::CellCycleScoring(
        x,
        g2m.features = g2m.genes,
        s.features = s.genes)

        if (!run.parallel)
        {
        t = future::plan()
        future::plan("sequential")
        }

        x <- Seurat::ScaleData(x, features = features, verbose = TRUE, assay="RNA", vars.to.regress = scale.regress)
        x <- Seurat::RunPCA(x, verbose = FALSE, reduction.name="pca", assay="RNA")

        if (!run.parallel)
        {
        future::plan(t)
        }

        x$project = x@project.name

        return(x)
    })

    return(list("data"=objlist, "features"=features))
}






#' This method performs the integration of all Seurat objects in `objlist` into a single Seurat object
#'
#' @param objlist list of Seurat objects
#' @param intname output path for all plots
#' @param features.integration Features to perform the gene expression integration on
#' @param gex.method.normalization Which normalization strategy to use for normalization (SCT or LogNormalize)
#' @param gex.assay assay to use for integration
#' @param gex.runpca whether to run PCA on the gene expression assay prior to integration
#' @param gex.dims number of dimensions for finding gene expression integration
#' @param gex.method.integration Which integration method to use (SCTransform, cca or rpca) for the gene expression data
#' @param gex.k.filter How many neighbors (k) to use when filtering anchors for the gene expression data
#' @param gex.k.anchor How many neighbors (k) to use when picking anchors for the gene expression data
#' @param gex.k.weight Number of neighbors to consider when weighting anchors for the gene expression data
#' @param add.do Whether an additional assay should be used for integration (e.g. CITE)
#' @param add.assay Which assay to use as additional assay
#' @param add.dims number of dimensions for finding additional integration
#' @param add.method.integration Which integration method to use (SCTransform, cca or rpca) for the additional data
#' @param add.k.filter How many neighbors (k) to use when filtering anchors for the additional data
#' @param add.k.anchor How many neighbors (k) to use when picking anchors for the additional data
#' @param add.features.integration Features to perform the additional data integration on
#' @param add.k.weight Number of neighbors to consider when weighting anchors for the additional data
#' @param run.parallel whether to run ScaleData sequentially
#'
#' @return list(integrated, multimodal) where integrated is the gene expression only integrated version, and multimodal the gex+add version
#' @export
#'
performIntegration = function(objlist, intname, features.integration = 3000, 
 gex.method.normalization="LogNormalize", gex.assay="RNA", gex.runpca=TRUE,
gex.dims=30, gex.method.integration="rpca", gex.k.filter=200, gex.k.anchor = 5,gex.k.weight=100, 
add.do=FALSE,add.assay="ABA",
add.dims=10, add.method.integration="rpca", add.k.filter=200, add.k.anchor = 5, add.features.integration=100, add.k.weight=5,
run.parallel=TRUE) 
{

     dir.create(intname, recursive = TRUE)

    if (add.do)
    {
      #
      # integrate based on ADDitional assay
      #
      objSamples = objlist

      objSamples = lapply(objSamples, function(x) {
        print(paste("Object", x@project.name))
        Seurat::DefaultAssay(x) <- add.assay
        Seurat::VariableFeatures(x) <- rownames(x) # all HTOs

        if (dim(x@assays[[add.assay]]@scale.data)[1] == 0)
        {
          if (!run.parallel)
          {
            t = future::plan()
            future::plan("sequential")
          }
          x = Seurat::ScaleData(x, assay=add.assay)
          if (!run.parallel)
          {
            future::plan(t)
          }

        }
        
        if (("pca" %in% names(x@reductions)) && (x@reductions$pca@assay.used == add.assay))
        {
          print("PCA ALREADY THERE")
        } else {
          x = Seurat::RunPCA(x,features = rownames(x),verbose = FALSE, reduction.name="pca", approx=FALSE, npcs=add.dims, assay=add.assay)
        }

        return(x)
      })

      print(objSamples)

      features_add <- rownames(objSamples[[1]][[add.assay]])

      if (!run.parallel)
      {
        t = future::plan()
        future::plan("sequential")
      }
      assayData = rep(add.assay, length(objSamples))
      objlist.anchors.add <- Seurat::FindIntegrationAnchors(object.list = objSamples, assay=assayData, normalization.method = "LogNormalize",
                                                    anchor.features = add.features.integration, dims = 1:add.dims, reduction = add.method.integration, k.filter = add.k.filter, k.anchor=add.k.anchor)

      add.list.integrated <- Seurat::IntegrateData(new.assay.name = "integrated_add", anchorset = objlist.anchors.add, normalization.method = "LogNormalize", dims=1:(add.dims-1), k.weight=add.k.weight)

      if (!run.parallel)
      {
        future::plan(t)
      }


      print("ADD integration done")
    }


    #
    # integrate based on RNA/GEX assay
    #
    objSamples = objlist
    print(objSamples)

    print("GEX integration features")
    print(objSamples)
    print(paste("Current integration mode:", gex.method.integration))

    objSamples = lapply(objSamples, function(x) {
        print(paste("Object", x@project.name))
        Seurat::DefaultAssay(x) <- gex.assay
        return(x)
    })


    if (gex.method.integration!="SCT")
    {


      if (gex.method.normalization == "SCT")
      {
        print("SCTransform")
        objSamples <- lapply(X = objSamples, FUN = sctransform::SCTransform, method = "glmGamPoi")

      }

      print("SelectIntegrationFeatures")
      if (is.numeric(features.integration))
      {
        print("RunPCA")
        objSamples <- lapply(X = objSamples, FUN = Seurat::RunPCA, npcs=min(c(50, gex.dims)), verbose = FALSE, reduction.name="pca", assay=gex.assay)
        print("Select Integration Features")
        features_gex <- Seurat::SelectIntegrationFeatures(object.list = objSamples, nfeatures = features.integration, assay=rep(gex.assay, length(objSamples)))
      } else {
        features_gex = features.integration
        print("RunPCA on given features")

        objSamples <- lapply(X = objSamples, FUN = function(x)
        {
          if (gex.runpca)
          {
            x = Seurat::RunPCA(x,npcs=min(c(50, gex.dims)),verbose = FALSE, reduction.name="pca", features=features_gex, assay=gex.assay)
          }

          return(x)
      })

      }


      if (gex.method.normalization == "SCT")
      {
        print("PrepSCTIntegration")
        objSamples <- Seurat::PrepSCTIntegration(object.list = objSamples, anchor.features = features_gex)
        print("Calculating PCAs on SCT")
        objSamples <- lapply(X = objSamples, FUN = Seurat::RunPCA, features = features_gex)
      }

      print("FindIntegrationAnchors")

      if (!run.parallel)
      {
        t = future::plan()
        future::plan("sequential")
      }
      objlist.anchors <- Seurat::FindIntegrationAnchors(object.list = objSamples,  reduction = gex.method.integration, dims = 1:gex.dims, anchor.features = features_gex, normalization.method=gex.method.normalization, k.anchor=gex.k.anchor, k.filter = gex.k.filter)
      obj.list.integrated <- Seurat::IntegrateData(new.assay.name = "integrated_gex", anchorset = objlist.anchors, dims = 1:gex.dims, verbose=T, normalization.method = gex.method.normalization, k.weight=gex.k.weight)
      
      if (!run.parallel)
      {
        future::plan(t)
      }


      print("IntegrateData")


    } else {

      objSamples = lapply(objSamples, function(x) {
          Seurat::DefaultAssay(x) <- gex.assay

          x <- Seurat::RunPCA(x, npcs=max(c(50, gex.dims)), verbose = FALSE, reduction.name="pca",  assay=gex.assay)
          suppressWarnings(x <- sctransform::SCTransform(x,vars.to.regress = c('percent.mt', 'percent.rp'), verbose = T))

          return(x)
      })

      if (!run.parallel)
      {
        t = future::plan()
        future::plan("sequential")
      }

      features_gex <- Seurat::SelectIntegrationFeatures(object.list = objSamples, nfeatures = features.integration)#, assay=rep("RNA", length(objSamples)))
      objSamples <- Seurat::PrepSCTIntegration(object.list = objSamples, anchor.features = features_gex)

      objlist.anchors <- Seurat::FindIntegrationAnchors(object.list = objSamples, normalization.method = "SCT", anchor.features = features_gex, k.anchor=gex.k.anchor,k.filter = add.k.filter)
      obj.list.integrated <- Seurat::IntegrateData(anchorset = objlist.anchors, normalization.method = "SCT", new.assay.name = "integrated_gex",verbose=T, k.weight=gex.k.weight)

      if (!run.parallel)
      {
        future::plan(t)
      }

    }
    print("GEX integration done")

    #
    # integrated GEX viz
    #
    if (gex.method.normalization != "SCT")
    {
      if (!run.parallel)
      {
        t = future::plan()
        future::plan("sequential")
      }
      obj.list.integrated = Seurat::ScaleData(obj.list.integrated, assay="integrated_gex")
    
      if (!run.parallel)
      {
        future::plan(t)
      }
    }
    obj.list.integrated <- Seurat::RunPCA(obj.list.integrated, npcs = gex.dims, reduction.name="igpca", assay="integrated_gex")
    obj.list.integrated <- Seurat::RunUMAP(obj.list.integrated, reduction = "igpca", dims = 1:gex.dims, reduction.name="ig.umap", reduction.key = "UMAPig_",)
    p=Seurat::DimPlot(obj.list.integrated, group.by="orig_project", reduction="ig.umap", shuffle = TRUE, seed = 1)
    save_plot(p, paste(intname, "ig_dimplot", sep="/"), 8, 6)

    obj.list.gex_add = NULL

    if (add.do)
    {
      #
      # integrated ADT viz
      #
      obj.list.integrated[["integrated_add"]] = add.list.integrated[["integrated_add"]]

      if (!run.parallel)
      {
        t = future::plan()
        future::plan("sequential")
      }

      obj.list.integrated = Seurat::ScaleData(obj.list.integrated, assay="integrated_add")

      if (!run.parallel)
      {
        future::plan(t)
      }

      obj.list.integrated <- Seurat::RunPCA(obj.list.integrated, features = rownames(add.list.integrated[[add.assay]]), verbose = FALSE, approx=FALSE, npcs=add.dims, reduction.name="iapca", assay="integrated_add")
      obj.list.integrated <- Seurat::RunUMAP(obj.list.integrated, reduction = "iapca", dims = 1:add.dims, reduction.name="ia.umap", reduction.key = "UMAPia_",)

      p=Seurat::DimPlot(obj.list.integrated, group.by="orig_project", reduction="ia.umap", shuffle = TRUE, seed = 1)
      save_plot(p, paste(intname, "ia_dimplot", sep="/"), 8, 6)


      #
      # multi modal neighbors
      #

      obj.list.gex_add <- Seurat::FindMultiModalNeighbors(
        obj.list.integrated, reduction.list = list("igpca", "iapca"), 
        dims.list = list(1:gex.dims, 1:add.dims), prune.SNN=1/20
      )
      #
      # multi modal viz
      #
      obj.list.gex_add <- Seurat::RunUMAP(obj.list.gex_add, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
      obj.list.gex_add <- Seurat::FindClusters(obj.list.gex_add, graph.name = "wsnn", algorithm = 3, resolution = 1, verbose = FALSE)

      p <- Seurat::DimPlot(obj.list.gex_add, reduction = 'wnn.umap', label = TRUE, repel = FALSE, label.size = 2.5)
      save_plot(p, paste(intname, "wnn_cluster_dimplot", sep="/"), 8, 6)

      p <- Seurat::DimPlot(obj.list.gex_add, group.by="orig_project", reduction = 'wnn.umap', label = TRUE, repel = FALSE, label.size = 2.5)
      save_plot(p, paste(intname, "wnn_project_dimplot", sep="/"), 8, 6)

      p <- Seurat::DimPlot(obj.list.gex_add, reduction = 'ig.umap', label = TRUE, repel = FALSE, label.size = 2.5)
      save_plot(p, paste(intname, "wnn_cluster_ig_dimplot", sep="/"), 8, 6)

      p <- Seurat::DimPlot(obj.list.gex_add, reduction = 'ia.umap', label = TRUE, repel = FALSE, label.size = 2.5)
      save_plot(p, paste(intname, "wnn_cluster_ia_dimplot", sep="/"), 8, 6)


      obj.list.integrated$wnn_clusters = Seurat::Idents(obj.list.gex_add)
    }

    return(list("integrated"=obj.list.integrated, "multimodal"=obj.list.gex_add))

}








#' Preprocesses an integrated Seurat object for downstream analysis
#'
#' @param obj.in Seurat object
#' @param useAssay which assay to use
#' @param inname output folder of all plots
#' @param do.scale whether to scale the object
#' @param num.pcs number of PCs to calculate and use for UMAP/neighbors
#' @param resolution resolution for the clustering
#' @param plot.reduction which reduction to use for plotting
#' @param dim.reduction name of the reduction in which the PCA will be stored, and which is used for UMAP and Neighbors
#' @param with.hto whether also HTO plots should be prepared
#' @param run.parallel whether the ScaleData function should run in parallel or sequential
#'
#' @return preprocessed Seurat object
#' @export
#'
preprocessIntegrated = function(obj.in, useAssay, inname, do.scale=T, num.pcs=50, resolution=0.5, plot.reduction="umap", dim.reduction="umap", with.hto=TRUE, run.parallel=TRUE)
{
  
  if (!dir.exists(inname))
  {
    print(paste("Creating DIR", inname))
    dir.create(inname, recursive = TRUE)
  }


  Seurat::DefaultAssay(obj.in) <- useAssay

  # Run the standard workflow for visualization and clustering
  if (do.scale)
  {
    # Scale is not required for SCT!
    print("Scale Data")

    if (!run.parallel)
    {
      t = future::plan()
      future::plan("sequential")
    }
    obj.in <- Seurat::ScaleData(obj.in, verbose = FALSE)

    if (!run.parallel)
    {
      future::plan(t)
    }
  }
    
  if ((is.null(dim.reduction)) || (!dim.reduction %in% names(obj.in@reductions)))
  {
    print("RunPCA Data")
    dim.reduction = "pca"

    obj.in <- Seurat::RunPCA(obj.in, npcs = max(c(num.pcs, 50)), verbose = FALSE, reduction.name=dim.reduction)
    
  }

  print(paste("dim.reduction", dim.reduction))

  p=Seurat::ElbowPlot(obj.in, ndims=30, reduction = dim.reduction)
  save_plot(p, paste(inname, "elbowplot", sep="/"), 12, 6)

  
  print("RunUMAP Data")
  obj.in <- Seurat::RunUMAP(obj.in, reduction = dim.reduction, dims = 1:num.pcs)
  print("FindNeighbors Data")
  obj.in <- Seurat::FindNeighbors(obj.in, reduction = dim.reduction, dims = 1:num.pcs)
  print("FindClusters Data")
  obj.in <- Seurat::FindClusters(obj.in, resolution = resolution)


  obj.in$idents = Seurat::Idents(obj.in)

  p=Seurat::DimPlot(obj.in, pt.size = 0.001, label=T, reduction = plot.reduction)
  save_plot(p, paste(inname, "dimplot_umap", sep="/"), fig.width=12, fig.height=8)

  numProjects = length(unique(obj.in$orig_project))
  numRows = ceiling(numProjects/2)

  p=Seurat::DimPlot(obj.in, pt.size = 0.001, label=T, split.by="orig_project", reduction = plot.reduction, ncol=2)
  save_plot(p, paste(inname, "dimplot_umap_project", sep="/"), fig.width=24, fig.height=8*numRows)
  

  if (with.hto)
  {
    tryCatch({

      obj.in$libraryHTO = paste(obj.in$library, obj.in$HTO_classification, sep="_")

      numLibHTO = length(unique(obj.in$libraryHTO))
      numRows = ceiling(numLibHTO/3)

      p=Seurat::DimPlot(obj.in, split.by="libraryHTO", pt.size = 0.001, label=T, reduction = plot.reduction, ncol=3)
      save_plot(p, paste(inname, "dimplot_umap_libraryHTO", sep="/"), fig.width=16, fig.height=min(3*numRows, 60))
    })
  }


  return(obj.in)

}











#' Creates QC Plots for a Seurat objects and a given reduction
#'
#' @param inobj Seurat object
#' @param outfolder Output path/stump for the plots
#' @param reduction using which reduction for the plots
#'
#' @return Seurat object
#' @export
#'
makeQCPlots = function(inobj, outfolder, reduction="umap")
{
  
  p=Seurat::FeaturePlot(inobj, "nCount_RNA", reduction=reduction)
  save_plot(p, paste(outfolder, "fplot_ncount_rna", sep="/"), fig.width=8, fig.height=6)

  p=Seurat::FeaturePlot(inobj, "nFeature_RNA", reduction=reduction)
  save_plot(p, paste(outfolder, "fplot_nfeature_rna", sep="/"), fig.width=8, fig.height=6)

  inobj$log_ncount = log(inobj$nCount_RNA)
  inobj$log_nfeature = log(inobj$nFeature_RNA)

  p=Seurat::FeaturePlot(inobj, "log_ncount", reduction=reduction)
  save_plot(p, paste(outfolder, "fplot_logncount_rna", sep="/"), fig.width=8, fig.height=6)

  p=Seurat::FeaturePlot(inobj, "log_nfeature", reduction=reduction)
  save_plot(p, paste(outfolder, "fplot_lognfeature_rna", sep="/"), fig.width=8, fig.height=6)

  p=Seurat::FeaturePlot(inobj, "percent.mt", reduction=reduction)
  save_plot(p, paste(outfolder, "fplot_percent_mt", sep="/"), fig.width=8, fig.height=6)

  p=Seurat::FeaturePlot(inobj, "percent.rp", reduction=reduction)
  save_plot(p, paste(outfolder, "fplot_percent_rp", sep="/"), fig.width=8, fig.height=6)

  p=Seurat::FeaturePlot(inobj, "G2M.Score", reduction=reduction)
  save_plot(p, paste(outfolder, "fplot_g2mscore", sep="/"), fig.width=8, fig.height=6)

  p=Seurat::FeaturePlot(inobj, "S.Score", reduction=reduction)
  save_plot(p, paste(outfolder, "fplot_sscore", sep="/"), fig.width=8, fig.height=6)




  p=Seurat::VlnPlot(inobj, "log_ncount", group.by="idents")
  save_plot(p, paste(outfolder, "vplot_logncount_rna", sep="/"), fig.width=12, fig.height=4)

  p=Seurat::VlnPlot(inobj, "log_nfeature", group.by="idents")
  save_plot(p, paste(outfolder, "vplot_lognfeature_rna", sep="/"), fig.width=12, fig.height=4)

  p=Seurat::VlnPlot(inobj, "percent.mt", group.by="idents")
  save_plot(p, paste(outfolder, "vplot_percent_mt", sep="/"), fig.width=12, fig.height=4)

  p=Seurat::VlnPlot(inobj, "percent.rp", group.by="idents")
  save_plot(p, paste(outfolder, "vplot_percent_rp", sep="/"), fig.width=12, fig.height=4)

  p=Seurat::VlnPlot(inobj, "G2M.Score", group.by="idents")
  save_plot(p, paste(outfolder, "vplot_g2mscore", sep="/"), fig.width=12, fig.height=4)

  p=Seurat::VlnPlot(inobj, "S.Score", group.by="idents")
  save_plot(p, paste(outfolder, "vplot_sscore", sep="/"), fig.width=12, fig.height=4)


  if ((reduction == "umap") && ("wnn_clusters" %in% names(inobj@meta.data)))
  {
    p=Seurat::VlnPlot(inobj, "log_ncount", group.by="wnn_clusters")
    save_plot(p, paste(outfolder, "vplot_wnn_clusters_logncount_rna", sep="/"), fig.width=12, fig.height=4)

    p=Seurat::VlnPlot(inobj, "log_nfeature", group.by="wnn_clusters")
    save_plot(p, paste(outfolder, "vplot_wnn_clusters_lognfeature_rna", sep="/"), fig.width=12, fig.height=4)

    p=Seurat::VlnPlot(inobj, "percent.mt", group.by="wnn_clusters")
    save_plot(p, paste(outfolder, "vplot_wnn_clusters_percent_mt", sep="/"), fig.width=12, fig.height=4)
  }


  return(inobj)
}











#' Creates a summary of a vector of expression values
#'
#' @param a vector of expression values
#' @param suffix column suffix
#'
#' @return list of summarised expression data
#' @export
#'
makesummary_getPop = function(a, suffix)
{
  out = {}
  out["num"] = length(a)
  
  if (length(a) == 0)
  {
      f = c(0,0,0,0,0)
      meanA = 0
  } else {
      f = stats::fivenum(a)
      meanA = mean(a)
  }
  
  out["min"] = f[1]
  out["lower_hinge"] = f[2]
  out["median"] = f[3]
  out["upper_hinge"] = f[4]
  out["max"] = f[5]
  out["mean"] = meanA
  
  names(out) = paste(names(out), suffix, sep=".")
  
  return(out)
}



#' Calculates expression values for a set of cells
#'
#' @param markerObj Seurat object
#' @param markerCells cells for which gene expression data is to be calculated
#' @param sampleSuffix suffix of the column names
#' @param slot data slot to pull expression values from
#' @param assay the assay from which to pull the expression values
#'
#' @return list of summarised expression data
#' @export
#'
getExprData_getPop = function(markerObj, markerCells, sampleSuffix, slot="data", assay="RNA")
{
  expTable = Seurat::GetAssayData(object = subset(x=markerObj, cells=markerCells), slot = slot, assay=assay)
  allgenes = rownames(expTable)
  cellnames = colnames(expTable)
  
  expt.r = as(expTable, "TsparseMatrix") # was dgTMatrix
  expt.df = data.frame(r = expt.r@i + 1, c = expt.r@j + 1, x = expt.r@x)
  
  DT <- data.table::data.table(expt.df)
  res = DT[, as.list(makesummary_getPop(x, sampleSuffix)), by = r]
  anumCol = paste("anum", sampleSuffix, sep=".")
  res[[anumCol]] = length(cellnames)
  res$gene = allgenes[res$r]
  
  res = res[,r:=NULL]
  
  return(res)
}




#' Annotate a marker gene data frame with gene expression values
#'
#' @param scdata Seurat object
#' @param markers data frame of marker genes per cluster
#' @param assay assay to pull expression values from
#' @param group.by grouping of the clusters, same as the one used for markers
#'
#' @return data frame with gene expression values
#' @export
#'
getDEXpressionDF = function ( scdata, markers, assay="SCT", group.by=NULL)
{
  outDF = NULL
  Seurat::DefaultAssay(object=scdata) = assay
  print(group.by)
  if (is.null(group.by))
  {
    clusterIDs = as.character(sort(unique(Seurat::Idents(scdata))))
  } else {
    clusterIDs = as.character(sort(unique(scdata[[group.by]][,])))
  }
  scCells = Seurat::Idents(scdata)
  scCells = names(scCells)
  scCells = unlist(as.character(scCells))
  for (clusterID in clusterIDs){
      
    print(clusterID)
    
    cellIdents = Seurat::Idents(scdata)
    
    if (is.null(group.by))
    {
      cellIdents.c = names(cellIdents[cellIdents == clusterID])
    } else {
      cellIdents.c = colnames(scdata[,scdata[[group.by]] == clusterID])
    }
    
    cellIdents.c = unlist(lapply(cellIdents.c, as.character))  
    
    cellIdents.bg = setdiff(unlist(lapply(names(cellIdents), as.character)), cellIdents.c)

    if (length(cellIdents.c) < 3)
    {
      print(paste("Skipping cluster", clusterID, "due to < 3 cells!"))  
      next
    }
    
    expvals = getExprData_getPop(scdata, cellIdents.c, "cluster", assay=assay)
    expvals.bg = getExprData_getPop(scdata, cellIdents.bg, "bg", assay=assay)
    modmarkers = markers[[clusterID]]
    modmarkers$gene = rownames(modmarkers)
    
    markerdf = as.data.frame(modmarkers)
    
    if ((nrow(markerdf) > 0) && (nrow(expvals) > 0))
    {
      expvals = merge(markerdf, expvals, all.x=T, by.x="gene", by.y = "gene")  
    }
    
    if ((nrow(expvals) > 0) && (nrow(expvals.bg) > 0))
    {
      expvals = merge(expvals, expvals.bg, all.x=T, by.x="gene", by.y = "gene")  
    }
    
    expvals = as.data.frame(cbind(clusterID, expvals))
    
    if (!is.data.frame(outDF) || nrow(outDF)==0)
    {
    outDF = expvals
    } else {
    outDF = as.data.frame(rbind(outDF, expvals))
    }
      
  }
  return(outDF)
}



#' Performs Marker analysis for all clusters defined by the group.by grouping. Analysis is performed in the given assay using the given test.
#'
#' @param inobj Seurat object
#' @param group.by column in the meta.data to perform the analysis on
#' @param assay assay on which the analysis is performed on
#' @param test the test used to perform for finding marker genes
#'
#' @return data frame with marker genes for each group/cluster
#' @export
#'
makeDEResults = function(inobj, group.by=NULL, assay="SCT", test="wilcox")
{
  if (is.null(group.by))
  {
    clusterIDs = as.character(sort(unique(Seurat::Idents(inobj))))
  } else {
    clusterIDs = as.character(sort(unique(inobj[[group.by]][,])))
  }
  
  retList = list()
  
  for(clusterID in clusterIDs)
  {
  
      cellIdents = Seurat::Idents(inobj)
      
      if (is.null(group.by))
      {
        cellIdents.c = names(cellIdents[cellIdents == clusterID])
  
      } else {
        cellIdents.c = colnames(inobj[,inobj[[group.by]] == clusterID])
      }
      
      cellIdents.c = unlist(lapply(cellIdents.c, as.character))
  
      print(paste("Processing cluster", clusterID, "with a total of", length(cellIdents.c), "cells"))
      print(length(cellIdents.c) < 3)
  
      if (length(cellIdents.c) < 3)
      {
        print(paste("Skipping cluster", clusterID, "due to < 3 cells!"))  
        next
      }
      deMarkers = Seurat::FindMarkers(inobj, assay=assay, ident.1 = cellIdents.c, test.use=test) 
      retList[[clusterID]] = deMarkers
  
  }

  retList = getDEXpressionDF(inobj, retList, assay=assay, group.by=group.by)
  
  return(retList)
}


#' Performs a comparison between cellsID1 and cellsID2
#'
#'
#' @param inobj Seurat object
#' @param cellsID1 cell names of cells to include for group1
#' @param cellsID2 cell names of cells to include for group2
#' @param suffix1 name of group1
#' @param suffix2 name of group2
#' @param assay assay to use for comparison (should be RNA)
#' @param test test to perform for differential expression
#' @param outfolder output folder/stump where to store the results
#' @param fcCutoff minimal logFC to report
#' @param heatmap.plot whether to plot a heatmap for the comparison
#' @param heatmap.pcutoff p-value cut off for genes to show in the heatmap
#' @param heatmap.addgenes additional genes to include in heatmap
#'
#' @return differential expression data frame
#'
#' @export
compareClusters = function(scdata, cellsID1, cellsID2, suffix1, suffix2, prefix="cluster", test="MAST", assay="RNA", outfolder="./", fcCutoff=0.25, heatmap.plot=FALSE, heatmap.pcutoff=0.05, heatmap.addgenes=NULL)
{
    logfc.threshold = fcCutoff

    if (!dir.exists(outfolder)){
        dir.create(outfolder, recursive = TRUE)
        print(outfolder)
        print("Dir created!")

    } else {
        print(outfolder)
        print("Dir already exists!")
    }
    
    markers = Seurat::FindMarkers(scdata, assay=assay, ident.1 = cellsID1, ident.2 = cellsID2, test.use=test, logfc.threshold=logfc.threshold)
    
    outvalues1 = getExprData_getPop(scdata, cellsID1, suffix1, assay=assay)
    outvalues2 = getExprData_getPop(scdata, cellsID2, suffix2, assay=assay) 
    
    
    markers$gene = rownames(markers)
    joinedData = merge(markers, outvalues1, by="gene", all=T)
    joinedData = merge(joinedData, outvalues2, by="gene", all=T)  
    
    joinedData = joinedData[!is.na(joinedData$p_val),]
    
    suffix1=stringr::str_replace_all(str_replace_all(suffix1, "\\/", "_"), " ", "_")
    suffix2=stringr::str_replace_all(str_replace_all(suffix2, "\\/", "_"), " ", "_")

    
    outfile = paste(outfolder, "/", prefix, ".", suffix1, "_", suffix2, ".tsv", sep="")
    
    message(outfile)
    write.table(joinedData, file=outfile, row.names = F,  quote=FALSE, sep='\t')
    
    outfile = paste(outfolder, "/", prefix, ".", suffix1, "_", suffix2, ".xlsx", sep="")
    
    message(outfile)
    writexl::write_xlsx(joinedData, path=outfile)

    if (heatmap.plot)
    {
      genes.interest = (joinedData[joinedData$p_val_adj < heatmap.pcutoff,] %>% arrange(p_val_adj) %>% head(40) %>% top_n(40, avg_log2FC))$gene

      if (!is.null(heatmap.addgenes))
      {
        genes.interest = unique(c(genes.interest, heatmap.addgenes))
      }

      if (length(genes.interest) > 0)
      {
        scaleColors = c("#8b0000", "grey", "#008b2b")

        obj.rel = subset(scdata, cells=c(cellsID1, cellsID2))
        cellAnnot = colnames(obj.rel)
        names(cellAnnot) = cellAnnot

        cellAnnot[cellsID1] = suffix1
        cellAnnot[cellsID2] = suffix2
        obj.rel$heatmap_annot = cellAnnot

        obj.rel = Seurat::ScaleData(obj.rel, features=genes.interest)

        p=Seurat::DoHeatmap(obj.rel, genes.interest, group.by="heatmap_annot")+ ggplot2::scale_fill_gradientn(colors = scaleColors)
        p=p + ggplot2::ggtitle("Heatmap of DE genes; Data scaled by shown genes.")
        save_plot(p, paste(outfolder, paste("hplot_", prefix, ".", suffix1, "_", suffix2, sep=""), sep="/"), fig.width=7, fig.height=0.3*length(genes.interest))
      }
    }
    
      
    return(joinedData)
}




#' Splits all Seurat objects in the input list by a specific meta tag; useful for sample integration
#'
#'
#' @param inobj Seurat object
#' @param cellsID1 cell names of cells to include for group1
#' @param cellsID2 cell names of cells to include for group2
#' @param suffix1 name of group1
#' @param suffix2 name of group2
#' @param group.by meta-data column to perform the comparisons for
#' @param assay assay to use for comparison (should be RNA)
#' @param test test to perform for differential expression
#' @param outfolder output folder/stump where to store the results
#' @param fcCutoff minimal logFC to report
#' @param heatmap.plot whether to plot a heatmap for the comparison
#' @param heatmap.pcutoff p-value cut off for genes to show in the heatmap
#' @param heatmap.addgenes additional genes to include in heatmap
#'
#' @return list of differential expression data frame
#'
#' @export
compareCellsByCluster = function(inobj, cellsID1, cellsID2, suffix1, suffix2, group.by=NULL, assay="RNA", test="t", outfolder="./", fcCutoff=0.25, heatmap.plot=FALSE, heatmap.pcutoff=0.05, heatmap.addgenes=NULL)
{
  if (is.null(group.by))
  {
    clusterIDs = as.character(sort(unique(Seurat::Idents(inobj))))  
  } else {
    clusterIDs = as.character(sort(unique(inobj[[group.by]][,])))
  }
  
  retList = list()
  
  for(clusterID in clusterIDs)
  {
  
      cellIdents = Seurat::Idents(inobj)
      
      if (is.null(group.by))
      {
        cellIdents.c = names(cellIdents[cellIdents == clusterID])
  
      } else {
        cellIdents.c = colnames(inobj[,inobj[[group.by]] == clusterID])
      }
      
      cellIdents.c = unlist(lapply(cellIdents.c, as.character))
  
      print(paste("Processing cluster", clusterID, "with a total of", length(cellIdents.c), "cells"))
      
      cells1 = intersect(cellsID1, cellIdents.c)
      cells2 = intersect(cellsID2, cellIdents.c)
      
      print(c(length(cells1), length(cells2)))
      
      if (length(cells1) < 3)
      {
        print("Cells1 too few")
        next
      }
      
      if (length(cells2) < 3)
      {
        print("Cells2 too few")
        next
      }

      clusterID_file=stringr::str_replace_all(str_replace_all(clusterID, "\\/", "_"), " ", "_")

  
      deMarkers = compareClusters(scdata=inobj,
                                  cellsID1=cells1,
                                  cellsID2=cells2,
                                  prefix= paste("cluster", clusterID_file, sep="_"),
                                  suffix1=suffix1,#paste(suffix1, clusterID, sep="_"),
                                  suffix2=suffix2, #paste(suffix2, clusterID, sep="_"),
                                  test=test, fcCutoff=fcCutoff, assay=assay, outfolder=outfolder,
                                  heatmap.plot=heatmap.plot, heatmap.pcutoff=heatmap.pcutoff, heatmap.addgenes=heatmap.addgenes)
  
      retList[[clusterID]] = deMarkers
  
  }
  
  return(retList)
}












#' creates volcano plots for DE results list
#'
#' creates volcano plots for DE results list
#'
#' @param loMG list of dataframes from differential expression testing
#' @param titlePrefix Prefix of the title of the plots
#' @param outname path-stump for output (path + prefix of outname)
#' @param restrict_labels list of genes for each data frame to plot
#' @param turnExpression multiply logFC values by -1
#' @param colors list of list (neg(sig, nosig), pos(sig, nosig)) for colors
#' @param FCcutoff minimal logFC to include in plot
#' @param pCutoff minimal adjusted p-value to include in plot
#' @param highlightGene features to always include in the volcano plots
#'
#' @return None
#'
#'
#' @export
makeVolcanos = function(loMG, titlePrefix, outname, restrict_labels=NULL, turnExpression=F, colors = list(neg=list(sig="#448CCA", nosig="#B4D1E9"), pos=list(sig="#F47B78", nosig="#FACAC9")), FCcutoff=0.5, pCutoff = 0.05, highlightGene=NULL)
{

    outfolder = dirname(outname)[1]

    if (!dir.exists(outfolder)){
        dir.create(outfolder, recursive = TRUE)
        print(outfolder)
        print("Dir created!")

    } else {
        print(outfolder)
        print("Dir already exists!")
    }
  
  for (cName in names(loMG))
  {
    print(cName)
    
    title=paste(titlePrefix, "(", cName, ")", sep=" ")
    subtitle=""
    
    
    indf = loMG[[cName]]

    print(dim(indf))

    if (dim(indf)[1] == 0)
    {
      print("Skipping sample for 0 rows:")
      print(cName)
      print(dim(indf))
      next()
    }
    
    cName = stringr::str_replace_all(str_replace_all(cName, "\\/", "_"), " ", "_")
    popName = stringr::str_to_lower( stringr::str_replace_all(str_replace_all(str_replace_all( cName, "\\(|\\)| ", "_"), "__", "_"), "_$", "") )

    plotlabels = NULL
    
    if (!is.null(restrict_labels))
    {
      if (cName %in% names(restrict_labels))
      {
          cInfoElem = restrict_labels[[cName]]
        
          popName = cInfoElem$fname
          plotlabels = cInfoElem$genes
      } else {
        selDF = indf[indf$p_val_adj < pCutoff,] %>% arrange(p_val_adj) %>% head(40)
        plotlabels = selDF$gene
      }
    }    

    if (!is.null(highlightGene))
    {

      possibleGenes = intersect(indf$gene, highlightGene)
      print("possibleGenes")
      print(possibleGenes)

      if (length(possibleGenes) > 0)
      {
        subtitle = "including highlight gene"
        
        plotlabels = c(plotlabels, possibleGenes)
      }      
    }


    if (turnExpression)
    {
      indf$avg_log2FC = -indf$avg_log2FC
    }
    

    
    keyvals <- ifelse(indf$avg_log2FC > 0,
                      
                      ifelse(indf$p_val_adj < pCutoff & indf$avg_log2FC > FCcutoff, colors$pos$sig, colors$pos$nosig)
                      ,
                      ifelse(indf$avg_log2FC <= 0,
                             ifelse(indf$p_val_adj < pCutoff & indf$avg_log2FC < -FCcutoff, colors$neg$sig, colors$neg$nosig)
                             ,
                             'black'))
    


    keyvals[is.na(keyvals)] <- 'black'
    names(keyvals)[keyvals == '#F47B78'] <- 'UpReg sig'
    names(keyvals)[keyvals == '#FACAC9'] <- 'UpReg non-sig'
    names(keyvals)[keyvals == '#448CCA'] <- 'DownReg sig'
    names(keyvals)[keyvals == '#B4D1E9'] <- 'DownReg non-sig'
    
    txtLabelSize = 5
    axisLabelSize=12
    legendLabSize=12
    legendIconSize=5.0
  
    filename=paste(outname, popName, "png", sep=".")
    print(filename)
    png(filename=filename,width = 1200, height = 700)
    p=EnhancedVolcano::EnhancedVolcano(indf[,c('avg_log2FC', 'p_val_adj')],
      lab = indf$gene,
      selectLab=plotlabels,
      x = 'avg_log2FC',
      y = 'p_val_adj',
      colCustom = keyvals,
      pCutoff = pCutoff,
      FCcutoff = FCcutoff,
      pointSize = 2.0,
      legendPosition = 'right',
      drawConnectors = TRUE,
      widthConnectors = 0.5,
      title = title,
      subtitle = subtitle,
      axisLabSize=axisLabelSize,
      labSize = txtLabelSize,
      legendLabSize = legendLabSize,
      legendIconSize = legendIconSize,
      ylab = bquote(~-Log[10]~italic(adj~P-Value)),
      gridlines.major = FALSE,
      gridlines.minor = FALSE,
      max.overlaps=100
     )
    plot(p)
    dev.off()
    
    
        filename=paste(outname, popName, "pdf", sep=".")
    print(filename)
        pdf(filename,width = 12, height = 7)
    p=EnhancedVolcano::EnhancedVolcano(indf,
      lab = indf$gene,
      selectLab=plotlabels,
      x = 'avg_log2FC',
      y = 'p_val_adj',
      colCustom = keyvals,
      pCutoff = pCutoff,
      FCcutoff = FCcutoff,
      pointSize = 2.0,
      legendPosition = 'right',
      drawConnectors = TRUE,
      widthConnectors = 0.5,
      title = title,
      subtitle = subtitle,
      axisLabSize=axisLabelSize,
      labSize = txtLabelSize,
      legendLabSize = legendLabSize,
      legendIconSize = legendIconSize,
      ylab = bquote(~-Log[10]~italic(adj~P-Value)),
      gridlines.major = FALSE,
      gridlines.minor = FALSE,
      max.overlaps=100
     )
    plot(p)
    dev.off()
    
    filename=paste(outname, popName, "svg", sep=".")
    print(filename)
    svglite::svglite(file = filename, width = 12, height = 7)
    p=EnhancedVolcano::EnhancedVolcano(indf,
      lab = indf$gene,
      selectLab=plotlabels,
      x = 'avg_log2FC',
      y = 'p_val_adj',
      colCustom = keyvals,
      pCutoff = pCutoff,
      FCcutoff = FCcutoff,
      pointSize = 2.0,
      legendPosition = 'right',
      drawConnectors = TRUE,
      widthConnectors = 0.5,
      title = title,
      subtitle = subtitle,
      axisLabSize=axisLabelSize,
      labSize = txtLabelSize,
      legendLabSize = legendLabSize,
      legendIconSize = legendIconSize,
      ylab = bquote(~-Log[10]~italic(adj~P-Value)),
      gridlines.major = FALSE,
      gridlines.minor = FALSE,
      max.overlaps=100
     )
    plot(p)
    dev.off()
      
  }
  
}
