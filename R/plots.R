



#' Combines several ggplot2 objects
#'
#' Combines several ggplot2 objects with cowplot and also stores its data
#'
#' @param ... ggplot2 objects
#'
#' @return summarized ggplot2 object
#'
#' @export
combine_plot_grid = function(...)
{
  inplot = list(...)
  
  dataList = list()
  for (i in 1:length(inplot))
  {
    dataList[[length(dataList)+1]] = inplot[[i]]$data
  }
  
  p=cowplot::plot_grid(...)
  
  p$data = dataList
  
  return(p)
}


#' Combines several ggplot2 objects
#'
#' Combines several ggplot2 objects with cowplot and also stores its data
#'
#' @param plotlist ggplot2 objects
#' @param ... additional parameters forwarded to cowplot's plot_grid
#'
#' @return summarized ggplot2 object
#'
#' @export
combine_plot_grid_list = function(plotlist, ...)
{
 
  dataList = list()
  for (i in 1:length(plotlist))
  {
    dataList[[i]] = plotlist[[i]]$data
  }
  
  p=cowplot::plot_grid(plotlist=plotlist, ...)
  
  p$data = dataList
  
  return(p)
}


#' Saves a ggplot object to png, pdf and svg
#'
#' Also saves a data file containing all data relevant to reproduce the plot
#'
#' @param plotobj ggplot2 object to plot
#' @param outname path (without prefix) to where the plot should be stored, e.g. folder1/plot1
#' @param fig.width width of the plot
#' @param fig.height height of the plot
#' @param save.data whether to also save the data
#'
#' @return None
#'
#' @export
save_plot = function(plotobj, outname, fig.width, fig.height, save.data=TRUE)
{
  print(paste(outname, fig.width, fig.height))
  
  fname=paste(outname, "png", sep=".")
  print(paste("Saving to file", fname))
  png(filename=fname, width = fig.width, height = fig.height, units = 'in', res = 300)#width = fig.width*100, height=fig.height*100)
  plot(plotobj)
  dev.off()
  
  fname=paste(outname, "pdf", sep=".")
  print(paste("Saving to file", fname))
  pdf(file=fname, width = fig.width, height=fig.height)
  plot(plotobj)
  dev.off()
  

  fname=paste(outname, "svg", sep=".")
  print(paste("Saving to file", fname))
  svglite::svglite(file = fname, width = fig.width, height = fig.height)
  plot(plotobj)
  dev.off()
  

  if (save.data)
  {
    if (class(plotobj$data) %in% c("list"))
    {
      print("list case")
      for (i in 1:length(plotobj$data))
      {
        fname = paste(outname,i, "data", sep=".")
        print(paste("Saving to file", fname))
        
        if (class(plotobj$data[[i]]) %in% c("list"))
        {
            print("multi list case")
            for (j in 1:length(plotobj$data[[i]]))
            {
                fname = paste(outname,i, j, "data", sep=".")
                print(paste("Saving to file", fname, class(plotobj$data[[i]][[j]])))

                if (class(plotobj$data[[i]][[j]]) %in% c("list", "waiver"))
                {
                  next()
                }
                write.table(plotobj$data[[i]][[j]], fname, row.names = TRUE, sep="\t")    

            }
        } else {
            
            tryCatch(write.table(plotobj$data[[i]], fname, row.names = TRUE, sep="\t"), error = function(e) NULL)
        }
        
      }
    } else {
      
        fname = paste(outname,"data", sep=".")
        print(paste("Saving to file", fname))

        write.table(plotobj$data, paste(outname, "data", sep="."), row.names = TRUE, sep="\t")
    }
  }
  
  return(plotobj)
}



#' Creates a combined Violin- and Box-Plot
#'
#' Creates a combined Violin- and Box-Plot
#'
#' @param obj.sc Seurat object for plotting
#' @param gene genes to plot
#' @param group.by name of the meta.data column used for grouping cells
#' @param split.by name of the meta.data column used for dodging violins
#' @param pt.size size of the single value dots
#' @param assay assay to use for plotting
#' @param min.threshold minimum number of cells in group
#' @param col named list of colors
#' @param per.sample performs the min.threshold check per split-by-group
#' 
#'
#' @return ggplot2
#'
#'
#' @export
VlnBoxPlot = function( obj.sc, gene, group.by="idents", split.by=NULL, pt.size=0, assay="RNA", min.threshold=3, col=NULL, per.sample=FALSE)
{

    remainingCells = NULL

    if (is.null(split.by))
    {

        cellGroupDF = obj.sc[[ group.by ]]

        countDF = table(cellGroupDF)
        countDF = countDF[countDF >= min.threshold]

        remainingClusters = rownames(countDF)
        remainingCells = rownames(cellGroupDF)[ cellGroupDF[[group.by]] %in% remainingClusters ]

    } else {

      cellGroupDF = obj.sc[[ c(group.by, split.by) ]]

      countDF = as.data.frame(table(cellGroupDF))
      countDF = countDF[countDF$Freq >= min.threshold,]

      if (per.sample)
      {

        snames = unique(obj.sc[[ split.by ]][[split.by]])
        remainingCells = c()

        for (sname in snames)
        {
          snameDF = countDF[ countDF[split.by] == sname,]
          snameDF = snameDF[snameDF$Freq >= min.threshold,]

          sCells = rownames(obj.sc[[split.by]])[ obj.sc[[split.by]][[split.by]] == sname ]
          gCells = rownames(obj.sc[[group.by]])[ obj.sc[[group.by]][[group.by]] %in% snameDF[[group.by]]]

          remainingCells = c(remainingCells, intersect(sCells, gCells))
        }

      } else {
        remainingClusters = countDF[duplicated(countDF[group.by]),c(group.by)]

        cellGroupDF = obj.sc[[ group.by ]]
        remainingCells = rownames(cellGroupDF)[ cellGroupDF[[group.by]] %in% remainingClusters ]
      }



    }
    

    obj.sc.subset = subset(obj.sc, cells=remainingCells)
    print(obj.sc.subset)


    if (is.null(split.by))
    {

      ucols = NULL
      if (!is.null(col))
      {
        ncolors = length(unique(obj.sc[[group.by]][[group.by]]))
        ucols = rep(col, ncolors)
      }
      

      plots=Seurat::VlnPlot(obj.sc.subset, gene, group.by=group.by, split.by=split.by, pt.size=pt.size, assay=assay, col=ucols)
      plots = plots + ggplot2::geom_boxplot(color="grey", alpha=0.4) + ggplot2::stat_summary(fun=mean, geom="point", color="black", size=4)

    } else {

      plots=Seurat::VlnPlot(obj.sc.subset, gene, group.by=group.by, split.by=split.by, pt.size=pt.size, assay=assay, col=col)
      plots = plots + ggplot2::geom_boxplot(color="grey", alpha=0.4, position =position_dodge(width = 0.9)) + ggplot2::stat_summary(fun=mean, geom="point", aes(group=split), position=position_dodge(.9), color="black", size=4)

    }
    

    return(plots)
}

#' Creates a combined Violin- and Box-Plot with split violins
#'
#' Creates a combined Violin- and Box-Plot with split violins
#'
#' @param obj.sc Seurat object for plotting
#' @param gene genes to plot
#' @param group.by name of the meta.data column used for grouping cells
#' @param split.by name of the meta.data column used for dodging violins
#' @param pt.size size of the single value dots
#' @param assay assay to use for plotting
#' @param min.threshold minimum number of cells in group
#' @param col named list of colors
#' @param per.sample performs the min.threshold check per split-by-group
#'
#' @return ggplot2
#'
#' @export
SplitVlnBoxPlot = function( obj.sc, gene, group.by="idents", split.by=NULL, pt.size=0, assay="RNA", min.threshold=3, col=NULL, per.sample=TRUE)
{

    splitValues = unique(obj.sc[[split.by]][[split.by]])

    vplots = list()

    for (sname in splitValues)
    {
      print(sname)
      subsetcellnames = rownames(obj.sc[[split.by]])[ obj.sc[[split.by]][split.by] == sname ]

      if (is.null(col))
      {
        scol = NULL
      } else {
        scol = col[[sname]]
      }

      p=VlnBoxPlot( subset(obj.sc, cells=subsetcellnames), gene, group.by=group.by, split.by=NULL, pt.size=pt.size, assay=assay, min.threshold=min.threshold, col=scol, per.sample=per.sample)
      p = p + ggplot2::ggtitle(paste(gene, "-", sname))
      vplots[[sname]] = p
    }

    ps = combine_plot_grid_list(plotlist=vplots, nrow=1)

    return(ps)
}



#' Creates a combined Violin- and Box-Plot with split violins and significance testing
#'
#' Creates a combined Violin- and Box-Plot with split violins and significance testing
#'
#' @param obj.sc Seurat object for plotting
#' @param feature genes to plot
#' @param group.by name of the meta.data column used for grouping cells
#' @param split.by name of the meta.data column used for dodging violins
#' @param pt.size size of the single value dots
#' @param min.threshold minimum number of cells in group
#' @param dsrCols named list of colors
#' @param onelineLabel significance-test results in one line
#' @param min.threshold minimum number of cells in group
#' @param split.values values in split.by to consider
#' @param yStepIncrease y distance between statistics
#'
#' @return ggplot2 object 
#' 
#' @export
plotPValueViolentBoxPlot = function( obj.sc, feature, group.by, split.by=NULL, split.values=NULL, dsrCols=NULL, onelineLabel=FALSE, dot.size=0, min.threshold=3, yStepIncrease=0.25)
{
  # dsrCols = list("Thrombus"="grey", "Blood="red")
  


  remainingCells = NULL

  if (is.null(split.by))
  {

      cellGroupDF = obj.sc[[ group.by ]]

      countDF = table(cellGroupDF)
      countDF = countDF[countDF >= min.threshold]

      remainingClusters = rownames(countDF)
      remainingCells = rownames(cellGroupDF)[ cellGroupDF[[group.by]] %in% remainingClusters ]

  } else {

    cellGroupDF = obj.sc[[ c(group.by, split.by) ]]

    countDF = as.data.frame(table(cellGroupDF))

    groupValues = unique(cellGroupDF[, group.by])
    splitValues = unique(cellGroupDF[, split.by])
    remainingCells = c() 
    for (gv in groupValues)
    {
      
      keepGroup = TRUE
      for (sv in splitValues)
      {
        counts = (countDF %>% dplyr::filter(.data[[group.by]] == gv) %>% dplyr::filter(.data[[split.by]] == sv))$Freq

        if (counts < min.threshold)
        {
          keepGroup = FALSE
        }
      }

      if (keepGroup)
      {
        sCells = cellIDForClusters(obj.sc, split.by, splitValues)
        gCells = cellIDForClusters(obj.sc, group.by, gv)
        remainingCells = c(remainingCells, intersect(sCells, gCells))
      } else {
        print(paste("Removing group", gv))
      }
    }

  }


  print(paste("Existing Cells", length(colnames(obj.sc))))
  print(paste("New Cells", length(remainingCells)))

  obj.sc = subset(obj.sc, cells=remainingCells)

  if ((is.null(split.values) && !is.null(split.by)))
  {
    split.values = as.character(unique(obj.sc@meta.data[[split.by]]))
    print(split.values)
  }

  if (feature %in% colnames(obj.sc@meta.data))
  {
    dataDF = obj.sc@meta.data[,c(feature, group.by, split.by)]
  } else {
    dataDF = obj.sc@meta.data[,c(group.by, split.by)]
    dataDF[[feature]] = obj.sc@assays$RNA@data[feature, rownames(dataDF)]
  }
  
  if (!is.null(split.by))
  {
    dataDF = dataDF[dataDF[[split.by]] %in% split.values, ]
  }
  

  if (!is.null(levels(dataDF[[group.by]])))
  {
    dataDF[[group.by]]= factor(dataDF[[group.by]], levels=intersect( levels(dataDF[[group.by]]), as.character(dataDF[[group.by]]) ))
  } else {
    dataDF[[group.by]] = as.factor(dataDF[[group.by]])
  }

  if ((!is.null(split.by)) && (is.null(levels(dataDF[[split.by]]))))
  {
    dataDF[[split.by]] = as.factor(dataDF[[split.by]])
  }


  keepAC = c()
  for (ac in unique(dataDF[[group.by]]))
  {
    subdf = dataDF %>% filter(.data[[group.by]] == ac)

    if (!is.null(split.by))
    {
      numDiffGroups = length(unique(subdf[[split.by]]))
    } else {
      numDiffGroups = 2
    }
    

    print(paste(ac, numDiffGroups))

    if (numDiffGroups > 1)
    {
      keepAC = c(keepAC, ac)
    }
  }

  stat.test = dataDF %>% filter((!!as.symbol(group.by)) %in% keepAC) %>% group_by_at(group.by, .add=T)

  if (!is.null(split.by))
  {
   stat.test = stat.test %>%
  pairwise_t_test(
      as.formula(paste(feature, " ~ ", split.by, sep="")), paired = FALSE, 
      p.adjust.method = "BH"
      )
  } else {

        print(head(as.data.frame(stat.test)))


    stat.test = as.data.frame(stat.test) %>%
  pairwise_t_test(
      as.formula(paste(feature, " ~ ", group.by, sep="")), paired = FALSE, 
      p.adjust.method = "BH"
      )

    print(stat.test)
  }
  


  if (onelineLabel)
  {
    stat.test$label = paste("n1=",stat.test$n1," ", "n2=",stat.test$n2,", p < ", formatC(stat.test$p.adj, format = "e", digits = 2), sep="")
  } else {
    stat.test$label = paste("n1=",stat.test$n1,",\n", "n2=",stat.test$n2,"\n","p < ", formatC(stat.test$p.adj, format = "e", digits = 2), sep="")
  }

  maxValue = max(dataDF[,c(feature)])
  minValue = min(dataDF[,c(feature)])


  if (!is.null(split.by))
  {

    groupValues = unique(dataDF[, group.by])
    splitValues = unique(dataDF[, split.by])

    remainingCells = c() 
    for (gv in groupValues)
    {
      
      keepGroup = TRUE
      for (sv in splitValues)
      {
        subsetExprDF = (dataDF %>% dplyr::filter(.data[[group.by]] == gv) %>% dplyr::filter(.data[[split.by]] == sv))

        featureSD = sd(subsetExprDF[ , c(feature)])

        if (featureSD == 0)
        {
          dataDF[rownames(subsetExprDF)[1], c(feature)] = dataDF[rownames(subsetExprDF)[1], c(feature)] + 0.000001
        }
        
      }
    }
  }

  print("Creating Plot")

  if (is.null(split.by))
  {
    split.by=group.by
  }

  bxp <- ggplot(dataDF, aes_string(x=group.by, y=feature)) +
        geom_violin(aes_string(fill=split.by), position =position_dodge(width = 0.9), trim=TRUE)

  if (dot.size > 0)
  {
    bxp = bxp + geom_jitter(aes_string(fill=split.by), position =position_jitterdodge(dodge.width=0.9), size=dot.size)
  }

  bxp = bxp +
        geom_boxplot(aes_string(fill=split.by),color="grey", alpha=0.4, position =position_dodge(width = 0.9)) +
        theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1.0),
              panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"),panel.background = element_blank())+labs(fill='Group')
        
  if (!is.null(dsrCols))
  {
    bxp = bxp + scale_fill_manual(values = dsrCols[names(dsrCols) %in% split.values])
  }

  stat.test = stat.test[stat.test$p.adj < 0.05,]

  stat.test <- stat.test %>% rstatix::add_xy_position(x = group.by, dodge = 0.8, step.increase=yStepIncrease)#, y.trans=function(x){x*1.1})
  bxp = bxp + stat_pvalue_manual(stat.test,  label = "label", tip.length = 0)+ylim(c(floor(minValue*1.1), ceiling(maxValue*1.2)))
  return(bxp)
}



#' Creates an enhanced Dot Plot
#'
#' Creates an enhanced DotPlot of the RNA-assay expression
#'
#' @param obj.in Seurat object for plotting
#' @param plot_gois list of list(clusters,genes) where each entry defines the cluster to be shown, and genes the corresponding genes
#' @param group.by name of the meta.data column used for grouping cells
#' @param include_all_clusters whether to fill up the shown groups with the remaining ones from the category
#' @param title title of the plot
#' @param scale.by how to show/scale the expression values (NULL; ALL; GLOBAL)
#' @param scale.limits lower and upper bounds of the scale for gene expression
#' @param log.expression whether to log1p the expression values (only applies to GLOBAL and ALL scale)
#'
#' @return ggplot2
#'
#'
#' @export
enhancedHeatMap = function( obj.in, plot_gois, group.by="idents", include_all_clusters=FALSE, title=NULL, scale.by="GLOBAL", scale.limits=NULL, log.expression=TRUE)
{

combine_gois=TRUE
stopifnot(is.null(scale.by) || scale.by %in% c("GLOBAL", "ALL"))

if (!is.null(scale.by))
{
  cat("Scaling data", scale.by)
}


if (is.null(scale.by))
{
    print("Fetching average expression")
    avgexp = Seurat::AverageExpression(obj.in, assays=c("RNA"), group.by=group.by, slot="data")$RNA

} else {

  if (scale.by %in% c("ALL"))
  {
    print("Fetching average expression")
    avgexp = Seurat::AverageExpression(obj.in, assays=c("RNA"), group.by=group.by, slot="data")$RNA
  } else {
    print("Fetching global scaled average expression")
    avgexp = Seurat::AverageExpression(obj.in, assays=c("RNA"), group.by=group.by, slot="scale.data")$RNA
  }
}

if (combine_gois)
{
  new_clusters = c()
  new_genes = c()
  for (gname in names(plot_gois))
  {

    clusters = plot_gois[[gname]]$clusters
    genes = plot_gois[[gname]]$genes

    new_clusters = c(new_clusters, clusters)
    new_genes = c(new_genes, genes)
  }

  plot_gois = list("combined_lists"=list(clusters=new_clusters, genes=new_genes))
}

goi2cluster_genes = list()
allClusters = c()
allGenes = c()

for (gname in names(plot_gois))
{

    clusters = plot_gois[[gname]]$clusters
    genes = plot_gois[[gname]]$genes

    print(gname)
    print(clusters)
    print(genes)

    print("Clusters Missing")
    print(setdiff(clusters, colnames(avgexp)))

    print("Genes Missing")
    print(setdiff(genes, rownames(avgexp)))
   
    clusters = as.character(clusters)

    if (include_all_clusters)
    {
      missingClusters = setdiff( colnames(avgexp), clusters )
      print("Adding missing clusters")
      print(missingClusters)
      clusters = c(clusters, missingClusters)
    }

    allClusters = c(allClusters, clusters)
    allGenes = c(allGenes, genes)

    goi2cluster_genes[[gname]] = list("clusters"=clusters, "genes"=genes)
}
allClusters = unique(allClusters)
allGenes = unique(allGenes)


allScaleMatrices = list()
if (!is.null(scale.by) && (scale.by=="ALL"))
{

  if (length(allClusters) == 1)
  {
    print("1 cluster case")
    mat = avgexp[genes, c(allClusters[1], allClusters[1])]
    mat[1:length(genes), 1] = rep(NA, length(genes))

    origColnames = colnames(mat)
    origColnames[1] = "NA"
    colnames(mat) = origColnames

  } else {
    mat = avgexp[allGenes,allClusters]
  }

  if (log.expression)
  {
    mat = log1p(mat)
  }
  

  matsd = sd(as.vector(mat))
  matmean = mean(mat)
  mat = (mat-matmean) / matsd

  for (gname in names(goi2cluster_genes))
  {

      clusters = goi2cluster_genes[[gname]]$clusters
      genes = goi2cluster_genes[[gname]]$genes

      allScaleMatrices[[gname]] = mat[genes, clusters]

  }
  
}



  clusters = goi2cluster_genes[["combined_lists"]]$clusters
  genes = goi2cluster_genes[["combined_lists"]]$genes

  if (is.null(scale.by) || scale.by=="GLOBAL")
  {
    if (length(clusters) == 1)
    {
        print("1 cluster case")
        mat = avgexp[genes, c(clusters[1], clusters[1])]
        mat[1:length(genes), 1] = rep(NA, length(genes))
        print(mat)

        origColnames = colnames(mat)
        origColnames[1] = "NA"
        colnames(mat) = origColnames

    } else {
        mat = avgexp[genes,clusters]
    }

    if (is.null(scale.by))
    {
      if (log.expression)
      {
        mat = log1p(mat)
      }
    }
    
  
  } else {
    # ALL case!

    mat = allScaleMatrices[[gname]]
  }



  valueTitle = "Average Expression"

  if (!is.null(scale.by))
  {
      valueTitle = paste("Average Scaled Expression", " (",scale.by, ")", sep="")
  }

  if (is.null(title))
  {
    title = paste("Heatmap", gname)
  }

  if (!is.null(scale.limits))
  {

    mat[mat > scale.limits[2]] = scale.limits[2]
    mat[mat < scale.limits[1]] = scale.limits[1]

  }

  p = ComplexHeatmap::Heatmap(mat, name=valueTitle, rect_gp = gpar(col = "white", lwd = 2), column_title=title, cluster_rows = FALSE, row_order = rownames(mat), column_order = colnames(mat))

  return(p)

}


#' Creates an enhanced Dot Plot
#'
#' Creates an enhanced DotPlot of the RNA-assay expression
#'
#' @param obj.in Seurat object for plotting
#' @param plot_gois list of list(clusters,genes) where each entry defines the cluster to be shown, and genes the corresponding genes
#' @param group.by name of the meta.data column used for grouping cells
#' @param split.by name of the meta.data column used for splitting cells
#' @param include_all_clusters whether to fill up the shown groups with the remaining ones from the category
#' @param title title of the plot
#' @param scale.by how to show/scale the expression values (NULL; ALL; GLOBAL)
#' @param scale.limits lower and upper bounds of the scale for gene expression
#' @param log.expression whether to log1p the expression values (only applies to GLOBAL and ALL scale)
#'
#' @return ggplot2
#'
#'
#' @export
makeComplexExprHeatmapSplit = function( obj.in, plot_gois, split.by="condition", group.by="idents", include_all_clusters=FALSE, title=NULL, scale.by="GLOBAL", scale.limits=c(-2, 0, 2), log.expression=TRUE)
{

  combine_gois=TRUE
  stopifnot(is.null(scale.by) || scale.by %in% c("GLOBAL", "GROUP", "ALL"))

  if (!is.null(scale.by))
  {
    cat("Scaling data", scale.by)
  }

  if (combine_gois)
  {

    orig_gois = plot_gois

    new_clusters = c()
    new_genes = c()
    for (gname in names(plot_gois))
    {

      clusters = plot_gois[[gname]]$clusters
      genes = plot_gois[[gname]]$genes

      new_clusters = c(new_clusters, clusters)
      new_genes = c(new_genes, genes)
    }

    plot_gois = list()
    plot_gois[["combined_lists"]] = list(clusters=unique(new_clusters), genes=new_genes)

  }

  clusters = plot_gois[["combined_lists"]]$clusters
  genes = plot_gois[["combined_lists"]]$genes

  print(clusters)
  print(genes)

  valueTitle = "Average Expression"

  if (!is.null(scale.by))
  {
    valueTitle = paste("Average Scaled\nExpression (", scale.by, ")", sep="")
  }






  processedMats = list()
  splitByValues = unique(obj.in@meta.data[[split.by]])
  plot = NULL
  for (splitName in splitByValues)
  {
    cells.sel = cellIDForClusters(obj.in, split.by, c(splitName))
    #
    ## Fetching AVG Expression per subset
    #
    if (is.null(scale.by) || scale.by %in% c("ALL", "GROUP"))
    {
      print("Fetching average expression")
      avgexp = Seurat::AverageExpression(subset(obj.in, cells=cells.sel), assays=c("RNA"), group.by=group.by, slot="data")$RNA

    } else {

      print("Fetching global scaled average expression")
      avgexp = Seurat::AverageExpression(subset(obj.in, cells=cells.sel), assays=c("RNA"), group.by=group.by, slot="scale.data")$RNA
    }  

    #
    ## Checking all clusters and genes there!
    #
    print("Clusters Missing")
    print(setdiff(clusters, colnames(avgexp)))

    print("Genes Missing")
    print(setdiff(genes, rownames(avgexp)))

    p.clusters = as.character(clusters)

    if (include_all_clusters)
    {
      missingClusters = setdiff( colnames(avgexp), p.clusters )
      print("Adding missing clusters")
      print(missingClusters)

      if (length(missingClusters) > 0)
      {
        p.clusters = c(p.clusters, missingClusters)
      }
    }

    #
    ## Subsetting mat as required
    #
    if (length(p.clusters) == 1)
    {
        print("1 cluster case")
        mat = avgexp[genes, c(p.clusters[1], p.clusters[1])]
        mat[1:length(genes), 1] = rep(NA, length(genes))
        print(mat)

        origColnames = colnames(mat)
        origColnames[1] = "NA"
        colnames(mat) = origColnames

    } else {

      print("Removing clusters because they're not represented in subset")
      print(setdiff(p.clusters, colnames(avgexp)))

      p.clusters = intersect(p.clusters, colnames(avgexp))
      mat = avgexp[genes,p.clusters]
    }

    if (is.null(scale.by) || scale.by%in%c("GROUP", "ALL"))
    {
      if (log.expression)
      {
        mat = log1p(mat)
      }
    }

    if (!is.null(scale.by) && scale.by == "GROUP")
    {
      matsd = sd(as.vector(mat))
      matmean = mean(mat)
      mat = (mat-matmean) / matsd
    }

    processedMats[[splitName]] = mat
  }

  if (!is.null(scale.by) && scale.by=="ALL")
  {


    valueVec = c()
    for (splitName in names(processedMats))
    {
      valueVec = c(valueVec, as.vector(processedMats[[splitName]]))
    }

    allsd = sd(valueVec)
    allmean = mean(valueVec)

    for (splitName in names(processedMats))
    {
      processedMats[[splitName]] = (processedMats[[splitName]]-allmean)/allsd

      print(processedMats[[splitName]])
    }

  }
  
  
  for (splitName in names(processedMats))
  {

    showlegend = splitName == splitByValues[1]

    if (is.null(title))
    {
      plottitle = paste("Heatmap", gname, splitName)
    } else {
      plottitle = paste(title, " (", splitName, ")", sep="")
    }

    mat = processedMats[[splitName]]

    if (!is.null(scale.limits))
    {

      mat[mat > scale.limits[3]] = scale.limits[3]
      mat[mat < scale.limits[1]] = scale.limits[1]

    }
  
    col_fun = circlize::colorRamp2(scale.limits, c("blue", "white", "red"))
    print(col_fun)

    p = ComplexHeatmap::Heatmap(mat, col=col_fun,show_heatmap_legend=showlegend, name=valueTitle, rect_gp = gpar(col = "white", lwd = 2), column_title=plottitle, cluster_rows = FALSE, row_order = rownames(mat), column_order = colnames(mat))
    
    if (showlegend)
    {
        plot = p
    } else {
        plot = plot + p
    }

  }


  return(plot)


}





















#' Creates an enhanced DotPlot
#'
#' @param scobj Seurat object for plotting
#' @param plotElems list of list(cells,label) where each entry defines one condition/split of the plot
#' @param featureGenes genes to plot
#' @param group.by name of the meta.data column used for grouping cells
#' @param scale.by how to show/scale the expression values (GROUP; FEATURE; ALL; GLOBAL)
#' @param col.min lower bound of the scale
#' @param col.max upper bound of the scale
#' @param cols color for the expression values
#' @param title title of the plot
#' @param rotate.x whether to rotate x-axis labels
#' @param abundance.perelem whether the cell group abundance is to be calculated on the global Seurat object, or per group
#' @param assay which assay to use for retrieving the expression values
#'
#' @return ggplot2 object
#'
#'
#' @export
enhancedDotPlot = function(scobj, plotElems, featureGenes = c(""), group.by="cellnames_manual", col.min = -3, col.max = 3, cols = c("grey", "blue"), title="", scale.by="GROUP", rotate.x=F, abundance.perelem=FALSE, assay="RNA")
{

  stopifnot(is.null(scale.by) ||scale.by %in% c("GROUP", "FEATURE", "ALL", "GLOBAL"))

  use.slot="data"

  if (!is.null(scale.by) && scale.by == "GLOBAL")
  {
    use.slot="scale.data"
  }

  plotData = list()
  allFoundFeatures = c() 
  ctFractions = list()
  
  scFullTable = table(scobj[[group.by]])
  scFullDf = as.data.frame(scFullTable)

  fillLimitsMax = 0
  elemCount = 0
  for (plotName in names(plotElems))
  {
    plotDef = plotElems[[plotName]]
    plotCells = plotDef$cells
    plotDescr = plotDef$label
    
    scobj_subset = subset(scobj, cells=plotCells)

    avgexpMat = Seurat::AverageExpression(scobj_subset, features=featureGenes, group.by=group.by, assay=assay, slot=use.slot)$RNA
    
    if (use.slot=="data")
    {
      avgexpMat = log1p(avgexpMat)
    }
    

    avgExpr = as.data.frame(data.table::data.table(features.plot = rownames(avgexpMat), id = colnames(avgexpMat), avg.exp = c(as.matrix(avgexpMat))))

    p=Seurat::DotPlot(scobj_subset, features=featureGenes, group.by=group.by, assay=assay)

    avgExpr = merge(avgExpr, p$data[, c("id", "features.plot", "pct.exp")], by=c("id", "features.plot"))
    
    scTable = table(scobj_subset[[group.by]])
    scDf = as.data.frame(scTable)
    scobj_subset = NULL
    
    if (abundance.perelem)
    {
      scDf$perc = scDf$Freq / sum(scDf$Freq) # per elem abundance!
    } else {
      scDf$perc = scDf$Freq / sum(scFullDf$Freq) # total abundance!
    }

    fillLimitsMax = max(c(fillLimitsMax,max(scDf$perc)))

    ctFractions[[plotName]] = scDf   
    plotData[[plotName]] = avgExpr
  }

  idLevels = levels(scobj@meta.data[, c(group.by)])
  featureLevels = featureGenes

  
  # initialize Values
  for (plotName in names(plotData))
  {
    #plotData[[plotName]]$avg.exp.scaled2 = 0
    plotData[[plotName]]$id = as.character(plotData[[plotName]]$id)
    plotData[[plotName]]$features.plot = as.character(plotData[[plotName]]$features.plot)
  }

  #print(plotData[["Young"]])
  
  allIDs = c()
  allFeatures = c()
  for (plotName in names(plotData))
  {
      allIDs = c(allIDs, plotData[[plotName]]$id)
      allFeatures = c(allFeatures, plotData[[plotName]]$features.plot)
  }
  allIDs = unique(allIDs)
  allFeatures = unique(allFeatures)

  if ((is.null(scale.by) || scale.by == "GLOBAL"))
  {

    for (plotName in names(plotData))
    {

      plotDF = plotData[[plotName]]

      #plotDF$id = factor(plotDF$id, levels = mixedsort(as.character(unique(plotDF$id))))
      #plotDF$idn= as.numeric(plotDF$id)
      #plotDF$features.plot = factor(plotDF$features.plot, levels=unique(plotDF$features.plot))
      
      #just copy old values
      plotDF[, "avg.exp.scaled2"] = plotDF$avg.exp

      plotData[[plotName]] = plotDF
    }

  } else if (scale.by == "GROUP")
  {

    for (plotName in names(plotData))
    {

      plotDF = plotData[[plotName]]

      #just copy old values
      plotDF[, "avg.exp.scaled2"] = scale(plotDF$avg.exp)

      plotData[[plotName]] = plotDF
    }


  } else if (scale.by == "FEATURE")
  {
    print("SCALING BY FEATURE")
    # calculate avg.exp.scaled2 for each feature
    for (featureName in allFoundFeatures)
    {
      print(featureName)
      allUnscaledValues = NULL
      for (plotName in names(plotData))
      {
        pData = plotData[[plotName]]
        missingIDs = setdiff(allIDs, unique(pData$id))
                
        for (mid in missingIDs)
        {
          for (feature in allFeatures)
          {
            dfrow = data.frame(avg.exp=0.0,pct.exp=0.0,features.plot=feature, id=mid, avg.exp.scaled=0.0, avg.exp.scaled2=0.0)
            pData = rbind.data.frame(pData, dfrow, stringsAsFactors = F)
          }
        }
        
        #pData$id = factor(pData$id, levels = allIDs)
        #pData$features.plot = factor(pData$features.plot, levels=allFeatures)
        
        plotData[[plotName]] = pData
        
        
        if (is.null(allUnscaledValues))
        {
          allUnscaledValues = data.frame(plotName=pData[ pData$features.plot==featureName, ]$avg.exp)
          allUnscaledValues[[plotName]] = pData[ pData$features.plot==featureName, ]$avg.exp  
          allUnscaledValues[["plotName"]] = NULL
          
        } else {
          allUnscaledValues[[plotName]] = pData[ pData$features.plot==featureName, ]$avg.exp  
        }
        
      }

      allUnscaledValues$rnames = as.numeric(rownames(allUnscaledValues))
      allUnscaledLong = allUnscaledValues %>% gather(Type, Value, names(plotData))
      allUnscaledLong$Value = scale(allUnscaledLong$Value)
      allScaledValues = allUnscaledLong %>% spread(Type, Value) %>% arrange( order(rnames))
      
      for (plotName in names(plotData))
      {

        pData = plotData[[plotName]]
               
        origData = pData[pData$features.plot==featureName, ]
        pData[pData$features.plot==featureName, "avg.exp.scaled2"] = pData[pData$features.plot==featureName, "avg.exp"]
        #pData$idn=as.numeric(pData$id)

        plotData[[plotName]] = pData
      }


      
    }

  } else if (scale.by == "ALL") {

    combinedDataDF = data.frame()

    for (plotName in names(plotData))
    {
      pData = plotData[[plotName]]
      pData["plotpart"] = plotName

      combinedDataDF = rbind(combinedDataDF, pData)
    }

    combinedDataDF[, c("avg.exp.scaled2")] = scale(combinedDataDF$avg.exp)


    # reorder
    originalGroups = scobj@meta.data[,group.by]
    if (is.factor(originalGroups))
    {
      originalSorting = levels(originalGroups)
      combinedDataDF$id = factor(combinedDataDF$id, levels = originalSorting)
    } else {
      combinedDataDF$id = factor(combinedDataDF$id, levels = gtools::mixedsort(as.character(unique(combinedDataDF$id))))
    }


    for (plotName in names(plotData))
    {
      subDF = combinedDataDF[combinedDataDF$plotpart == plotName, c("features.plot", "id", "avg.exp.scaled2", "plotpart")]

      pData = plotData[[plotName]]

      pData = merge(pData, subDF, by.x=c("features.plot", "id"), by.y=c("features.plot", "id"))

      #pData$avg.exp.scaled2 = subDF[pData$features.plot,]$avg.exp.scaled2
      #pData$id = plotData$id
      #pData$idn= as.numeric(pData$id)
      #pData$features.plot = factor(pData$features.plot, levels=unique(pData$features.plot))
      
      plotData[[plotName]] = pData

    }
      
  } else {
    stopifnot(FALSE)
  }


  for (plotName in names(plotData))
  {
    pData = plotData[[plotName]]

    pData$id = factor(pData$id, levels=idLevels)
    pData$idn = as.numeric(pData$id)
    pData$features.plot = factor(pData$features.plot, levels=featureLevels)

    plotData[[plotName]] = pData
  }


  # prepare ID DF
  idDF = data.frame()
  for (plotName in names(plotData))
  {
    pData = plotData[[plotName]]
    idDF = rbind(idDF, pData[, c("id", "idn")])
  }

  rownames(idDF) = NULL
  idDF = unique(idDF)

  if (!is.null(scale.by))
  {
    title.expression = paste("Avg. Expression\n(scaled by ", scale.by, ")", sep="")
  } else {
    title.expression = "Avg. Expression"
  }

  if (abundance.perelem)
  {
    title.cellabundance = "Cell Abundance\n(per condition)"
  } else {
    title.cellabundance = "Cell Abundance\n(all object cells)"
  }

  plotList = list()
  
  for (plotName in names(plotData))
  {
    pData = plotData[[plotName]]

    print(pData)

    pData2 = merge(x=pData,y=ctFractions[[plotName]],by.x="id", by.y=group.by,all.x=TRUE)
    pData <-pData2 %>% dplyr::mutate(featuren=as.numeric(features.plot), percn=100*perc)  

    print(pData)
    
    pData[pData$avg.exp.scaled2>col.max, 'avg.exp.scaled2'] = col.max
    pData[pData$avg.exp.scaled2<col.min, 'avg.exp.scaled2'] = col.min

    fillLimits = c(0, ceiling(fillLimitsMax*10)*10)

    print(pData)
    
    plotElem <- ggplot2::ggplot(pData) +
      ggplot2::scale_x_continuous(breaks=pData$featuren, labels=pData$features.plot) +
      ggplot2::scale_y_continuous(breaks=idDF$idn, labels=idDF$id, limits = c(min(idDF$idn)-0.6, max(idDF$idn)+0.6))+
      ggplot2::geom_rect(ggplot2::aes(xmin=featuren-.5, xmax=featuren+.5, ymin = idn-0.5, ymax = idn+0.5, fill=percn), alpha = 0.4, linetype="blank") +
      ggplot2::scale_fill_distiller(palette='Spectral', limits = fillLimits)+
      ggplot2::scale_size_continuous(range = c(0, 10))+
      ggplot2::geom_point(ggplot2::aes(x=featuren, y=idn, colour = avg.exp.scaled2, size = pct.exp)) +
      ggplot2::scale_color_gradient2(limits=c(col.min, col.max), low = cols[1], mid=cols[2], high = cols[3])+
      ggplot2::guides(color=ggplot2::guide_colourbar(title=title.expression, order = 1),
      size=ggplot2::guide_legend(title="Percent Expressing", order = 2),
      fill=ggplot2::guide_colourbar(title=title.cellabundance, order = 3))
    
    plotElem = plotElem + ggplot2::theme(axis.title.x=ggplot2::element_blank(), axis.title.y=ggplot2::element_blank(), legend.position = "none", panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line())

    if (rotate.x)
    {
      plotElem = plotElem + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
    }

    plotList[[plotName]] = plotElem
  }
    
  print("Final plot checks")

  plot_list=list()
  for (plotName in names(plotList))
  {
    
    plotElem = plotList[[plotName]]      
    
    print("descr label")
    pe = make_descr_label(plotElem, plotName) #plotElems[[plotName]]$label
  
    plot_list[[plotName]] = pe
 
  }

  print("Preparing Legend")
  
  # extract a legend that is laid out horizontally #
  legend_b <- ggpubr::get_legend(
    plotElem + 
      ggplot2::guides(color = ggplot2::guide_colorbar(title = title.expression, direction="horizontal"))+ #guide_legend(nrow = 1, override.aes = list(size=6)))+
      ggplot2::theme(legend.position = "bottom", legend.box = "vertical")
  )
  
  title <- cowplot::ggdraw() + cowplot::draw_label(title, fontface='bold')
  print("Preparing plot")

  dataList = list()
  for (i in 1:length(plot_list))
  {
    dataList[[i]] = plot_list[[i]]$data
  }
  ap=combine_plot_grid_list(
    plotlist = plot_list,
    labels = NULL,
    nrow=1,
    align = "v", axis="bt"
  )

  finalPlotList = list(title, ap, legend_b)

  fp = combine_plot_grid_list(plotlist = finalPlotList, ncol = 1, rel_heights = c(.05, 1, 0.2) )
  fp$data = dataList

  return(fp)
}



#
##
###
#### Plotting tools
###
##
#


cbind.fill <- function(...){
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

make_descr_label = function(plot, descr)
{
  descrLabel <- cowplot::ggdraw() + cowplot::draw_label(descr, fontface='bold', angle = 0)
  
  pe = combine_plot_grid_list(plotlist=list("a"=descrLabel, "b"=plot), ncol=1, nrow=2, labels=NULL,rel_heights = c(0.1, 1), align = "h", axis = "l")
  
  return(pe)
}


