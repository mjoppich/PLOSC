





#' performEnrichtmentAnalysis - performs enrichment analysis on differential expression results provided by compareCellsByCluster
#'
#' @param scobj Seurat object with active assay from which the DE results are (used to get gene universe for set enrichment)
#' @param deList differential expression results provided by compareCellsByCluster
#' @param organismName human or mouse
#' @param outFile where to store the enrichment results as Rds file
#' @param reverseLogFC set true if you want to reverse logFC values (logFC=-logFC)
#'
#' @return list with enrichment results for each comparison given in deList
#' @export
#'
#' @examples
performEnrichtmentAnalysis = function(scobj, deList, organismName, outFile, reverseLogFC=FALSE) {
  
  if (organismName == "human")
  {
    print("Human Org DB")
    orgShort = "hsa"
    #BiocManager::install("org.Hs.eg.db")
    globOrgDB = org.Hs.eg.db::org.Hs.eg.db
  } else {
    orgShort = "mmu"
    #BiocManager::install("org.Mm.eg.db")
    globOrgDB = org.Mm.eg.db::org.Mm.eg.db
  }
  
  
  allDEResultIDs = names(deList)
  gseResults = list()
  
  universeSym = rownames(scobj)
  
  geneNames = clusterProfiler::bitr(universeSym, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = globOrgDB)
  universeEntrez = geneNames[universeSym %in% geneNames$SYMBOL,]$ENTREZID
  universeEntrez = universeEntrez[!is.na(universeEntrez)]
  
  for (resID in allDEResultIDs)
  {
    print(resID)
    resDF = deList[[resID]]
    
    
    #
    #
    #
    #
    
    print("Selecting Genes")
    
    if ((is.null(resDF)) || (nrow(resDF) == 0))
    {
      print("NULL genes")
      gseResults[[resID]] = list(valid=F, rao=NULL, rao_up=NULL, rao_down=NULL)
      next()
    }
    
    # select significantly regulated genes
    siggenes = resDF[resDF$p_val_adj < 0.05,]
    
    print(head(siggenes))  
    if (! "gene" %in% colnames(siggenes))
    {
      siggenes$gene = rownames(siggenes)
    }
    
    # reverse avg logFC
    if (reverseLogFC)
    {
      siggenes$avg_log2FC = -siggenes$avg_log2FC
    }
    
    #geneRegVec = siggenes$avg_logFC
    #names(geneRegVec) = siggenes$gene
    
    
    # convert gene symbols to entrez ID for gseGO
    geneNames = clusterProfiler::bitr(siggenes$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = globOrgDB)
    print("Unresolved gene symbols")
    print(siggenes[!(siggenes$gene %in% geneNames$SYMBOL),]$gene)
    
    siggenes = siggenes[siggenes$gene %in% geneNames$SYMBOL,]
    
    print(head(siggenes))
    
    geneVector = as.vector(siggenes$avg_log2FC)
    names(geneVector) = siggenes$gene
    
    geneNameMap = geneNames$ENTREZID
    names(geneNameMap) = geneNames$SYMBOL
    names(geneVector) = as.numeric(geneNameMap[names(geneVector)])
    
    # sort decreasing ...
    # we sort by decreasing logFC: geneList = sort(geneList, decreasing = TRUE) https://github.com/YuLab-SMU/DOSE/wiki/how-to-prepare-your-own-geneList
    geneVector = sort(geneVector, decreasing = TRUE)
    geneNameMap = sort(geneNameMap, decreasing = TRUE)
    
    keepElems = !is.na(geneVector)
    geneVector <- geneVector[keepElems]
    geneNameMap <- geneNameMap[keepElems]
    
    print("Got geneRegVec")
    print(length(geneVector))
    
    
    if (length(geneVector) == 0)
    {
      gseResults[[resID]] = list(valid=F, rao=NULL, rao_up=NULL, rao_down=NULL)
      next()
    }
    
    
    upGenes=geneVector[geneVector > 0]
    upGenesSym = geneNameMap[geneVector > 0]
    print(paste("UpGenes", length(upGenes)))
    
    downGenes = geneVector[geneVector < 0]
    downGenesSym = geneNameMap[geneVector < 0]
    print(paste("downGenes", length(downGenes)))
    
    
    #
    ## REACTOME
    #
    rao = ReactomePA::enrichPathway(gene=names(geneVector), organism=organismName, universe=universeEntrez, pvalueCutoff = 1, qvalueCutoff=1, readable=TRUE, pAdjustMethod="BH") 
    rao_up = ReactomePA::enrichPathway(gene=names(upGenes), organism=organismName, universe=universeEntrez, pvalueCutoff = 1, qvalueCutoff=1, readable=TRUE, pAdjustMethod="BH")
    rao_down = ReactomePA::enrichPathway(gene=names(downGenes), organism=organismName, universe=universeEntrez, pvalueCutoff = 1, qvalueCutoff=1, readable=TRUE, pAdjustMethod="BH")
    
    #
    ## KEGG
    #
    keggo = clusterProfiler::enrichKEGG(gene=names(geneVector), organism=orgShort, universe=universeEntrez, pvalueCutoff = 1, qvalueCutoff=1, pAdjustMethod="BH")
    if (!is.null(keggo)) keggo <- clusterProfiler::setReadable(keggo, OrgDb = globOrgDB, keyType="ENTREZID")
    
    keggo_up = clusterProfiler::enrichKEGG(gene=names(upGenes), organism=orgShort, universe=universeEntrez, pvalueCutoff = 1, qvalueCutoff=1, pAdjustMethod="BH")
    if (!is.null(keggo_up)) keggo_up <- clusterProfiler::setReadable(keggo_up, OrgDb = globOrgDB, keyType="ENTREZID")
    
    keggo_down = clusterProfiler::enrichKEGG(gene=names(downGenes), organism=orgShort, universe=universeEntrez, pvalueCutoff = 1, qvalueCutoff=1, pAdjustMethod="BH")
    if (!is.null(keggo_down)) keggo_down <- clusterProfiler::setReadable(keggo_down, OrgDb = globOrgDB, keyType="ENTREZID")
    
    #
    ## GO
    #
    goo = clusterProfiler::enrichGO(gene=names(geneVector), OrgDb=globOrgDB, universe=universeEntrez, pvalueCutoff = 1, qvalueCutoff=1, pAdjustMethod="BH", readable=TRUE)
    goo_up = clusterProfiler::enrichGO(gene=names(upGenes), OrgDb=globOrgDB, universe=universeEntrez, pvalueCutoff = 1, qvalueCutoff=1, pAdjustMethod="BH", readable=TRUE)
    goo_down = clusterProfiler::enrichGO(gene=names(downGenes), OrgDb=globOrgDB, universe=universeEntrez, pvalueCutoff = 1, qvalueCutoff=1, pAdjustMethod="BH", readable=TRUE)
    
    #
    ## MKEGG
    #
    mkeggo = clusterProfiler::enrichMKEGG(gene=names(geneVector), organism=orgShort, universe=universeEntrez, pvalueCutoff = 1, qvalueCutoff=1, pAdjustMethod="BH")
    if (!is.null(mkeggo)) mkeggo <- clusterProfiler::setReadable(mkeggo, OrgDb = globOrgDB, keyType="ENTREZID")
    
    mkeggo_up = clusterProfiler::enrichMKEGG(gene=names(upGenes), organism=orgShort, universe=universeEntrez, pvalueCutoff = 1, qvalueCutoff=1, pAdjustMethod="BH")
    if (!is.null(mkeggo_up)) mkeggo_up <- clusterProfiler::setReadable(mkeggo_up, OrgDb = globOrgDB, keyType="ENTREZID")
    
    mkeggo_down = clusterProfiler::enrichMKEGG(gene=names(downGenes), organism=orgShort, universe=universeEntrez, pvalueCutoff = 1, qvalueCutoff=1, pAdjustMethod="BH")
    if (!is.null(mkeggo_down)) mkeggo_down <- clusterProfiler::setReadable(mkeggo_down, OrgDb = globOrgDB, keyType="ENTREZID")
    
    
    gseRes = list(valid=TRUE,
                  rao=rao, rao_up=rao_up, rao_down=rao_down,
                  keggo=keggo, keggo_up=keggo_up, keggo_down=keggo_down,
                  mkeggo=mkeggo, mkeggo_up=mkeggo_up, mkeggo_down=mkeggo_down,
                  goo=goo, goo_up=goo_up, goo_down=goo_down,
                  expression=list(all_genes_entrez=geneVector, up_genes_entrez=upGenes, down_genes_entrez=downGenes,
                                  all_genes_symbol=geneNameMap, up_genes_symbol=upGenesSym, down_genes_symbol=downGenesSym)
    )
    
    
    #
    #
    #
    #
    
    gseResults[[resID]] = gseRes
    
  }
  
  if (!is.null(outFile))
  {
    saveRDS(gseResults, outFile) 
  }
  
  return(gseResults)
  
}

#' createDirIfNotThere - creates directory if it does not exist yet
#'
#' @param which directory to create
#'
#' @return
#' @export
#'
#' @examples
createDirIfNotThere = function(infolder)
{
  if (!dir.exists(infolder))
  {
    dir.create(infolder, recursive=TRUE)
  }
}






#' makeEnrichmentPlots - makes plots for the enrichment results generated by performEnrichtmentAnalysis
#'
#' @param enrichResult list of enrichment results as generated by performEnrichtmentAnalysis
#' @param outfolder folder where to store plots from enrichment analysis
#' @param qvalue_threshold qvalue threshold up to where to accept gene sets
#'
#' @return
#' @export
#'
#' @examples
makeEnrichmentPlots = function(enrichResult, outfolder, qvalue_threshold=0.2)
{
  
  createDirIfNotThere(paste(outfolder, "enrichment", sep="/"))
  
  
  for (comparison in names(enrichResult))
  {
    
    createDirIfNotThere(paste(outfolder, "enrichment", comparison, sep="/"))
    
    
    all_tests = names(enrichResult[[comparison]])
    all_tests = setdiff(all_tests, c("valid", "expression"))
    
    geneExpression = enrichResult[[comparison]][["expression"]]
    
    geneExpression_vec = geneExpression[["all_genes_entrez"]]
    geneExpression_vec_up = geneExpression[["up_genes_entrez"]]
    geneExpression_vec_down = geneExpression[["down_genes_entrez"]]
    
    
    for (perf_test in all_tests)
    {
      
      print(paste(comparison, perf_test))
      result_orig = enrichResult[[comparison]][[perf_test]]
      
      if (is.null(result_orig))
      {
        next()
      }
      
      result = clusterProfiler::filter(result_orig, Count > 0)
      result = clusterProfiler::filter(result, qvalue < qvalue_threshold)
      
      
      if (nrow(result) == 0)
      {
        next()
      }
      
      createDirIfNotThere(paste(outfolder,"enrichment", comparison, perf_test, sep="/"))
      
      
      write.table(result_orig@result, paste(outfolder,"enrichment", comparison, perf_test, "enrichment_result.tsv", sep="/"), sep="\t", row.names=F, quote = F)
      writexl::write_xlsx(result_orig@result, paste(outfolder,"enrichment", comparison, perf_test, "enrichment_result.xlsx", sep="/"))
      
      p= enrichplot::dotplot(result, showCategory=30) + ggplot2::ggtitle(paste("dotplot for", perf_test))
      save_plot(p, paste(outfolder,"enrichment", comparison, perf_test, "dotplot_30", sep="/"), fig.width=6, fig.height=12)
      
      p= barplot(dplyr::mutate(result, qscore = -log(qvalue, base=10)), x="qscore", fill="qvalue", showCategory=30)
      save_plot(p, paste(outfolder,"enrichment", comparison, perf_test, "barplot_qvalue_30", sep="/"), fig.width=6, fig.height=12)
      
      testExpression = NULL
      if (endsWith(perf_test, "_up"))
      {
        testExpression = geneExpression_vec_up
      } else if (endsWith(perf_test, "_down"))
      {
        testExpression = geneExpression_vec_down
      } else {
        testExpression = geneExpression_vec
      }
      
      cnetResult = clusterProfiler::filter(result, Count > 1)
      
      if (nrow(cnetResult) > 0)
      {
        p <- enrichplot::cnetplot(cnetResult, color.params = list(foldChange = testExpression))
        save_plot(p, paste(outfolder,"enrichment", comparison, perf_test, "cnetplot", sep="/"), fig.width=12, fig.height=12, save.data=FALSE)
      }
      
      result2 <- enrichplot::pairwise_termsim(result)
      
      ncluster = min(c(5, nrow(result)))
      
      if (ncluster  <= 1)
      {
        next()
      }
      
      p <- enrichplot::treeplot(result2, nCluster=ncluster)
      save_plot(p, paste(outfolder,"enrichment", comparison, perf_test, "treeplot", sep="/"), fig.width=8, fig.height=12, save.data=FALSE)
      
      p <- enrichplot::emapplot(result2)
      save_plot(p, paste(outfolder,"enrichment", comparison, perf_test, "emapplot", sep="/"), fig.width=8, fig.height=12, save.data=FALSE)
      
    }
  }
  
}



