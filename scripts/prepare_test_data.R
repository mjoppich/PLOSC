devtools::install_github('satijalab/seurat-data')
library(SeuratData)
library(Seurat)
library(patchwork)

InstallData("ifnb")

ifnb <- LoadData("ifnb")
ifnb = UpdateSeuratObject(object = ifnb)
ifnb

ifnb[["RNA"]] <- split(ifnb[["RNA"]], f = ifnb$stim)

ifnb <- NormalizeData(ifnb)
ifnb <- FindVariableFeatures(ifnb)
ifnb <- ScaleData(ifnb, features=rownames(ifnb))
ifnb <- RunPCA(ifnb)

ifnb <- IntegrateLayers(object = ifnb, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)
# re-join layers after integration
ifnb[["RNA"]] <- JoinLayers(ifnb[["RNA"]])

ifnb <- FindNeighbors(ifnb, reduction = "integrated.cca", dims = 1:30)
ifnb <- FindClusters(ifnb, resolution = 1)

ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "integrated.cca")
DimPlot(ifnb, reduction = "umap", group.by = c("stim", "seurat_annotations"))


#ifnb = NormalizeData(ifnb)
#ifnb = PLOSC::preprocessIntegrated(ifnb, useAssay = "RNA", with.hto = FALSE)

saveRDS(ifnb, "./ifnb.Rds")
ifnb = readRDS("./ifnb.Rds")

library(tidyr)
library(dplyr)
markers = PLOSC::makeDEResults(ifnb, group.by="stim", assay="RNA", test="wilcox")

head(markers[markers$clusterID == "STIM" & markers$avg_log2FC > 0 & markers$mean.cluster > 3 & markers$mean.bg > 2,] %>% arrange(desc(pct.1-pct.2)) )
head(markers[markers$clusterID == "CTRL" & markers$avg_log2FC > 0 & markers$mean.cluster > 3 & markers$mean.bg > 2,] %>% arrange(desc(pct.1-pct.2)) )



relGenes = c("IL8", "S100A8", "LGALS1", "VIM")

plotElems = list()
plotElems[["CTRL"]] = list("label"="CTRL", "cells"=PLOSC::cellIDForClusters(ifnb, "stim", c("CTRL")))
plotElems[["STIM"]] = list("label"="STIM", "cells"=PLOSC::cellIDForClusters(ifnb, "stim", c("STIM")))

p=PLOSC::enhancedDotPlot(ifnb, plotElems, featureGenes = relGenes, scale.by = "ALL", group.by = "seurat_annotations")
PLOSC::save_plot(p, "plots/dotplot_split_all", fig.width = 14, fig.height=8)

p=PLOSC::enhancedDotPlot(ifnb, plotElems, featureGenes = relGenes, scale.by = "GLOBAL", group.by = "seurat_annotations")
PLOSC::save_plot(p, "plots/dotplot_split_global", fig.width = 14, fig.height=8)

p=PLOSC::enhancedDotPlot(ifnb, plotElems, featureGenes = relGenes, scale.by = "GROUP", group.by = "seurat_annotations")
PLOSC::save_plot(p, "plots/dotplot_split_group", fig.width = 14, fig.height=8)


p=PLOSC::enhancedDotPlot(ifnb, plotElems, featureGenes = relGenes, scale.by = "FEATURE", group.by = "seurat_annotations")
PLOSC::save_plot(p, "plots/dotplot_split_feature", fig.width = 14, fig.height=8)


groupby="seurat_annotations"
allGroups = ifnb@meta.data[ , groupby]
allGroups = unique(as.character(allGroups))

plotGOIs = list()
plotGOIs[["all"]] = list(clusters=allGroups, genes=relGenes)

p=PLOSC::enhancedHeatMap( ifnb, plotGOIs, group.by=groupby, include_all_clusters=TRUE, scale.by="GLOBAL")
PLOSC::save_plot(p, "plots/heatmap_global", fig.width = 14, fig.height=8, save.data=FALSE)

p=PLOSC::enhancedHeatMap( ifnb, plotGOIs, group.by=groupby, include_all_clusters=TRUE, scale.by="ALL")
PLOSC::save_plot(p, "plots/heatmap_all", fig.width = 14, fig.height=8, save.data=FALSE)


p=PLOSC::makeComplexExprHeatmapSplit(obj.in = ifnb, plot_gois = plotGOIs, split.by = "stim", group.by = groupby, scale.by = "GLOBAL",include_all_clusters=TRUE)
PLOSC::save_plot(p, "plots/heatmap_split_global", fig.width = 14, fig.height=8, save.data=FALSE)


p=PLOSC::makeComplexExprHeatmapSplit(obj.in = ifnb, plot_gois = plotGOIs, split.by = "stim", group.by = groupby, scale.by = "ALL",include_all_clusters=TRUE)
PLOSC::save_plot(p, "plots/heatmap_split_all", fig.width = 14, fig.height=8, save.data=FALSE)


DimPlot(ifnb)

PLOSC::comparativeVioBoxPlot(ifnb, "S100A8", group.by="seurat_annotations", split.by = "stim")


