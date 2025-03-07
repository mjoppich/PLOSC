# PLO(SC)²
[![CircleCI](https://circleci.com/gh/mjoppich/PLOSC.svg?style=shield)](https://circleci.com/gh/mjoppich/PLOSC)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/mjoppich/PLOSC)

## Plots and Scripts for scRNA-seq analysis

Background scRNA-seq analysis has become a standard technique for studying biological systems.
As costs decrease, scRNA-seq experiments become increasingly complex. While typical scRNA-seq
analysis frameworks provide basic functionality to analyze such data sets, downstream analysis and
visualization become a bottleneck. Standard plots are not always suitable to provide specific insight into
such complex data sets and should be extended to provide camera-ready, meaningful plots.

Results With PLO(SC)², a collection of plotting and analysis scripts for use in Seurat-based scRNA-seq
data analyses is presented, which are accessible for custom script-based analyses or within an R shiny
app. The analysis scripts mainly provide a collection of code blocks which enable a comfortable basic
analysis of scRNA-seq data from Seurat object creation, filtering, and over data set integration in less
than 10 function calls. Subsequently, code blocks for performing differential and enrichment analyses and
corresponding visualizations are provided. Finally, several enhanced visualizations are provided, such as
the enhanced Heatmap, DotPlot and comparative Box-/Violin plots. These, particularly, allow the user to
specify how the shown values should be scaled, allowing the accurate creation of condition-wise plots.

Conclusion With the PLO(SC)² framework data analysis of scRNA-seq experiments is performed
more comfortable and stream-lined, while visualizations are enhanced to be suitable for interpreting
complex datasets. The PLO(SC)² scripts are available from GitHub, including a notebook showing how
PLO(SC)² is applied within a script-based analysis, and an R shiny app.

![plosc_app](https://github.com/user-attachments/assets/5fb25448-b6fb-464e-9f11-02a13212d6d2)


## Install

Before installing PLO(SC)², make sure the following Bioconductor dependencies are installed:

    BiocManager::install(c("biomaRt", "clusterProfiler", "ReactomePA", "org.Hs.eg.db", "org.Mm.eg.db", "ComplexHeatmap", "enrichplot", "EnhancedVolcano"))
  
Then you can install PLO(SC)² using remotes:

    remotes::install_github("mjoppich/PLOSC")
