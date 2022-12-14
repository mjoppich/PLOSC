% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{enhancedDotPlot}
\alias{enhancedDotPlot}
\title{Creates an enhanced DotPlot}
\usage{
enhancedDotPlot(
  scobj,
  plotElems,
  featureGenes = c(""),
  group.by = "cellnames_manual",
  col.min = -3,
  col.max = 3,
  cols = c("grey", "blue"),
  title = "",
  scale.by = "GROUP",
  rotate.x = F,
  abundance.perelem = FALSE,
  assay = "RNA"
)
}
\arguments{
\item{scobj}{Seurat object for plotting}

\item{plotElems}{list of list(cells,label) where each entry defines one condition/split of the plot}

\item{featureGenes}{genes to plot}

\item{group.by}{name of the meta.data column used for grouping cells}

\item{col.min}{lower bound of the scale}

\item{col.max}{upper bound of the scale}

\item{cols}{color for the expression values}

\item{title}{title of the plot}

\item{scale.by}{how to show/scale the expression values (GROUP; FEATURE; ALL; GLOBAL)}

\item{rotate.x}{whether to rotate x-axis labels}

\item{abundance.perelem}{whether the cell group abundance is to be calculated on the global Seurat object, or per group}

\item{assay}{which assay to use for retrieving the expression values}
}
\value{
ggplot2 object
}
\description{
Creates an enhanced DotPlot
}
