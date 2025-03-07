% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{comparativeVioBoxPlot}
\alias{comparativeVioBoxPlot}
\title{Creates a combined Violin- and Box-Plot with split violins and significance testing}
\usage{
comparativeVioBoxPlot(
  obj.sc,
  feature,
  group.by,
  adj.pval.threshold = 0.05,
  split.by = NULL,
  split.values = NULL,
  dsrCols = NULL,
  onelineLabel = FALSE,
  dot.size = 0,
  min.threshold = 3,
  yStepIncrease = 0.5,
  override = FALSE,
  boxplot_grey = "grey",
  verbose = FALSE
)
}
\arguments{
\item{obj.sc}{Seurat object for plotting}

\item{feature}{genes to plot}

\item{group.by}{name of the meta.data column used for grouping cells}

\item{adj.pval.threshold}{adjusted p-value threshold to show result}

\item{split.by}{name of the meta.data column used for dodging violins}

\item{split.values}{values in split.by to consider}

\item{dsrCols}{named list of colors}

\item{onelineLabel}{significance-test results in one line}

\item{min.threshold}{minimum number of cells in group}

\item{yStepIncrease}{y distance between statistics}

\item{override}{ignore min.threshold}

\item{boxplot_grey}{fill color of boxplot}

\item{pt.size}{size of the single value dots}
}
\value{
ggplot2 object
}
\description{
Creates a combined Violin- and Box-Plot with split violins and significance testing
}
