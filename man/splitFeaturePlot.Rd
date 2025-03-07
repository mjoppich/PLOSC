% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{splitFeaturePlot}
\alias{splitFeaturePlot}
\title{SplitFeaturePlot - creates a Featureplot split by a given condition and equal scales over all subplots}
\usage{
splitFeaturePlot(
  obj,
  feature,
  split.by,
  title = NULL,
  limits = c(-1, 1),
  reduction = "umap",
  ncol = NULL,
  low = "lightgrey",
  high = "blue",
  mid = "white",
  mirrorLimits = TRUE
)
}
\arguments{
\item{obj}{Seurat object}

\item{feature}{which feature to show}

\item{split.by}{by which condition/meta-data column to split}

\item{title}{Title of the plot}

\item{limits}{limits of the features}

\item{reduction}{which reduction to plot}

\item{ncol}{how many columns the plot should have}

\item{low}{color for low expression}

\item{high}{color for high expression}

\item{mid}{medium color}

\item{mirrorLimits}{whether limits should be mirrored. useful for scaled data.}
}
\description{
SplitFeaturePlot - creates a Featureplot split by a given condition and equal scales over all subplots
}
