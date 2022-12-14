% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scripts.R
\name{getDEXpressionDF}
\alias{getDEXpressionDF}
\title{Annotate a marker gene data frame with gene expression values}
\usage{
getDEXpressionDF(scdata, markers, assay = "SCT", group.by = NULL)
}
\arguments{
\item{scdata}{Seurat object}

\item{markers}{data frame of marker genes per cluster}

\item{assay}{assay to pull expression values from}

\item{group.by}{grouping of the clusters, same as the one used for markers}
}
\value{
data frame with gene expression values
}
\description{
Annotate a marker gene data frame with gene expression values
}
