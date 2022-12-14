% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scripts.R
\name{makeSeuratObj}
\alias{makeSeuratObj}
\title{Turns a feature-cell-matrix into a Seurat objects with called mt/rp/rps/rpl-content}
\usage{
makeSeuratObj(matrix, proj, pl)
}
\arguments{
\item{matrix}{feature-cell-matrix of the sample}

\item{proj}{project name for the seurat object}

\item{pl}{patternlist for mt-content and rp/rps/rpl-content}
}
\value{
Seurat object
}
\description{
Turns a feature-cell-matrix into a Seurat objects with called mt/rp/rps/rpl-content
}
