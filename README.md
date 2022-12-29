# PLO(SC)²
[![CircleCI](https://circleci.com/gh/mjoppich/PLOSC.svg?style=shield)](https://circleci.com/gh/mjoppich/PLOSC)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/mjoppich/PLOSC)

## Plots and Scripts for scRNA-seq analysis

scRNA-seq analysis becomes a standard technique to study biological systems.
With decreasing costs for scRNA-seq experiments, experiments also become more and more complex.
While the typical scRNA-seq analysis frameworks provide functionalities for the analysis of even such data sets, the workflow involves the use of a chain of provided functions.
Moreover, default plots already provide specific insight into such complex data sets, but should also be enhanced, such that camera-ready fully interpretable plots are provided.

We thus describe here a collection of plotting and analysis scripts for use in Seurat-based scRNA-seq data analyses.
We first provide a collection of script blocks which allows for an easy basic analysis of scRNA-seq from Seurat object creation, filtering over dataset integration.
Subsequently, we provide code blocks for the easy differential analysis of the obtained data sets, including visualizations.
Finally, several visualizations not common to any scRNA-seq analysis framework are presented, such as the enhanced Heatmap and DotPlot.
These, particularly, allow the user to specify how the shown values should be scaled, allowing the creation of condition-wise plots.

With the PLO(SC)² framework the data analysis of scRNA-seq experiments becomes more stream-lined.
This way, fellow researchers can directly apply the methods on their data.
