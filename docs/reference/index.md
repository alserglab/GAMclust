# Package index

## Data Preparation & Clustering

Core functions for configuring networks and executing the main GAMclust
algorithms.

- [`prepareData()`](https://github.com/alserglab/GAMclust/reference/prepareData.md)
  : Prepare gene expression matrix
- [`prepareNetwork()`](https://github.com/alserglab/GAMclust/reference/prepareNetwork.md)
  : Prepare network
- [`preClustering()`](https://github.com/alserglab/GAMclust/reference/preClustering.md)
  : Defining initial patterns
- [`gamClustering()`](https://github.com/alserglab/GAMclust/reference/gamClustering.md)
  : GAM-clustering analysis
- [`seq_batch_solver()`](https://github.com/alserglab/GAMclust/reference/seq_batch_solver.md)
  : Function for batch solving SGMWCS problems

## Visualization & Downstream Analysis

Functions for exploring network graphs, generating annotation tables,
and plotting heatmaps.

- [`getAnnotationHeatmap()`](https://github.com/alserglab/GAMclust/reference/getAnnotationHeatmap.md)
  : Heatmap for annotated pathways
- [`getAnnotationTables()`](https://github.com/alserglab/GAMclust/reference/getAnnotationTables.md)
  : Get files with modules' annotations
- [`getGeneTables()`](https://github.com/alserglab/GAMclust/reference/getGeneTables.md)
  : Get files with modules' genes
- [`getGraphs()`](https://github.com/alserglab/GAMclust/reference/getGraphs.md)
  : Get files with modules' graphs
- [`modulesSimilarity()`](https://github.com/alserglab/GAMclust/reference/modulesSimilarity.md)
  : Compare two GAM-clustering runs
