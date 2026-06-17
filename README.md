# GAM-clustering

`GAMclust` is an R package for identifying transcriptionally regulated metabolic modules in complex RNA sequencing datasets.

`GAMclust` implements a flexible data processing pipeline designed to analyze three major types of RNA-seq data: bulk, single-cell, and spatial. The pipeline takes as an input a gene expression matrix and a global metabolic network. After initial preprocessing, both of these inputs are integrated by the GAM-clustering algorithm, which identifies metabolic subnetworks (modules) with correlated gene expression behaviors. To aid in the interpretation of the results, several post-processing procedures are carried out, including visualization and annotation of the modules.

Interactive results of GAM-clustering analysis of ImmGen Open Source and Tabula Muris data can be explored [here](http://artyomovlab.wustl.edu/immgen-met/), details are in [Gainullina et al, 2023](https://www.cell.com/cell-reports/fulltext/S2211-1247(23)00057-8).

# Installation

`GAMclust` package can be installed from GitHub:

```
devtools::install_github("alserglab/GAMclust")
```

# Usage

Examples of applying `GAMclust` to different type of data:

-   bulk RNA-seq data ([tutorial](https://rpubs.com/anastasiiaNG/GAMclust_BULK)),
-   single cell RNA-seq data ([tutorial](https://rpubs.com/anastasiiaNG/GAMclust_SC)),
-   spatial RNA-seq data ([tutorial](https://rpubs.com/anastasiiaNG/GAMclust_SPAT)).



