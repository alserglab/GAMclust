# Defining initial patterns

Generates an initial set of transcriptional patterns that serve as
starting points for module detection. The function first restricts the
expression matrix to metabolic genes represented in the network and then
performs clustering to identify preliminary expression signatures. By
default, cluster centers are obtained using k-means clustering. These
initial patterns are subsequently refined by the GAM-clustering
algorithm.

## Usage

``` r
preClustering(
  E.prep,
  network.prep,
  initial.number.of.clusters = 32,
  network.annotation,
  use.ICA = FALSE
)
```

## Arguments

- E.prep:

  Expression matrix after the
  [`prepareData()`](https://github.com/alserglab/GAMclust/reference/prepareData.md)
  function.

- network.prep:

  Network edge table driven from
  [`prepareNetwork()`](https://github.com/alserglab/GAMclust/reference/prepareNetwork.md)
  function.

- initial.number.of.clusters:

  The number of clusters for the initial approximation of modules.

- network.annotation:

  Metabolic network annotation.

## Value

Initial patterns.
