# Defining initial patterns

Defining initial patterns

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
