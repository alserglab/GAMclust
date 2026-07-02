# Prepare network

Prepare network

## Usage

``` r
prepareNetwork(
  E,
  network,
  topology = c("metabolites", "atoms"),
  met.to.filter = data.table::fread(system.file("mets2mask.lst", package =
    "GAMclust"))$ID,
  network.annotation,
  gene2reaction.extra = NULL
)
```

## Arguments

- E:

  Expression matrix after the
  [`prepareData()`](https://github.com/alserglab/GAMclust/reference/prepareData.md)
  function.

- network:

  Metabolic network.

- topology:

  Vertices can be represented either as `metabolites`, either as
  `atoms`.

- met.to.filter:

  Metabolites that should not be used as connections in the module.

- network.annotation:

  Metabolic network annotation.

- gene2reaction.extra:

  For a combined network: supplementary file with genes that either come
  from proteome or are not linked to a specific enzyme.

## Value

Edges of the final network.
