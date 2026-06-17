# Prepare network

Constructs a metabolic interaction network tailored to the analyzed
dataset. The function integrates network topology with
gene–enzyme–reaction annotations, removes filtered metabolites, and
retains only genes present in the processed expression matrix. Networks
can be represented at either the metabolite or atom level. The final
output corresponds to the largest connected component of the metabolic
network and serves as the graph structure for GAM-clustering.

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

  Expression matrix after the \`prepareData()\` function.

- network:

  Metabolic network.

- topology:

  Vertices can be represented either as \`metabolites\`, either as
  \`atoms\`.

- met.to.filter:

  Metabolites that should not be used as connections in the module.

- network.annotation:

  Metabolic network annotation.

- gene2reaction.extra:

  For a combined network: supplementary file with genes that either come
  from proteome or are not linked to a specific enzyme.

## Value

Edges of the final network.
