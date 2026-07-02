# Prepare gene expression matrix

Prepare gene expression matrix

## Usage

``` r
prepareData(
  E,
  gene.id.type = NULL,
  keep.top.genes = 12000,
  use.PCA = TRUE,
  use.PCA.n = 50,
  repeats = seq_len(ncol(E)),
  network.annotation
)
```

## Arguments

- E:

  Expression matrix with rownames as gene symbols.

- gene.id.type:

  Gene ID type.

- keep.top.genes:

  Which top of the most expressed genes to keep for the further
  analysis.

- use.PCA:

  Whether to reduce matrix dimensionality by PCA or not.

- repeats:

  Here you may collapse biological replicas by providing vector with
  repeated sample names

- network.annotation:

  Metabolic network annotation.

## Value

Expression matrix prepared for the analysis.
