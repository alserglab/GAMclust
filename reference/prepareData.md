# Prepare gene expression matrix

Prepares a gene expression matrix for GAM-clustering analysis. The
function optionally converts gene identifiers to the annotation format
used by the metabolic network, removes duplicates, retains the most
highly expressed genes, and standardizes expression values. Biological
replicates can be collapsed by averaging, and dimensionality reduction
can be applied using PCA. The resulting matrix contains only processed
features suitable for downstream network-based module discovery.

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
