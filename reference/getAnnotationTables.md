# Get files with modules' annotations

Get files with modules' annotations

## Usage

``` r
getAnnotationTables(network.annotation, nets, work.dir, padj.threshold = Inf)
```

## Arguments

- network.annotation:

  Metabolic network annotation.

- nets:

  Scored networks.

- work.dir:

  Working directory where files with module genes are (results will be
  saved here as well).

- padj.threshold:

  Threshold, adjusted p-value, for pathways.

## Value

Results of this function can be seen in work.dir (.tsv files).
