# Get files with modules' annotations

Annotates metabolic modules using pathway enrichment analysis. For each
module, over-representation analysis is performed against the pathway
collection contained in the network annotation, followed by pathway
redundancy reduction. The function exports pathway tables containing
significantly enriched biological processes and their associated module
genes. These annotations help connect discovered modules to known
metabolic functions.

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
