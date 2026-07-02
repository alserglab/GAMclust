# Heatmap for annotated pathways

Builds a heatmap summarizing pathway annotations across all identified
modules. The visualization displays the proportion of module genes
associated with each enriched pathway, enabling rapid comparison of
functional themes between modules. Filtering options allow users to
restrict displayed pathways based on adjusted p-values and gene overlap
thresholds. The resulting heatmap provides a compact overview of the
functional landscape revealed by GAM-clustering.

## Usage

``` r
getAnnotationHeatmap(
  work.dir,
  padj.threshold = Inf,
  threshold = 0.05,
  file_name = "Modules_heatmap.png"
)
```

## Arguments

- work.dir:

  Directory with gene and pathways tables

- padj.threshold:

  Set threshold for p-adjusted

- threshold:

  Minimal percent of genes in the pathway

- file_name:

  Name of file with heatmap

## Value

Results of this function can be seen in work.dir (png files)
