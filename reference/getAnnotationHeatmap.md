# Heatmap for annotated pathways

Heatmap for annotated pathways

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
