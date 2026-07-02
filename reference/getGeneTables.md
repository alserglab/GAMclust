# Get files with modules' genes

Get files with modules' genes

## Usage

``` r
getGeneTables(
  modules,
  nets,
  patterns,
  gene.exprs,
  network.annotation,
  work.dir = work.dir
)
```

## Arguments

- modules:

  Metabolic modules.

- nets:

  Scored networks.

- patterns:

  Patterns of metabolic modules.

- gene.exprs:

  Gene expression.

- network.annotation:

  Metabolic network annotation.

- work.dir:

  Working directory where results should be saved.

## Value

m.\*.genes.tsv – module genes.

m.\*.notInModule.genes.tsv – genes with positive score not included into
module.

m.\*.complete.genes.tsv – top 300 of all genes sorted by correlation
value.

Results of this function can be seen in work.dir (three .tsv files for
each module with gene lists).
