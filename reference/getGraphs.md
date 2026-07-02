# Get files with modules' graphs

Get files with modules' graphs

## Usage

``` r
getGraphs(
  modules,
  network.annotation,
  metabolites.annotation,
  seed.for.layout = 22,
  work.dir = work.dir
)
```

## Arguments

- modules:

  Metabolic modules.

- network.annotation:

  Metabolic network annotation.

- metabolites.annotation:

  Metabolites annotation.

- seed.for.layout:

  Random seed required for gatom to draw metabolic graph. each atom to
  corresponding metabolite's name

- work.dir:

  Working directory where results should be saved.

## Value

Results of this function can be seen in work.dir (.pdf, .png and .xgmml
files).
