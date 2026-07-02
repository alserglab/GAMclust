# Compare two GAM-clustering runs

Compare two GAM-clustering runs

## Usage

``` r
modulesSimilarity(
  dir1,
  dir2,
  name1 = basename(dir1),
  name2 = basename(dir2),
  same.data = FALSE,
  use.genes.with.pos.score = FALSE,
  work.dir = getwd(),
  plot.height = 10,
  plot.width = NULL,
  file.name = sprintf("%s_VS_%s.png", name1, name2)
)
```

## Arguments

- dir1:

  Folder with GAM-clustering results of the first run.

- dir2:

  Folder with GAM-clustering results of the second run.

- name1:

  Name of the first run.

- name2:

  Name of the second run.

- same.data:

  Whether two runs were made on the same data.

- use.genes.with.pos.score:

  Whether to build figure considering positive genes only.

- work.dir:

  Folder where final figure should be saved in.

- plot.height:

  Height of the resulting plot in inches.

- plot.width:

  Width of the plot in inches.

- file.name:

  Name of the final figure.

## Value

Results of this function can be seen in work.dir (.pdf, .png and .xgmml
files).
