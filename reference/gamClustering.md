# GAM-clustering analysis

Performs the core GAM-clustering procedure to identify transcriptionally
coordinated metabolic modules. Starting from initial expression
patterns, the algorithm iteratively scores genes, solves SGMWCS
optimization problems on the metabolic network, and updates
module-specific expression profiles until convergence. Additional
refinement steps merge highly similar modules, remove uninformative
patterns, and control module size. The function returns final metabolic
modules, scored networks, module expression patterns, and detailed
iteration statistics.

## Usage

``` r
gamClustering(
  E.prep,
  network.prep,
  cur.centers,
  start.base = 0.5,
  base.dec = 0.1,
  max.module.size = 50,
  cor.threshold = 0.8,
  p.adj.val.threshold = 0.001,
  batch.solver = seq_batch_solver(solver),
  work.dir,
  show.intermediate.clustering = TRUE,
  verbose = TRUE,
  collect.stats = TRUE,
  reference.patterns = NULL
)
```

## Arguments

- E.prep:

  Expression matrix after the
  [`prepareData()`](https://github.com/alserglab/GAMclust/reference/prepareData.md)
  function.

- network.prep:

  Network edge table driven from
  [`prepareNetwork()`](https://github.com/alserglab/GAMclust/reference/prepareNetwork.md)
  function.

- cur.centers:

  Initial patterns produced by
  [`preClustering()`](https://github.com/alserglab/GAMclust/reference/preClustering.md)
  function.

- start.base:

  The parameter which influences modules sizes.

- base.dec:

  The value controlling how strongly `base` parameter should be reduced
  if some module's size is bigger that `max.module.size`. The update
  rule is: `base <- base * (1 - base.dec)`. Detaulf: `0.1`.

- max.module.size:

  Maximal number of unique genes in the final module.

- cor.threshold:

  Threshold for correlation between module patterns.

- p.adj.val.threshold:

  Padj threshold of geseca score for final modules.

- batch.solver:

  Solver for SGMWCS problem.

- work.dir:

  Working directory where results should be saved.

- show.intermediate.clustering:

  Whether to show or not heatmap of intermediate clusters.

- verbose:

  Verbose running.

- collect.stats:

  Whether to save or not running statistics.

- reference.patterns:

  Matrix of reference patterns to track correlation of centers against.
  Pattern per row. Number of columns should be the sames as in E.prep.

## Value

results\$modules – Metabolic modules.

results\$nets – Scored networks.

results\$patterns.pos – Modules' patterns (genes with positive score
only considered).

results\$patterns.all – Modules' patterns (all genes considered).

results\$iter.stats – Statistics from iterations.

## Parameter tuning

There is a set of parameters which determine the size and number of the
final modules. We recommend starting with the default settings, but the
following adjustments may be useful:

- If you consider final modules to be too small or too big and it
  complicates interpretation for you, you can either increase or reduce
  by 10 units the max.module.size parameter.

- If among final modules you consider presence of any modules with too
  similar patterns, you can reduce by 0.1 units the cor.threshold
  parameter.

- If among final modules you consider presence of any uninformative
  modules, you can reduce by 10 times the p.adj.val.threshold parameter.
