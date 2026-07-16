# Algorithm overview

`GAMclust` identifies coordinated metabolic modules by combining gene
co-expression analysis with the topology of the KEGG metabolic network.
Rather than clustering genes directly, the algorithm iteratively refines
representative expression patterns and identifies connected subnetworks
whose genes best match these patterns.

## Initialization

The algorithm starts by clustering all metabolic genes in sample using
*k*-means with default `k = 32`. The resulting cluster centroids serve
as the initial expression patterns.

For every gene *g*_(*i*) and every pattern *c*_(*j*), a score is
calculated as

``` math
\mathrm{score}(g_i,c_j) =
-\log\left(\frac{d(g_i,c_j)}{d'(g_i,c_j)}\right)
```

where

``` math
d(g_i,c_j)=1-\mathrm{cor}(g_i,c_j)
```
is the correlation-based distance between the gene expression profile
*i* and pattern *j*,

``` math
d(g_i,c_0)=\mathrm{base}
```

is the initially fixed distance to a fictitious pattern, and

``` math
d'\left( g_{i}, c_{j} \right) = \min_{k \neq j,\; k \in (0, M)} \left( d\left( g_{i}, c_{k} \right) \right)
```

is the minimum distance to any other pattern, including the fictitious
one.

A positive score indicates that the gene *i* is more strongly associated
with pattern *j* than with any other pattern.

These scores are assigned as edge weights in the KEGG metabolic network.
For every expression pattern, the algorithm identifies a maximum-weight
connected subgraph using an SGMWCS (signal generalized maximum-weight
connected subgraph) solver. Each connected subgraph represents a
candidate metabolic module.

## Iterative refinement & convergence criteria

Having obtained these initial candidate modules from the SGMWCS solver,
the algorithm proceeds to iteratively refine them. At each iteration,
every expression pattern is recalculated as the averaged expression
profile of all genes contributing positive scores to its module. The
updated patterns are subsequently used to recompute gene scores and
identify new graph-based modules. This process continues until a stable
set of modules is achieved according to explicit convergence (a matrix
of scaled patterns is close to a corresponding matrix from one of the
previous iterations with at most 0.01 element-wise absolute difference).

After convergence, additional refinement steps are performed whenever
necessary to ensure that the final modules remain biologically
meaningful and easy to interpret:

- First, excessively large modules are avoided. If any module contains
  more genes than specified by `max.module.size`, the parameter `base`
  is decreased by `base.dec` (10% by default), making inclusion of genes
  into modules more strict.

- Second, highly similar modules are subject to adjustment. If the
  correlation between two module expression patterns exceeds
  `cor.threshold`, the two most similar modules are merged.

- Third, the answer should not include uninformative modules. Module
  informativeness is evaluated using a GESECA (GSEA-like) enrichment
  score which captures coordinated changes across the entire
  transcriptome. Module that does not pass the significance threshold
  (`p.adj.val.threshold`) is considered uninformative and is merged with
  the any other most similar module which it has the greatest gene
  intersection with. If there is no such similar module, the
  uninformative module is removed.

These refinement steps are repeated until all modules satisfy the above
quality criteria.

## Parameters tuning

Thus, the characteristics of the resulting modules and the rules
governing module refinement are determined by the set of pre-defined
parameters.

    cur.centers <- preClustering(...,
                                 initial.number.of.clusters = 32)

    results <- gamClustering(...,
                             start.base = 0.5,
                             base.dec = 0.1,
                             max.module.size = 50,
                             cor.threshold = 0.8,
                             p.adj.val.threshold = 1e-5)

We recommend starting with these default settings, as they provide a
good balance between module size, uniqueness of expression patterns, and
biological interpretability.

However, below we provide a set of instructions on how to tune
parameters depending on certain characteristics of the output modules:

- if the resulting modules are consistently too small or too large for
  convenient biological interpretation, correspondingly increase or
  decrease `max.module.size` by approximately 10 genes;

- if several final modules exhibit highly similar expression patterns,
  decrease `cor.threshold` by approximately 0.1 to encourage additional
  merging;

- if some final modules appear biologically uninformative, decrease
  `p.adj.val.threshold` by one order of magnitude (for example, from
  `1e-5` to `1e-6`) to apply a more stringent informativeness criterion.

## Algorithm pseudocode

Algorithm: GAM-clustering

Input:

• Gene expression matrix *E* of size $`n \times m`$  
• Network edge table *network*  
• List of $`k`$ initial patterns *cur.patterns* of length $`m`$  
• Parameters: *start.base*, *base.dec*, *max.module.size*,
*cor.threshold*, *p.adj.val.threshold*

Result:

• Final metabolic *modules*

  

*base* ← *start.base*  

while true do

for $`j \in \mathbb{N}`$do

*k’* ← \|*cur.patterns*\|  
*revisions* ← {}  
*dist* ← calculateCorrelationDistances(*E*, *cur.patterns*)  

for $`i \in \{1,\ldots,k'\}`$do

*scores*\[*i*\] ← calculateGeneScores(*dist*, *i*, *base*)  
*mwcs*\[*i*\] ← makeMwcsInstance(*network*, *scores*\[*i*\])  
*modules*\[*i*\] ← solveMwcs(*mwcs*\[*i*\])  
*new.patterns*\[*i*\] ← calculatePatterns(*modules*\[*i*\])

end  

*cur.patterns* ← *new.patterns*  
*revisions*\[*j*\] ← *cur.patterns*  

if not changedSignificantly(*cur.patterns*, *revisions*\[$`1..(j-1)`$\])
then

break

end

end  

*hasBig* ← areThereAnyTooBigModules(*modules*, *max.module.size*)  
*hasCorr*, *mc*₁, *mc*₂ ← areThereAnyCorrelatedModules(*modules*,
*cor.threshold*)  
*hasUninform* ← areThereAnyUninformativeModules(*modules*,
*p.adj.val.threshold*)  

if *hasBig* then

*base* ← *base* − *base.dec* × *base*

end  

if *hasCorr* then

*cur.patterns* ← mergeTwoModules(*cur.patterns*, *mc*₁, *mc*₂)

end  

if *hasUninform* then

*maxOverlap*, *mu*₁, *mu*₂ ← checkOverlap(*modules*)  

if *maxOverlap* \> 0 then

*cur.patterns* ← mergeTwoModules(*cur.patterns*, *mu*₁, *mu*₂)

else

*cur.patterns*\[*mu*₁\] ← NULL

end

end  

if not any(*hasBig*, *hasCorr*, *hasUninform*) then

break

end

end  

return *modules*
