# GAM-clustering Tutorial Notebook For Reanalysis Of \[Metabolic Tabula Muris Senis\](http://artyomovlab.wustl.edu/immgen-met/) Single Cell RNA-Seq Data

Install GAMclust package:

``` r
devtools::install_github("alserglab/GAMclust")
```

``` r
library(GAMclust)
library(gatom)
library(mwcsr)
library(fgsea)
library(data.table)
library(Seurat)
library(futile.logger)

set.seed(42)
```

### Preparing working environment

First, please load and initialize all objects required for
GAM-clustering analysis:

1.  load metabolic network and its metabolites annotation. We provide
    two networks: [KEGG](https://www.genome.jp/kegg/) and combined
    network that includes [KEGG](https://www.genome.jp/kegg/),
    [Rhea](https://www.rhea-db.org/) and transport reactions:

(1.1.) load [KEGG](https://www.genome.jp/kegg/) metabolic network
[`network.kegg.rds`](http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/network.kegg.rds)
and its metabolites annotation
[`met.kegg.db.rds`](http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/met.kegg.db.rds)

``` r
# KEGG network:
network <- readRDS(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/network.kegg.rds"))
metabolites.annotation <- readRDS(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/met.kegg.db.rds"))
```

(1.2.) or load combined metabolic network
[`network.combined.rds`](http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/network.combined.rds),
its metabolites annotation
[`met.combined.db.rds`](http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/met.combined.db.rds)
and species-specific list of genes that either come from proteome or are
not linked to a specific enzyme
[`gene2reaction.combined.mmu.eg.tsv`](http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/gene2reaction.combined.mmu.eg.tsv)
for mouse and
[`gene2reaction.combined.hsa.eg.tsv`](http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/gene2reaction.combined.hsa.eg.tsv)
for human data;

``` r
# combined network (KEGG+Rhea+transport reactions):
network <- readRDS(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/network.combined.rds"))
metabolites.annotation <- readRDS(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/met.combined.db.rds"))
gene2reaction.extra <- data.table::fread("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/gene2reaction.combined.mmu.eg.tsv", colClasses="character")
```

2.  load species-specific network annotation:
    [`org.Hs.eg.gatom.anno`](http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/org.Hs.eg.gatom.anno.rds)
    for human data or
    [`org.Mm.eg.gatom.anno`](http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/org.Mm.eg.gatom.anno.rds)
    for mouse data;

``` r
network.annotation <- readRDS(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/org.Mm.eg.gatom.anno.rds"))
```

3.  load provided list of metabolites that should not be considered
    during the analysis as connections between reactions (e.g., CO2,
    HCO3-, etc);

``` r
met.to.filter <- data.table::fread(system.file("mets2mask.lst", package="GAMclust"))$ID
```

4.  initialize SMGWCS solver:

(4.1.) we recommend to use here either heuristic relax-and-cut solver
`rnc_solver` from [`mwcsr`](https://github.com/ctlab/mwcsr) package,

``` r
solver <- mwcsr::rnc_solver()
```

(4.2.) either proprietary [CPLEX
solver](https://www.ibm.com/products/ilog-cplex-optimization-studio)
(free for academy);

``` r
cplex.dir <- "/opt/ibm/ILOG/CPLEX_Studio1271"
```

``` r
solver <- mwcsr::virgo_solver(cplex_dir = cplex.dir)
```

5.  set working directory where the results will be saved to.

``` r
work.dir <- "results_sc"
dir.create(work.dir, showWarnings = F, recursive = T)
```

6.  TEMPORARY: collecting logs while developing the tool.

``` r
stats.dir <- file.path(work.dir, "stats")
dir.create(stats.dir, showWarnings = F, recursive = T)

setup_logger <- function(log.file.path, logger.name = "stats.logger") {
  file.appender <- appender.file(log.file.path)
  console.appender <- appender.console()
  combined.appender <- function(line) {
    file.appender(line)
    console.appender(line)
  }
  flog.appender(combined.appender, name = logger.name)
  flog.threshold(TRACE, name = logger.name)
}

log.file <- file.path(stats.dir, "log.txt")
setup_logger(log.file.path = log.file, logger.name = "stats.logger")
```

### Preparing objects for the analysis

#### Preparing data

GAMclust works with bulk, single cell and spatial RNA-seq data.

This vignette shows how to process single cell RNA-seq data on the
example of [Tabula Muris
Senis](http://artyomovlab.wustl.edu/immgen-met/) data reanalysis.

For single cell data, take 6,000-12,000 genes for the GAM-clustering
analysis. To do this while [preprocessing data with Seurat
pipeline](https://satijalab.org/seurat/articles/sctransform_vignette),
set `variable.features.n = 12000` in
[`SCTransform()`](https://satijalab.org/seurat/reference/SCTransform.html)
function. In case of [preprocessing multi-sample
data](https://satijalab.org/seurat/articles/integration_introduction),
set `nfeatures=12000` in
[`SelectIntegrationFeatures()`](https://satijalab.org/seurat/reference/SelectIntegrationFeatures.html)).

Let’s load already preprocessed data.

``` r
seurat_object <- readRDS(url("http://artyomovlab.wustl.edu/publications/supp_materials/GAMclust/tms12k.rds"))

E <- as.matrix(Seurat::GetAssayData(object = seurat_object,
                                    assay = "SCT",
                                    layer = "scale.data"))

nrow(E) # ! make sure this value is in range from 6,000 to 12,000

E[1:3, 1:3]
```

Genes in your dataset may be named as Symbol, Entrez, Ensembl or RefSeq
IDs. One of these names should be specified as value of `gene.id.type`
parameter in
[`prepareData()`](https://github.com/alserglab/GAMclust/reference/prepareData.md).

If you analyse singe cell or spatial RNA-seq data, please set
`use.PCA=TRUE` in
[`prepareData()`](https://github.com/alserglab/GAMclust/reference/prepareData.md).

``` r
E.prep <- prepareData(E = E,
                      gene.id.type = "Symbol",
                      use.PCA = TRUE,
                      use.PCA.n = 50,
                      network.annotation = network.annotation)

E.prep[1:3, 1:3]
```

    #                PC1        PC2         PC3
    # 16176   0.07415203 2.20087518 -0.85846106
    # 278180  1.21894952 2.84602959  1.64928876
    # 17394  -5.31802334 0.07568864  0.08151284

#### Preparing network

The
[`prepareNetwork()`](https://github.com/alserglab/GAMclust/reference/prepareNetwork.md)
function defines the structure of the final metabolic modules.

``` r
network.prep <- prepareNetwork(E = E.prep,
                               network = network,
                               topology = "metabolites",
                               met.to.filter = met.to.filter,
                               network.annotation = network.annotation,
                               gene2reaction.extra = gene2reaction.extra) # for combined network
```

    # INFO [2026-06-15 21:33:34] Global metabolite network contains 6754 edges.
    # INFO [2026-06-15 21:33:34] Largest connected component of this global network contains 1397 nodes and 5139 edges.

#### Preclustering

The
[`preClustering()`](https://github.com/alserglab/GAMclust/reference/preClustering.md)
function defines initial patterns using k-medoids clustering on gene
expression matrix. It is strongly recommended to do initial clustering
with no less than 32 clusters (`initial.number.of.clusters = 32`).

You can visualize the initial heatmap as shown below.

``` r
cur.centers <- preClustering(E.prep = E.prep,
                             network.prep = network.prep,
                             initial.number.of.clusters = 32,
                             network.annotation = network.annotation)
```

    # INFO [2026-06-15 21:33:34] 1123 metabolic genes from the analysed dataset mapped to this component.

``` r
cur.centers[1:3, 1:3]
```

    #         PC1         PC2       PC3
    # 1  4.777154 -0.04738601  1.088851
    # 2  2.852732  0.71344306 -2.404469
    # 3 -1.788987  0.91377131 -1.362991

``` r
pheatmap::pheatmap(
      GAMclust:::normalize.rows(cur.centers),
      cluster_rows=F, cluster_cols=F,
      show_rownames=F, show_colnames=T)
```

![plot of chunk pre-clustering](figure/pre-clustering-1.png)

plot of chunk pre-clustering

### GAM-clustering analysis

Now you have everything prepared for the GAM-clustering analysis.

Initial patterns will be now refined in an iterative process. The output
of
[`gamClustering()`](https://github.com/alserglab/GAMclust/reference/gamClustering.md)
function presents a set of specific subnetworks (also called metabolic
modules) that reflect metabolic variability within a given
transcriptional dataset.

Note, that it may take a long time to derive metabolic modules by
[`gamClustering()`](https://github.com/alserglab/GAMclust/reference/gamClustering.md)
function (tens of minutes).

There is a set of parameters which determine the size and number of your
final modules. We recommend you to start with the default settings,
however you can adjust them based on your own preferences:

1.  If you consider final modules to bee too small or too big and it
    complicates interpretation for you, you can either increase or
    reduce by 10 units the `max.module.size` parameter.

2.  If among final modules you consider presence of any modules with too
    similar patterns, you can reduce by 0.1 units the `cor.threshold`
    parameter.

3.  If among final modules you consider presence of any uninformative
    modules, you can reduce by 10 times the `p.adj.val.threshold`
    parameter.

``` r
results <- gamClustering(E.prep = E.prep,
                         network.prep = network.prep,
                         cur.centers = cur.centers,
                         
                         start.base = 0.5,
                         base.dec = 0.05,
                         max.module.size = 50,

                         cor.threshold = 0.8,
                         p.adj.val.threshold = 1e-5,

                         batch.solver = seq_batch_solver(solver),
                         work.dir = work.dir,
                         
                         show.intermediate.clustering = FALSE,
                         verbose = FALSE,
                         collect.stats = TRUE)
```

    # INFO [2026-06-15 21:33:35] GAM-CLUSTERING starts here.
    # INFO [2026-06-15 21:33:35] [*] Iteration 1
    # INFO [2026-06-15 21:34:59] [*] Iteration 2
    # INFO [2026-06-15 21:36:26] >> max cor exceeded 0.8: 0.86
    # INFO [2026-06-15 21:36:26] [*] Iteration 3
    # INFO [2026-06-15 21:37:49] [*] Iteration 4
    # INFO [2026-06-15 21:39:09] >> max cor exceeded 0.8: 0.82
    # INFO [2026-06-15 21:39:09] [*] Iteration 5
    # INFO [2026-06-15 21:40:28] [*] Iteration 6
    # INFO [2026-06-15 21:41:48] >> max cor exceeded 0.8: 0.82
    # INFO [2026-06-15 21:41:48] [*] Iteration 7
    # INFO [2026-06-15 21:43:03] [*] Iteration 8
    # INFO [2026-06-15 21:44:17] >> max cor exceeded 0.8: 0.83
    # INFO [2026-06-15 21:44:17] [*] Iteration 9
    # INFO [2026-06-15 21:45:27] [*] Iteration 10
    # INFO [2026-06-15 21:46:37] >> max cor exceeded 0.8: 0.85
    # INFO [2026-06-15 21:46:37] [*] Iteration 11
    # INFO [2026-06-15 21:47:46] [*] Iteration 12
    # INFO [2026-06-15 21:48:54] >> max cor exceeded 0.8: 0.82
    # INFO [2026-06-15 21:48:54] [*] Iteration 13
    # INFO [2026-06-15 21:50:00] [*] Iteration 14
    # INFO [2026-06-15 21:51:05] >> max cor exceeded 0.8: 0.82
    # INFO [2026-06-15 21:51:05] [*] Iteration 15
    # INFO [2026-06-15 21:52:07] [*] Iteration 16
    # INFO [2026-06-15 21:53:09] >> max cor exceeded 0.8: 0.84
    # INFO [2026-06-15 21:53:09] [*] Iteration 17
    # INFO [2026-06-15 21:54:08] [*] Iteration 18
    # INFO [2026-06-15 21:55:07] >> max cor exceeded 0.8: 0.82
    # INFO [2026-06-15 21:55:07] [*] Iteration 19
    # INFO [2026-06-15 21:56:04] [*] Iteration 20
    # INFO [2026-06-15 21:56:59] >> max cor exceeded 0.8: 0.84
    # INFO [2026-06-15 21:56:59] [*] Iteration 21
    # INFO [2026-06-15 21:57:51] [*] Iteration 22
    # INFO [2026-06-15 21:58:49] [*] Iteration 23
    # INFO [2026-06-15 21:59:39] [*] Iteration 24
    # INFO [2026-06-15 22:00:31] [*] Iteration 25
    # INFO [2026-06-15 22:01:20] [*] Iteration 26
    # INFO [2026-06-15 22:02:11] [*] Iteration 27
    # INFO [2026-06-15 22:02:57] [*] Iteration 28
    # INFO [2026-06-15 22:03:46] [*] Iteration 29
    # INFO [2026-06-15 22:04:32] [*] Iteration 30
    # INFO [2026-06-15 22:05:12] [*] Iteration 31
    # INFO [2026-06-15 22:05:55] [*] Iteration 32
    # INFO [2026-06-15 22:06:35] [*] Iteration 33
    # INFO [2026-06-15 22:07:12] [*] Iteration 34
    # INFO [2026-06-15 22:07:52] [*] Iteration 35
    # INFO [2026-06-15 22:08:26] [*] Iteration 36
    # INFO [2026-06-15 22:09:02] [*] Iteration 37
    # INFO [2026-06-15 22:09:36] [*] Iteration 38
    # INFO [2026-06-15 22:10:09] [*] Iteration 39
    # INFO [2026-06-15 22:10:40] [*] Iteration 40
    # INFO [2026-06-15 22:11:12] [*] Iteration 41
    # INFO [2026-06-15 22:11:40] [*] Iteration 42
    # INFO [2026-06-15 22:12:10] [*] Iteration 43
    # INFO [2026-06-15 22:12:36] [*] Iteration 44
    # INFO [2026-06-15 22:13:03] [*] Iteration 45
    # INFO [2026-06-15 22:13:27] [*] Iteration 46
    # INFO [2026-06-15 22:13:47] [*] Iteration 47
    # INFO [2026-06-15 22:14:08] [*] Iteration 48
    # INFO [2026-06-15 22:14:26] [*] Iteration 49
    # INFO [2026-06-15 22:14:46] [*] Iteration 50
    # INFO [2026-06-15 22:15:01] [*] Iteration 51
    # INFO [2026-06-15 22:15:18] [*] Iteration 52
    # INFO [2026-06-15 22:15:31] [*] Iteration 53
    # INFO [2026-06-15 22:15:48] GAM-CLUSTERING ends here.

### Visualizing and exploring the GAM-clustering results

Each metabolic module is a connected piece of metabolic network whose
genes expression is correlated across all dataset.

The following functions will help you to visualize and explore the
obtained results.

##### Get graphs of modules:

``` r
getGraphs(modules = results$modules,
          network.annotation = network.annotation,
          metabolites.annotation = metabolites.annotation,
          seed.for.layout = 42,
          work.dir = work.dir)
```

    # Graphs for module 1 are built

    # Graphs for module 2 are built

    # Graphs for module 3 are built

    # Graphs for module 4 are built

    # Graphs for module 5 are built

Example of the graph of the third module:

![plot of chunk unnamed-chunk-5](results_sc/m.3.svg)

plot of chunk unnamed-chunk-5

##### Get gene tables:

The table contains gene list. Each gene has two descriptive values: i)
gene’s correlation value with the modules pattern and ii) gene’s score.
High score means that this gene’s expression is similar to the module’s
pattern and not similar to other modules’ patterns.

``` r
m.gene.list <- getGeneTables(modules = results$modules,
                             nets = results$nets,
                             patterns = results$patterns.pos,
                             gene.exprs = E.prep,
                             network.annotation = network.annotation,
                             work.dir = work.dir)
```

    # Gene tables for module 1 are produced

    # Gene tables for module 2 are produced

    # Gene tables for module 3 are produced

    # Gene tables for module 4 are produced

    # Gene tables for module 5 are produced

Example of the gene table of the third module:

``` r
head(GAMclust:::read.tsv(file.path(work.dir, "m.3.genes.tsv"))) |>
  kableExtra::kable() |>
  kableExtra::kable_styling()
```

| symbol  | Entrez |    score |       cor |
|:--------|-------:|---------:|----------:|
| Msra    | 110265 | 2.857254 | 0.9795676 |
| Elovl7  |  74559 | 2.350795 | 0.9718903 |
| Slc27a4 |  26569 | 2.342289 | 0.9713593 |
| Dgat1   |  13350 | 2.188936 | 0.9680066 |
| Acpp    |  56318 | 1.803953 | 0.9571438 |
| Tst     |  22117 | 1.735225 | 0.9541673 |

##### Get plots of patterns:

``` r
for(i in 1:length(m.gene.list)){
  
  print(fgsea::plotCoregulationProfileReduction(m.gene.list[[i]], 
                                               seurat_object, 
                                               title = paste0("module ", i),
                                               raster = TRUE,
                                               reduction = "pumap"))
}
```

![plot of chunk figures-side](figure/figures-side-1.png)![plot of chunk
figures-side](figure/figures-side-2.png)![plot of chunk
figures-side](figure/figures-side-3.png)![plot of chunk
figures-side](figure/figures-side-4.png)![plot of chunk
figures-side](figure/figures-side-5.png)

plot of chunk figures-side

##### Get plots of individual genes expression (example for the third module):

``` r
Seurat::DefaultAssay(seurat_object) <- "SCT"

i <- 3

Seurat::FeaturePlot(seurat_object, 
                    reduction = "pumap",
                    order = T,
                    features = m.gene.list[[i]],
                    ncol = 6,
                    raster = TRUE,
                    combine = TRUE)
```

![plot of chunk plot-genes](figure/plot-genes-1.png)

plot of chunk plot-genes

##### Get tables and plots with annotation of modules:

Functional annotation of obtained modules is done based on KEGG and
Reactome canonical metabolic pathways.

``` r
getAnnotationTables(network.annotation = network.annotation,
                    nets = results$nets,
                    work.dir = work.dir)
```

    # Pathway annotation for module 1 is produced

    # Pathway annotation for module 2 is produced

    # Pathway annotation for module 3 is produced

    # Pathway annotation for module 4 is produced

    # Pathway annotation for module 5 is produced

Example of the annotation table of the third module:

``` r
head(GAMclust:::read.tsv(file.path(work.dir, "m.3.pathways.tsv"))) |>
  kableExtra::kable() |>
  kableExtra::kable_styling(font_size = 8) |>
  kableExtra::column_spec(1, width = "1.6in") |>
  kableExtra::column_spec(2:6, width = "0.5in") |>
  kableExtra::column_spec(7, width = "1.2in")
```

| pathway | pval | padj | foldEnrichment | overlap | size | overlapGenes |
|:---|---:|---:|---:|---:|---:|:---|
| mmu_M00086: beta-Oxidation, acyl-CoA synthesis | 0.0145161 | 1 | 68.88889 | 1 | 1 | Acsl1 |
| mmu_M00089: Triacylglycerol biosynthesis | 0.0288447 | 1 | 34.44444 | 1 | 2 | Dgat1 |
| mmu00920: Sulfur metabolism | 0.0288447 | 1 | 34.44444 | 1 | 2 | Tst |
| mmu_M00415: Fatty acid elongation in endoplasmic reticulum | 0.0429877 | 1 | 22.96296 | 1 | 3 | Elovl7 |

Annotation heatmap for all modules:

``` r
getAnnotationHeatmap(work.dir = work.dir)
```

    # Processing module 1

    # Module size: 14

    # Processing module 2

    # Module size: 15

    # Processing module 3

    # Module size: 9

    # Processing module 4

    # Module size: 8

    # Processing module 5

    # Module size: 18

![plot of chunk plot-anno](figure/plot-anno-1.png)

plot of chunk plot-anno

##### Compare modules obtained in different runs:

You may also compare two results of running GAM-clustering on the same
dataset (e.g. runs with different parameters) or compare two results of
running GAM-clustering on different datasets (then set
`same.data=FALSE`).

``` r
modulesSimilarity(dir1 = work.dir,
                  dir2 = "old.dir",
                  name1 = "new",
                  name2 = "old",
                  same.data = TRUE,
                  use.genes.with.pos.score = TRUE,
                  work.dir = work.dir,
                  file.name = "comparison.png")
```
