# GAM-clustering Tutorial Notebook For Spatial RNA-Seq Data

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
gene2reaction.extra <- data.table::fread("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/gene2reaction.combined.hsa.eg.tsv", colClasses="character")
```

2.  load species-specific network annotation:
    [`org.Hs.eg.gatom.anno`](http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/org.Hs.eg.gatom.anno.rds)
    for human data or
    [`org.Mm.eg.gatom.anno`](http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/org.Mm.eg.gatom.anno.rds)
    for mouse data;

``` r
network.annotation <- readRDS(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/org.Hs.eg.gatom.anno.rds"))
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
cplex.dir <-  "/opt/ibm/ILOG/CPLEX_Studio1271"
```

``` r
solver <- mwcsr::virgo_solver(cplex_dir = cplex.dir)
```

5.  set working directory where the results will be saved to.

``` r
work.dir <- "results_spat"
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

This vignette shows how to process spatial RNA-seq data.

Let’s load the data and take 10,000 genes for the GAM-clustering
analysis.

``` r
seurat_object <- readRDS(url("http://artyomovlab.wustl.edu/publications/supp_materials/GAMclust/275_T_seurat.rds"))

seurat_object <- Seurat::SCTransform(seurat_object, 
                                     assay = "Spatial", 
                                     variable.features.n = 10000,
                                     verbose = FALSE)

E <- as.matrix(Seurat::GetAssayData(object = seurat_object,
                                    assay = "SCT",
                                    slot = "scale.data"))

nrow(E) # ! make sure this value is in range from 6,000 to 12,000
```

    # [1] 10000

``` r
E[1:3, 1:3]
```

    #       AAACAAGTATCTCCCA-1 AAACACCAATAACTGC-1 AAACAGAGCGACTCCT-1
    # NOC2L           1.249130          0.8296517          3.7173470
    # ISG15           2.642288         -0.7594430          0.5659099
    # AGRN           -0.620876         -0.7318498          0.5300554

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

    #                 PC1       PC2        PC3
    # 4504       6.326338  2.513886  0.5230565
    # 22865     -1.635009 -2.005070  0.9369077
    # 100129792 -1.437753 -1.059253 -1.0455837

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

    # INFO [2026-06-15 22:20:39] Global metabolite network contains 5746 edges.
    # INFO [2026-06-15 22:20:39] Largest connected component of this global network contains 1276 nodes and 4403 edges.

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

    # INFO [2026-06-15 22:20:39] 1002 metabolic genes from the analysed dataset mapped to this component.

``` r
cur.centers[1:3, 1:3]
```

    #           PC1        PC2        PC3
    # 1  0.08031134 -0.3373824  2.2391351
    # 2  3.50078547 -0.7895749  3.3735482
    # 3 -2.35233343 -1.4231047 -0.8383281

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

                         cor.threshold = 0.5,
                         p.adj.val.threshold = 1e-5,

                         batch.solver = seq_batch_solver(solver),
                         work.dir = work.dir,
                         
                         show.intermediate.clustering = FALSE,
                         verbose = FALSE,
                         collect.stats = TRUE)
```

    # INFO [2026-06-15 22:20:39] GAM-CLUSTERING starts here.
    # INFO [2026-06-15 22:20:39] [*] Iteration 1
    # INFO [2026-06-15 22:21:51] [*] Iteration 2
    # INFO [2026-06-15 22:23:03] >> max cor exceeded 0.5: 0.76
    # INFO [2026-06-15 22:23:03] [*] Iteration 3
    # INFO [2026-06-15 22:24:12] [*] Iteration 4
    # INFO [2026-06-15 22:25:21] >> max cor exceeded 0.5: 0.73
    # INFO [2026-06-15 22:25:21] [*] Iteration 5
    # INFO [2026-06-15 22:26:27] [*] Iteration 6
    # INFO [2026-06-15 22:27:32] >> max cor exceeded 0.5: 0.7
    # INFO [2026-06-15 22:27:32] [*] Iteration 7
    # INFO [2026-06-15 22:28:37] [*] Iteration 8
    # INFO [2026-06-15 22:29:43] >> max cor exceeded 0.5: 0.66
    # INFO [2026-06-15 22:29:43] [*] Iteration 9
    # INFO [2026-06-15 22:30:46] [*] Iteration 10
    # INFO [2026-06-15 22:31:49] >> max cor exceeded 0.5: 0.65
    # INFO [2026-06-15 22:31:49] [*] Iteration 11
    # INFO [2026-06-15 22:32:49] >> max cor exceeded 0.5: 0.64
    # INFO [2026-06-15 22:32:49] [*] Iteration 12
    # INFO [2026-06-15 22:33:47] [*] Iteration 13
    # INFO [2026-06-15 22:34:46] >> max cor exceeded 0.5: 0.6
    # INFO [2026-06-15 22:34:46] [*] Iteration 14
    # INFO [2026-06-15 22:35:43] [*] Iteration 15
    # INFO [2026-06-15 22:36:41] >> max cor exceeded 0.5: 0.59
    # INFO [2026-06-15 22:36:41] [*] Iteration 16
    # INFO [2026-06-15 22:37:37] [*] Iteration 17
    # INFO [2026-06-15 22:38:32] >> max cor exceeded 0.5: 0.57
    # INFO [2026-06-15 22:38:32] [*] Iteration 18
    # INFO [2026-06-15 22:39:26] [*] Iteration 19
    # INFO [2026-06-15 22:40:19] >> max cor exceeded 0.5: 0.54
    # INFO [2026-06-15 22:40:19] [*] Iteration 20
    # INFO [2026-06-15 22:41:11] [*] Iteration 21
    # INFO [2026-06-15 22:42:02] >> max cor exceeded 0.5: 0.54
    # INFO [2026-06-15 22:42:02] [*] Iteration 22
    # INFO [2026-06-15 22:42:51] [*] Iteration 23
    # INFO [2026-06-15 22:43:38] >> max cor exceeded 0.5: 0.53
    # INFO [2026-06-15 22:43:39] [*] Iteration 24
    # INFO [2026-06-15 22:44:24] [*] Iteration 25
    # INFO [2026-06-15 22:45:10] >> max cor exceeded 0.5: 0.51
    # INFO [2026-06-15 22:45:10] [*] Iteration 26
    # INFO [2026-06-15 22:45:54] [*] Iteration 27
    # INFO [2026-06-15 22:46:38] >> max cor exceeded 0.5: 0.5
    # INFO [2026-06-15 22:46:38] [*] Iteration 28
    # INFO [2026-06-15 22:47:18] [*] Iteration 29
    # INFO [2026-06-15 22:47:59] >> max cor exceeded 0.5: 0.5
    # INFO [2026-06-15 22:47:59] [*] Iteration 30
    # INFO [2026-06-15 22:48:37] [*] Iteration 31
    # INFO [2026-06-15 22:49:19] [*] Iteration 32
    # INFO [2026-06-15 22:49:58] [*] Iteration 33
    # INFO [2026-06-15 22:50:35] [*] Iteration 34
    # INFO [2026-06-15 22:51:12] [*] Iteration 35
    # INFO [2026-06-15 22:51:45] [*] Iteration 36
    # INFO [2026-06-15 22:52:20] [*] Iteration 37
    # INFO [2026-06-15 22:52:51] [*] Iteration 38
    # INFO [2026-06-15 22:53:24] [*] Iteration 39
    # INFO [2026-06-15 22:53:53] [*] Iteration 40
    # INFO [2026-06-15 22:54:24] [*] Iteration 41
    # INFO [2026-06-15 22:54:51] [*] Iteration 42
    # INFO [2026-06-15 22:55:21] [*] Iteration 43
    # INFO [2026-06-15 22:55:47] [*] Iteration 44
    # INFO [2026-06-15 22:56:11] [*] Iteration 45
    # INFO [2026-06-15 22:56:35] [*] Iteration 46
    # INFO [2026-06-15 22:56:55] [*] Iteration 47
    # INFO [2026-06-15 22:57:17] [*] Iteration 48
    # INFO [2026-06-15 22:57:35] [*] Iteration 49
    # INFO [2026-06-15 22:57:54] [*] Iteration 50
    # INFO [2026-06-15 22:58:11] [*] Iteration 51
    # INFO [2026-06-15 22:58:29] [*] Iteration 52
    # INFO [2026-06-15 22:58:43] [*] Iteration 53
    # INFO [2026-06-15 22:58:59] [*] Iteration 54
    # INFO [2026-06-15 22:59:11] [*] Iteration 55
    # INFO [2026-06-15 22:59:27] GAM-CLUSTERING ends here.

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

Example of the graph of the fourth module:

![plot of chunk unnamed-chunk-7](results_spat/m.4.svg)

plot of chunk unnamed-chunk-7

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

Example of the gene table of the fourth module:

``` r
head(GAMclust:::read.tsv(file.path(work.dir, "m.4.genes.tsv"))) |>
  kableExtra::kable() |>
  kableExtra::kable_styling()
```

| symbol | Entrez |    score |       cor |
|:-------|-------:|---------:|----------:|
| SLC2A3 |   6515 | 3.313484 | 0.9539268 |
| SLC2A1 |   6513 | 2.705281 | 0.9306558 |
| EGLN3  | 112399 | 2.619330 | 0.9266640 |
| PGK1   |   5230 | 2.445532 | 0.9200964 |
| GLUL   |   2752 | 2.456103 | 0.9196962 |
| GAPDH  |   2597 | 2.481569 | 0.9178584 |

##### Get plots of patterns:

``` r
for(i in 1:length(m.gene.list)){
 
  print(fgsea::plotCoregulationProfileSpatial(m.gene.list[[i]], 
                                             seurat_object,
                                             title = paste0("module ", i)))
}
```

![plot of chunk figures-side](figure/figures-side-1.png)![plot of chunk
figures-side](figure/figures-side-2.png)![plot of chunk
figures-side](figure/figures-side-3.png)![plot of chunk
figures-side](figure/figures-side-4.png)

plot of chunk figures-side

##### Get plots of individual genes expression (example for the fourth module):

``` r
Seurat::DefaultAssay(seurat_object) <- "SCT"

i <- 4

Seurat::SpatialFeaturePlot(seurat_object, 
                           features = m.gene.list[[i]], 
                           pt.size.factor = 2, stroke = 0.01, alpha = 0.7,
                           ncol = 8)
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

Example of the annotation table of the fourth module:

``` r
head(GAMclust:::read.tsv(file.path(work.dir, "m.4.pathways.tsv"))) |>
  kableExtra::kable() |>
  kableExtra::kable_styling(font_size = 8) |>
  kableExtra::column_spec(1, width = "1.6in") |>
  kableExtra::column_spec(2:6, width = "0.5in") |>
  kableExtra::column_spec(7, width = "1.2in")
```

| pathway | pval | padj | foldEnrichment | overlap | size | overlapGenes |
|:---|---:|---:|---:|---:|---:|:---|
| R-HSA-70171: Glycolysis | 0.0000000 | 0.0000000 | 17.032258 | 12 | 13 | GNPDA1 ENO2 ALDOA GAPDH GPI HK2 PFKFB4 PFKP PGK1 PKM BPGM TPI1 |
| hsa00500: Starch and sucrose metabolism | 0.0000703 | 0.0041371 | 9.225806 | 5 | 10 | GPI GYS1 HK2 PGM1 PYGL |
| hsa_M00570: Isoleucine biosynthesis, threonine =\> 2-oxobutanoate =\> isoleucine | 0.0028474 | 0.0894088 | 18.451613 | 2 | 2 | SDSL BCAT2 |
| R-HSA-196836: Vitamin C (ascorbate) metabolism | 0.0082525 | 0.1850920 | 12.301075 | 2 | 3 | SLC2A1 SLC2A3 |
| R-HSA-71240: Tryptophan catabolism | 0.0159469 | 0.3129585 | 9.225806 | 2 | 4 | KYAT3 SLC3A2 |
| hsa00220: Arginine biosynthesis | 0.0256824 | 0.4222898 | 7.380645 | 2 | 5 | GLUL GPT2 |

Annotation heatmap for all modules:

``` r
getAnnotationHeatmap(work.dir = work.dir)
```

    # Processing module 1

    # Module size: 50

    # Processing module 2

    # Module size: 50

    # Processing module 3

    # Module size: 36

    # Processing module 4

    # Module size: 31

![plot of chunk plot-anno](figure/plot-anno-1.png)

plot of chunk plot-anno

##### Compare modules obtained in different runs:

You may also compare two results of running GAM-clustering on the same
dataset (e.g. runs with different parameters) or compare two results of
running GAM-clustering on different datasets (then set
`same.data=FALSE`).

``` r
modulesSimilarity(dir1 = work.dir,
                  dir2 = "old_dir",
                  name1 = "new",
                  name2 = "old",
                  same.data = TRUE,
                  use.genes.with.pos.score = TRUE,
                  work.dir = work.dir,
                  file.name = "comparison.png")
```
