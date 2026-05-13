options(timeout = 600)
data("bulkRNAdata_ex_good")

met.to.filter <- data.table::fread(system.file("mets2mask.lst", package="GAMclust"))$ID
# KEGG:
network.kegg <- readRDS(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/network.kegg.rds"))
metabolites.annotation.kegg <- readRDS(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/met.kegg.db.rds"))

network.annotation.mm <- readRDS(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/org.Mm.eg.gatom.anno.rds"))

test_that("Pipeline works for bulk data", {

  work.dir <- file.path(tempdir(), "testGAMclustering")
  dir.create(work.dir, showWarnings = F)
  
  cplex.dir <- Sys.getenv("CPLEX_HOME")
  
  if (cplex.dir == "") {
    cplex.dir <- NULL # using MST heuristic, no CPLEX
  }
  
  solver <- mwcsr::virgo_solver(cplex_dir = cplex.dir)

  E.prep <- prepareData(E = Biobase::exprs(data.good),
                        gene.id.type = "Symbol",
                        keep.top.genes = 12000,
                        use.PCA = TRUE,
                        use.PCA.n = ncol(Biobase::exprs(data.good)) - 1,
                        repeats = seq_len(ncol(Biobase::exprs(data.good))),
                        network.annotation = network.annotation.mm)
  expect_true(ncol(E.prep) == ncol(Biobase::exprs(data.good)) - 1)
  
  network.prep <- prepareNetwork(E = E.prep,
                                 topology = "metabolites",
                                 met.to.filter = met.to.filter,
                                 network = network.kegg,
                                 network.annotation = network.annotation.mm)
  expect_equal(colnames(network.prep), expected = c("from", "to", "gene"))
  
  cur.centers <- preClustering(E.prep = E.prep,
                               network.prep = network.prep,
                               initial.number.of.clusters = 8,
                               network.annotation = network.annotation.mm)
  expect_true(is.matrix(cur.centers))
  expect_gt(sum(dim(cur.centers)), expected = 2)
  
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
  expect_length(results, 5)
  expect_true(all(
    names(results) %in% c("modules","nets",
                          "patterns.pos", "patterns.all",
                          "iter.stats")))
  
  r <- evaluate_promise(getGraphs(modules = results$modules,
                                  network.annotation = network.annotation.mm,
                                  metabolites.annotation = metabolites.annotation.kegg,
                                  seed.for.layout = 42,
                                  work.dir = work.dir))
  expect_true(is.null(r$result))
  expect_true(all(grepl("Graphs for module \\d are built", r$messages)))
  
  r <- evaluate_promise(getGeneTables(modules = results$modules,
                                      nets = results$nets,
                                      patterns = results$patterns.pos,
                                      gene.exprs = E.prep,
                                      network.annotation = network.annotation.mm,
                                      work.dir = work.dir))
  expect_length(r$result, length(results$modules))
  
  r <- evaluate_promise(getAnnotationTables(network.annotation = network.annotation.mm,
                                            nets = results$nets,
                                            work.dir = work.dir))
  expect_true(is.null(r$result))
  expect_true(all(grepl("Pathway annotation for module \\d is produced", r$messages)))
  
  r <- evaluate_promise(getAnnotationHeatmap(work.dir = work.dir))
  expect_true(class(r$result) == "pheatmap")

  unlink(work.dir, recursive = T)

})
