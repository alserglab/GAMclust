options(timeout = 600)
data("bulkRNAdata_ex_good_prep")

network.annotation.mm <- readRDS(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/org.Mm.eg.gatom.anno.rds"))


test_that("preClustering works", {

  cur.centers <- preClustering(E.prep = E.prep,
                               network.prep = network.prep,
                               initial.number.of.clusters = 8,
                               network.annotation = network.annotation.mm)

  expect_true(is.matrix(cur.centers))
  expect_gt(sum(dim(cur.centers)), expected = 2)
  expect_message(preClustering(E.prep = E.prep,
                               network.prep = network.prep,
                               initial.number.of.clusters = 8,
                               network.annotation = network.annotation.mm),
                 "metabolic genes from the analysed dataset mapped to this component")

})


# test_that("preClustering works with ICA", {
# 
#   cur.centers <- preClustering(E.prep = E.prep,
#                                use.ICA = T,
#                                network.prep = network.prep,
#                                initial.number.of.clusters = 8,
#                                network.annotation = network.annotation.mm)
# 
#   expect_true(is.matrix(cur.centers))
#   expect_gt(sum(dim(cur.centers)), expected = 2)
#   expect_message(preClustering(E.prep = E.prep,
#                                use.ICA = T,
#                                network.prep = network.prep,
#                                initial.number.of.clusters = 8,
#                                network.annotation = network.annotation.mm),
#                  "metabolic genes from the analysed dataset mapped to this component")
# 
# })
# 
# 
# test_that("preClustering with ICA stops if no PCs", {
#   
#   expect_error(preClustering(E.prep = E.prep.expr,
#                              use.ICA = T,
#                              network.prep = network.prep,
#                              initial.number.of.clusters = 8,
#                              network.annotation = network.annotation.mm),
#                "To perform ICA, set ")
#   
# })



test_that("gamClustering works with Virgo solver", {

  cplex.dir <- Sys.getenv("CPLEX_HOME")
  
  if (cplex.dir == "") {
    cplex.dir <- NULL # using MST heuristic, no CPLEX
  }
  
  solver <- mwcsr::virgo_solver(cplex_dir = cplex.dir)

  cur.centers <- preClustering(E.prep = E.prep,
                               network.prep = network.prep,
                               initial.number.of.clusters = 8,
                               network.annotation = network.annotation.mm)

  work.dir <- file.path(tempdir(), "testGAMclustering")
  dir.create(work.dir)

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
   expect_true(all(sapply(results$modules, class) == "igraph"))
   expect_false(any(is.na(results$patterns.pos)))

   unlink(work.dir, recursive = T)

})
