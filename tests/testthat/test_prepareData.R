options(timeout = 600)
data("bulkRNAdata_ex_good")
data("bulkRNAdata_ex_bad")
data("scRNAdata_small_ex")

network.annotation.mm <- readRDS(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/org.Mm.eg.gatom.anno.rds"))
network.annotation.hs <- readRDS(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/org.Hs.eg.gatom.anno.rds"))


test_that("prepareData works for Mus musculus with PCA", {

  E.prep <- prepareData(E = Biobase::exprs(data.good),
                        gene.id.type = "Symbol",
                        keep.top.genes = 12000,
                        use.PCA = TRUE,
                        use.PCA.n = ncol(Biobase::exprs(data.good)) - 1,
                        repeats = seq_len(ncol(Biobase::exprs(data.good))),
                        network.annotation = network.annotation.mm)

  expect_true(ncol(E.prep) == ncol(Biobase::exprs(data.good)) - 1)
})


test_that("prepareData works for Mus musculus with expression", {
  
  E.prep <- prepareData(E = Biobase::exprs(data.good),
                        gene.id.type = "Symbol",
                        keep.top.genes = 12000,
                        use.PCA = FALSE,
                        repeats = seq_len(ncol(Biobase::exprs(data.good))),
                        network.annotation = network.annotation.mm)
  
  expect_true(ncol(E.prep) == ncol(Biobase::exprs(data.good)))
  expect_equal(colnames(E.prep), expected = colnames(Biobase::exprs(data.good)))
})


test_that("prepareData works for Homo sapiens with PCA", {
  
  E.prep <- prepareData(E = Biobase::exprs(sc.data),
                        gene.id.type = "Entrez",
                        keep.top.genes = 12000,
                        use.PCA = TRUE,
                        use.PCA.n = ncol(Biobase::exprs(sc.data)) - 1,
                        repeats = seq_len(ncol(Biobase::exprs(sc.data))),
                        network.annotation = network.annotation.hs)
  
  expect_true(ncol(E.prep) == ncol(Biobase::exprs(sc.data)) - 1)
})


test_that("prepareData stops with undefined gene annotation", {

  expect_error(
    prepareData(E = Biobase::exprs(data.good),
                gene.id.type = "Something",
                keep.top.genes = 12000,
                use.PCA = TRUE,
                use.PCA.n = ncol(Biobase::exprs(data.good)) - 1,
                repeats = seq_len(ncol(Biobase::exprs(data.good))),
                network.annotation = network.annotation.mm), 
    "Something")

})


test_that("prepareData no annotation when gene.id.type is Entrez", {
  
  expect_message(prepareData(E = Biobase::exprs(data.good),
                           gene.id.type = "Entrez",
                           keep.top.genes = 12000,
                           use.PCA = TRUE,
                           use.PCA.n = ncol(Biobase::exprs(data.good)) - 1,
                           repeats = seq_len(ncol(Biobase::exprs(data.good))),
                           network.annotation = network.annotation.mm),
               "No gene annotation was performed")
  
})


test_that("prepareData no annotation when gene.id.type is NULL and works", {
  
  expect_message(prepareData(E = Biobase::exprs(data.good),
                             gene.id.type = NULL,
                             keep.top.genes = 12000,
                             use.PCA = TRUE,
                             use.PCA.n = ncol(Biobase::exprs(data.good)) - 1,
                             repeats = seq_len(ncol(Biobase::exprs(data.good))),
                             network.annotation = network.annotation.mm),
                 "No gene annotation was performed")
  
})

test_that("prepareData works with duplicate gene names", {
  E <- Biobase::exprs(data.good)
  rownames(E)[1] <- "Smim36" # has two Entrez ids
  E.prep <- prepareData(E = E,
                        gene.id.type = "Symbol",
                        keep.top.genes = 12000,
                        use.PCA = FALSE,
                        repeats = seq_len(ncol(Biobase::exprs(data.good))),
                        network.annotation = network.annotation.mm)

  new2old <- attributes(E.prep)$original.gene.names 
  newName <- names(which(new2old == rownames(E)[1]))
  
  expect_gt(cor(E.prep[newName, ], E[1,]), 0.99) 
})
