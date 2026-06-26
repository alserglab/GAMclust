#' Function for batch solving SGMWCS problems
#'
#' @param mwcs_solver SGMWCS solver from mwcsr package
#' @export
seq_batch_solver <- function(mwcs_solver) {
  function(instances) { lapply(instances, mwcsr::solve_mwcsp, solver=mwcs_solver) }
}

#' Calculate correlation matrix from prepped expression matrices 
#' both a and b have the same number of columns, b could be a vector
#' data is precentered and prescaled
corFromPrep <- function(a, b) {
  if (is.vector(b)) {
    b <- t(b)
  }
  m <- a %*% t(b) # 32 x samples * samples x genes = 32 x genes
  # ensure that the similarities between genes & centroids in m are not biased by differences in the magnitudes of ...
  # ...gene expression values (normalizes the columns of m):
  m <- m / sqrt(max(rowSums(b**2))) # scales the rows of m by their Euclidean lengths -> m is [-1, 1]
  # ...centroid values (normalizes the rows of m):
  m <- sweep(x = m, MARGIN = 1, FUN = '/', STATS = sqrt(rowSums(a**2)))
  m
}

# calculate cosine similarity between rows of a and b (b could be a vector)
cosine <- function(a, b=NULL) {
  if (is.null(b)) {
    b <- a
  }
  if (is.vector(b)) {
    b <- t(b)
  }
  m <- a %*% t(b) 
  m <- sweep(x = m, MARGIN = 2, FUN = '/', STATS = sqrt(rowSums(b**2)))
  m <- sweep(x = m, MARGIN = 1, FUN = '/', STATS = sqrt(rowSums(a**2)))
  m
}

getCenter <- function(gene.exprs,
                      cluster.genes=seq_len(nrow(gene.exprs)),
                      cluster.genes.neg=c(),
                      method=c("pearson", "spearman")) {

  method <- match.arg(method)

  cluster.exprs <- rbind(
    gene.exprs[cluster.genes,, drop=F],
    -gene.exprs[cluster.genes.neg,, drop=F])

  if (method == "spearman") {
    cluster.exprs <- t(apply(cluster.exprs, 1, rank))
  }

  # cluster.exprs.znorm <- zScore(cluster.exprs)
  # center <- apply(cluster.exprs.znorm, 2, mean, na.rm=T)

  center <- apply(cluster.exprs, 2, mean, na.rm=T)
  # center <- center / sd(center, na.rm = T) # ???
  center[!is.finite(center)] <- NA
  center
}


updCenters <- function(cur.centers, m1, m2, E.prep, ms_mods) {
  
  genes <- unique(c(igraph::E(ms_mods[[m1]])[score > 0]$gene,
                    igraph::E(ms_mods[[m2]])[score > 0]$gene))
  
  if(length(genes) == 0){
    cur.centers[m1, ] <- colMeans(cur.centers[c(m1,m2), ])
  } else {
    cur.centers[m1, ] <- getCenter(E.prep, genes)
  }
  
  cur.centers <- cur.centers[-m2, , drop=F]
}


dualGraph <- function(m, what) {
  t1 <- data.table::as.data.table(igraph::as_data_frame(m)[, c("from", "to", what)])
  t2 <- data.table::copy(t1)

  data.table::setnames(t2, c("from", "to"), c("to", "from"))
  tt <- unique(rbind(t1, t2)[, c("from", what), with = FALSE])
  tt.g <- unique(merge(tt, tt, by="from", allow.cartesian = TRUE)
                 [, c(paste0(what, ".x"),
                      paste0(what, ".y")),
                   with = FALSE])
  dm <- igraph::graph_from_data_frame(d=tt.g, directed = FALSE)
  dm
}

doGeseca <- function(E.prep,
                     network.prep,
                     network.annotation,
                     modules,
                     scale = FALSE,
                     center = FALSE,
                     verbose = TRUE,
                     gesecaSeed = NULL) {
  
  E.prep_filtered <- E.prep[rownames(E.prep) %in% network.prep$gene, , drop = F]

  modules_path <- setNames(
    lapply(modules, function(cm) {igraph::E(cm)$gene}), 
    paste0("c.pos", seq_along(modules)))
  
  if (!is.null(gesecaSeed)) {
    gesecaRes <- rbindlist(lapply(seq_along(modules_path), 
                           function(i) {
                             if (exists(".Random.seed", envir = .GlobalEnv)) {
                               old_seed <- get(".Random.seed", envir = .GlobalEnv)
                               # Ensure it restores when the function finishes or errors out
                               on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv))
                             }
                             set.seed(gesecaSeed)
                             suppressWarnings(fgsea::geseca(pathways = modules_path[i],
                                           E = E.prep_filtered,
                                           scale = scale,
                                           center = center,
                                           sampleSize=1001))
                           }))
    gesecaRes <- gesecaRes[order(pval), ]
  } else {
    suppressWarnings(gesecaRes <- fgsea::geseca(pathways = modules_path,
                                                E = E.prep_filtered,
                                                scale = scale,
                                                center = center,
                                                sampleSize=1001)) # higher accuracy for better stability around threshold
    
  }
  
                              
  
  if(verbose){
    flog.info(">> geseca padjs were in range: %s", 
              paste(sprintf("%.2g", range(gesecaRes$padj)), collapse = " - "),
              name = "stats.logger")}
  
  gesecaRes
}

ulength <- function(x) length(unique(x))

zScore <- function(x) {
  x.means <- apply(x, 1, mean, na.rm=T)
  x.sds <- apply(x, 1, sd, na.rm=T)
  res <- sweep(sweep(x, 1, x.means), 1, x.sds, "/")
  return(res)
}

normalize.rows <- function(x) {
  x <- sweep(x, 1, apply(x, 1, min, na.rm=TRUE))
  sweep(x, 1, apply(x, 1, max, na.rm=TRUE), "/")
}

messagef <- function(...) {
  message(sprintf(...))
}

write.tsv <- function(table, dir, file=NULL, gzip=FALSE, row.names=NA, col.names=NA, ...) {
  name <- deparse(substitute(table))
  table <- as.data.frame(table)

  if (is.null(file)) {
    file <- file.path(dir, paste0(name, ".tsv", if (gzip) ".gz"))
  }

  if (is.na(row.names)) {
    row.names <- is.character(attr(table, "row.names"))
  }

  if (!row.names && is.na(col.names)) {
    col.names=T
  }

  for (c in colnames(table)) {
    if (is.character(table[[c]])) {
      table[[c]] <- sub("#", "", table[[c]])
    }
  }

  if (gzip) {
    file <- gzfile(file, "w")
  }
  write.table(table, file, quote=F,
              row.names=row.names, col.names=col.names, sep="\t")
  if (gzip) {
    close(file)
  }
}

read.tsv <- function(file, header=T, sep="\t", quote="", comment.char="", check.names=FALSE, ...) {
  args <- list(...)
  res <- read.table(file, header=header, sep=sep, quote=quote,
                    comment.char=comment.char, check.names=check.names,
                    stringsAsFactors=FALSE,
                    ...)
  if ((!"row.names" %in% names(args)) && (colnames(res)[1] == "")) {
    rownames(res) <- res[, 1]
    res[[1]] <- NULL
  }
  res
}

saveStats <- function(work.dir, rev, gesecaRes, iter.stats){
  
  dir.create(sprintf("%s/stats", work.dir), showWarnings=FALSE, recursive=TRUE)
  
  write.tsv(rev$centers.pos, file = sprintf("%s/stats/patterns_pos.tsv", work.dir))
  write.tsv(rev$centers.all, file = sprintf("%s/stats/patterns_all.tsv", work.dir))
  write.tsv(gesecaRes, file = sprintf("%s/stats/geseca_scores.tsv", work.dir))
  
  saveRDS(iter.stats, file = sprintf("%s/stats/iter.stats.rds", work.dir))
}

# Get names with number of genes in modules
getNamesWithLength <- function(dir, files, use.genes.with.pos.score.only=FALSE) {

  n_genes <- sapply(files, function(f) {

    len_m <- read.tsv(file.path(dir, f))
    if (use.genes.with.pos.score.only) {
      len_m <- len_m[which(len_m$score > 0), ]
    }

    len_m <- nrow(len_m)
    m <- sub(".genes.tsv", paste0("(", len_m, ")"), f)
  })

  return(n_genes)
}

# Number of intersected genes
compareModules <- function(dir1, dir2,
                           use.genes.with.pos.score.only=FALSE,
                           file.name="modules_similarity_intersection.png") {

  one.files <- list.files(dir1)
  one.files <- grep(one.files, pattern = "m.\\d+.genes.tsv", value = TRUE)
  two.files <- list.files(dir2)
  two.files <- grep(two.files, pattern = "m.\\d+.genes.tsv", value = TRUE)

  compare.df <- data.frame(matrix(nrow = length(one.files),
                                  ncol = length(two.files)))
  rownames(compare.df) <- getNamesWithLength(dir1, one.files, use.genes.with.pos.score.only)
  colnames(compare.df) <- getNamesWithLength(dir2, two.files, use.genes.with.pos.score.only)

  for (i in seq_along(one.files)) {
    # list of genes from dir1
    one.genes <- read.tsv(file.path(dir1, one.files[i]))
    if (use.genes.with.pos.score.only) {
      one.genes <- one.genes[which(one.genes$score > 0), ]
    }

    for (j in seq_along(two.files)) {
      # list of genes from dir2
      two.genes <- read.tsv(file.path(dir2, two.files[j]))
      if (use.genes.with.pos.score.only) {
        two.genes <- two.genes[which(two.genes$score > 0), ]
      }

      # intersection
      both.entrez <- length(intersect(one.genes$Entrez, two.genes$Entrez))
      compare.df[i, j] <- both.entrez
    }
  }

  # sort columns and rows
  compare.df <- compare.df[
    rownames(compare.df)[
      order(as.numeric(gsub(".*m\\.(\\d+).*", "\\1", rownames(compare.df))))],
    colnames(compare.df)[
      order(as.numeric(gsub(".*m\\.(\\d+).*", "\\1", colnames(compare.df))))],
    drop = F]

  if (any(dim(compare.df) == 1)) {
    compare.df.norm <- compare.df
  } else {
    compare.df.norm <- normalize.rows(compare.df) # 04.04.2022 Zhenya
    compare.df.norm[is.na(compare.df.norm)] <- 0  # 04.04.2022 Zhenya
  }

  all.same <- length(unique(unlist(as.vector(compare.df.norm)))) == 1
  if (all.same) {breaks <- c(0, 1)} else {breaks <- NA}
  # save heatmap
  pheatmap::pheatmap(
    compare.df.norm,  # 04.04.2022 Zhenya
    filename = file.name,
    cluster_rows=F, cluster_cols=F,
    angle_col = 315,
    display_numbers = compare.df, # 04.04.2022 Zhenya
    number_format = "%.f",
    fontsize_number = 16,
    fontsize_col = 18,
    fontsize_row = 18,
    width=8, height=5,
    breaks = breaks,
    show_rownames=T, show_colnames=T,
    main = paste0('Intersection',
                  ifelse(use.genes.with.pos.score.only, ": pos.genes", "")))
}

# Correlation between patterns (for same data)
compareModulesCor <- function(patterns.dir1, patterns.dir2,
                              file.name="modules_similarity_correlation.png") {

  one.patterns <- read.tsv(patterns.dir1)
  two.patterns <- read.tsv(patterns.dir2)
  two.patterns <- two.patterns[, colnames(one.patterns)]

  compare.df <- data.frame(matrix(nrow = nrow(one.patterns),
                                  ncol = nrow(two.patterns)))

  rownames(compare.df) <- paste("m.", seq(1:nrow(one.patterns)), sep="")
  colnames(compare.df) <- paste("m.", seq(1:nrow(two.patterns)), sep="")

  for (i in 1:nrow(one.patterns)) {
    # pattern from 1
    pattern1 <- as.numeric(one.patterns[i, ])

    for (j in 1:nrow(two.patterns)) {
      # pattern from 2
      pattern2 <- as.numeric(two.patterns[j,])

      patterns.cor <- cor(pattern1, pattern2, method = "spearman")

      # write fraction to the df
      compare.df[i, j] <- patterns.cor
    }
  }

  # sort columns and rows
  compare.df <- compare.df[
    rownames(compare.df)[
      order(as.numeric(gsub(".*m\\.(\\d+).*", "\\1", rownames(compare.df))))],
    colnames(compare.df)[
      order(as.numeric(gsub(".*m\\.(\\d+).*", "\\1", colnames(compare.df))))]]

  all.same <- length(unique(unlist(as.vector(compare.df)))) == 1
  if (all.same) {breaks <- c(0, 1)} else {breaks <- NA}
  # save heatmap
  pheatmap::pheatmap(
    compare.df,
    filename = file.name,
    cluster_rows=F, cluster_cols=F,
    angle_col = 315,
    display_numbers = T,
    fontsize_number = 16,
    fontsize_col = 18,
    fontsize_row = 18,
    width=8, height=5,
    breaks = breaks,
    show_rownames=T, show_colnames=T,
    main = paste0('Correlation of patterns',
                  ifelse(substr(patterns.dir1,
                                nchar(patterns.dir1) - 6,
                                nchar(patterns.dir1)) == "pos.tsv",
                         ": pos.genes", "")))
}

# Tversky index: Jaccard index / Sorensen-Dice coefficient
compareModulesIndex <- function(dir1, dir2,
                                jaccard=TRUE,
                                use.genes.with.pos.score.only=FALSE,
                                file.name="modules_similarity_index.png") {

  one.files <- list.files(dir1)
  one.files <- grep(one.files, pattern = "m.\\d+.genes.tsv", value = TRUE)
  two.files <- list.files(dir2)
  two.files <- grep(two.files, pattern = "m.\\d+.genes.tsv", value = TRUE)

  compare.df <- data.frame(matrix(nrow = length(one.files),
                                  ncol = length(two.files)))
  rownames(compare.df) <- getNamesWithLength(dir1, one.files, use.genes.with.pos.score.only)
  colnames(compare.df) <- getNamesWithLength(dir2, two.files, use.genes.with.pos.score.only)

  for (i in seq_along(one.files)) {
    # list of genes from dir1
    one.genes <- read.tsv(file.path(dir1, one.files[i]))
    if (use.genes.with.pos.score.only) {
      one.genes <- one.genes[which(one.genes$score > 0), ]
    }

    for (j in seq_along(two.files)) {
      # list of genes from dir2
      two.genes <- read.tsv(file.path(dir2, two.files[j]))
      if (use.genes.with.pos.score.only) {
        two.genes <- two.genes[which(two.genes$score > 0), ]
      }

      # intersection and set differences
      ab.intersection <- length(intersect(one.genes$Entrez, two.genes$Entrez))
      ab.diff <- length(setdiff(one.genes, two.genes))
      ba.diff <- length(setdiff(two.genes, one.genes))

      ab.coef <- ifelse(jaccard, 1, 0.5)

      # index
      index <- ab.intersection / (ab.intersection + ab.coef*ab.diff + ab.coef*ba.diff)
      # write fraction to the df
      compare.df[i, j] <- index
    }
  }

  # sort columns
  compare.df <- compare.df[rownames(compare.df)[
    order(as.numeric(gsub(".*m\\.(\\d+).*", "\\1", rownames(compare.df))))],
    colnames(compare.df)[
      order(as.numeric(gsub(".*m\\.(\\d+).*", "\\1", colnames(compare.df))))]]

  all.same <- length(unique(unlist(as.vector(compare.df)))) == 1
  if (all.same) {breaks <- c(0, 1)} else {breaks <- NA}
  # save heatmap
  pheatmap::pheatmap(
    compare.df,
    filename = file.name,
    cluster_rows=F, cluster_cols=F,
    angle_col = 315,
    display_numbers = T,
    fontsize_number = 16,
    fontsize_col = 18,
    fontsize_row = 18,
    width=8, height=5,
    breaks = breaks,
    show_rownames=T, show_colnames=T,
    # main = sprintf("Rows: %s / Columns: %s", rows.name, cols.name)
    main = paste0(
      ifelse(jaccard, 'Jaccard index', 'Sorensen-Dice index'),
      ifelse(use.genes.with.pos.score.only, ": pos.genes", "")))
}


# from rUtils package
collapseBy <- function(es, factor, FUN=median) {
  ranks <- apply(Biobase::exprs(es), 1, FUN)
  t <- data.frame(f=factor, i=seq_along(ranks), r=ranks)
  t <- t[order(t$r, decreasing=T), ]
  keep <- t[!duplicated(t$f) & !is.na(t$f),]$i
  res <- es[keep, ]
  Biobase::fData(res)$origin <- rownames(res)
  rownames(res) <- factor[keep]
  res
}

colors4heatmap <- c("snow2", 
                    "purple4", "#5E4FA2", "#486ecf", "#3288BD", "#66adc2", 
                    "#66C2A5", "#ABDDA4", "#E6F598", "#FEE08B", "#fad56b", 
                    "#f7bd16", "#FDAE61", "#fa912a", "#F46D43", "#D53E4F", 
                    "#9E0142", "#8B0001", "#830001", "#720000", "#690000")
