#' Prepare gene expression matrix
#'
#' Prepares a gene expression matrix for GAM-clustering analysis. The function
#' optionally converts gene identifiers to the annotation format used by the
#' metabolic network, removes duplicates, retains the most highly expressed
#' genes, and standardizes expression values.
#' Biological replicates can be collapsed by averaging, and dimensionality
#' reduction can be applied using PCA.
#' The resulting matrix contains only processed features suitable for
#' downstream network-based module discovery.
#'
#' @param E Expression matrix with rownames as gene symbols.
#' @param gene.id.type Gene ID type.
#' @param keep.top.genes Which top of the most expressed genes to keep for the further analysis.
#' @param use.PCA Whether to reduce matrix dimensionality by PCA or not.
#' @param repeats Here you may collapse biological replicas by providing vector with repeated sample names
#' @param network.annotation Metabolic network annotation.
#' @return Expression matrix prepared for the analysis.
#' @import data.table futile.logger
#' @export
prepareData <- function(
    E,
    gene.id.type = NULL,
    keep.top.genes = 12000,
    use.PCA = TRUE,
    use.PCA.n = 50,
    repeats = seq_len(ncol(E)),
    network.annotation){
  
  if(any(duplicated(repeats))){
    colnames(E) <- repeats
    E <- t(apply(E, 1, function(x) tapply(x, colnames(E), mean)))
  }
  
  # memory efficient, limits creating temp matrices in memory, except for t(E)
  mu <- matrixStats::rowMeans2(E)
  s  <- matrixStats::rowSds(E);   s[s == 0] <- 1
  
  if (use.PCA && use.PCA.n >= ncol(E) - 1)  {
    use.PCA <- FALSE
  }
  
  if (use.PCA) {
    pcaRev  <- irlba::irlba(t(E), nv = use.PCA.n,
                            center = mu, scale = s)
    E.red <- pcaRev$v %*% diag(pcaRev$d)
    rownames(E.red) <- rownames(E)
    colnames(E.red) <- paste0("PC", seq_len(use.PCA.n))
    E <- E.red
    # rescale to ignore everything beyond the top components
    E <- t(base::scale(t(E), center = FALSE, scale = TRUE))
    E[is.na(E)] <- 0 # for sd=0
  } else {
    E <- t(base::scale(t(E), center = mu, scale = s))  
  }

  new2old <- rownames(E)
  
  if(is.null(gene.id.type) || gene.id.type == network.annotation$baseId){
    flog.info("No gene annotation was performed.", name = "stats.logger")
  } else {
    
    if(gene.id.type %in% names(network.annotation$mapFrom)){
      
      rownames.dubl <- network.annotation$mapFrom[[gene.id.type]][rownames(E)]
      rownames(E) <- rownames.dubl$gene
      
    } else {
      valid.types.string <- paste(c(names(network.annotation$mapFrom), "Entrez"), collapse = ", ")
      message.string <- sprintf("Invalid `gene.id.type`: %s. Must be one of: %s.", gene.id.type, valid.types.string)
      flog.error(message.string, name = "stats.logger")
      stop(message.string)
    }
  }
  names(new2old) <- rownames(E)
  
  
  E <- E[order(mu, decreasing = T), ]
  E <- E[!is.na(rownames(E)), ]
  E <- E[!duplicated(rownames(E)), ]
  E <- head(E, keep.top.genes)

  attributes(E)$original.gene.names <- new2old
  
  E
}

#' Prepare network
#'
#' Constructs a metabolic interaction network tailored to the analyzed dataset.
#' The function integrates network topology with gene–enzyme–reaction
#' annotations, removes filtered metabolites, and retains only genes present
#' in the processed expression matrix.
#' Networks can be represented at either the metabolite or atom level.
#' The final output corresponds to the largest connected component of the
#' metabolic network and serves as the graph structure for GAM-clustering.
#'
#' @param E Expression matrix after the `prepareData()` function.
#' @param network Metabolic network.
#' @param topology Vertices can be represented either as `metabolites`, either as `atoms`.
#' @param met.to.filter Metabolites that should not be used as connections in the module.
#' @param network.annotation Metabolic network annotation.
#' @param gene2reaction.extra For a combined network: supplementary file with genes that either come from proteome or are not linked to a specific enzyme.
#' @return Edges of the final network.
#' @import data.table
#' @export
prepareNetwork <- function(
    E,
    network,
    topology = c("metabolites", "atoms"),
    met.to.filter = data.table::fread(system.file("mets2mask.lst", package="GAMclust"))$ID,
    network.annotation,
    gene2reaction.extra = NULL){
  
  topology <- match.arg(topology)

  globalEdgeTable_pre <- as.data.frame(network$reaction2align)
  globalEdgeTable_pre <- merge(globalEdgeTable_pre, network$enzyme2reaction)
  globalEdgeTable_pre <- merge(globalEdgeTable_pre, network.annotation$gene2enzyme)
  
  if(!is.null(gene2reaction.extra)){
    globalEdgeTable_pre.extra <- merge(gene2reaction.extra, network$reaction2align)
    globalEdgeTable_pre.extra <- cbind(enzyme = "-.-.-.-", globalEdgeTable_pre.extra)
    globalEdgeTable_pre <- rbind(globalEdgeTable_pre, 
                                 globalEdgeTable_pre.extra) 
    }
  colnames(globalEdgeTable_pre)[which(colnames(globalEdgeTable_pre) == "atom.x")] <- "from"
  colnames(globalEdgeTable_pre)[which(colnames(globalEdgeTable_pre) == "atom.y")] <- "to"
  globalEdgeTable_pre$from.m <- network$atoms$metabolite[match(globalEdgeTable_pre$from, network$atoms$atom)]
  globalEdgeTable_pre$to.m <- network$atoms$metabolite[match(globalEdgeTable_pre$to, network$atoms$atom)]
  globalEdgeTable_pre <- globalEdgeTable_pre[which(!globalEdgeTable_pre$from.m %in% met.to.filter), ]
  globalEdgeTable_pre <- globalEdgeTable_pre[which(!globalEdgeTable_pre$to.m %in% met.to.filter), ]
  globalEdgeTable_pre <- globalEdgeTable_pre[which(globalEdgeTable_pre$gene %in% rownames(E)), ]
  
  if(topology == "atoms"){
    
    globalEdgeTable_pre <- globalEdgeTable_pre[, c("from", "to", "gene")]
    globalEdgeTable_pre <- globalEdgeTable_pre[!duplicated(globalEdgeTable_pre), ]
    
    flog.info("Global atom network contains %s edges.", dim(globalEdgeTable_pre)[1], 
              name = "stats.logger")
  }
  if(topology == "metabolites"){
    
    globalEdgeTable_pre <- globalEdgeTable_pre[, c("from.m", "to.m", "gene")]
    globalEdgeTable_pre <- globalEdgeTable_pre[!duplicated(globalEdgeTable_pre), ]
    colnames(globalEdgeTable_pre)[which(colnames(globalEdgeTable_pre) == "from.m")] <- "from"
    colnames(globalEdgeTable_pre)[which(colnames(globalEdgeTable_pre) == "to.m")] <- "to"
    
    flog.info("Global metabolite network contains %s edges.", dim(globalEdgeTable_pre)[1], 
              name = "stats.logger")
  }

  if (dim(globalEdgeTable_pre)[1] == 0) {
    message.string <- "No metabolic genes from the analysed dataset mapped to the metabolic network.\n
      In this case GAM-clustering will not work. Please try another subset of genes if it is possible."
    flog.error(message.string, name = "stats.logger")
    stop(message.string)
  }
  
  globalEdgeTable_pre_graph <- igraph::graph_from_data_frame(globalEdgeTable_pre, directed=FALSE)
  globalEdgeTable_pre_graph_cc <- igraph::decompose.graph(globalEdgeTable_pre_graph)
  globalGraph <- globalEdgeTable_pre_graph_cc[[which.max(sapply(globalEdgeTable_pre_graph_cc, igraph::vcount))]]
  # multi-edges, loops
  
  flog.info("Largest connected component of this global network contains %s nodes and %s edges.", 
            igraph::vcount(globalGraph), igraph::ecount(globalGraph), 
            name = "stats.logger")

  x.1p <- paste(globalEdgeTable_pre$from, globalEdgeTable_pre$to)
  x.2p <- with(igraph::as_data_frame(globalGraph), paste(c(from, to), c(to, from)))
  globalEdgeTable <- globalEdgeTable_pre[x.1p %in% x.2p, ]
  
  globalEdgeTable
}

#' Defining initial patterns
#'
#' Generates an initial set of transcriptional patterns that serve as
#' starting points for module detection. The function first restricts the
#' expression matrix to metabolic genes represented in the network and then
#' performs clustering to identify preliminary expression signatures.
#' By default, cluster centers are obtained using k-means clustering.
#' These initial patterns are subsequently refined by the GAM-clustering algorithm.
#'
#' @param E.prep Expression matrix after the `prepareData()` function.
#' @param network.prep Network edge table driven from `prepareNetwork()` function.
#' @param initial.number.of.clusters The number of clusters for the initial approximation of modules.
#' @param network.annotation Metabolic network annotation.
#' @return Initial patterns.
#' @export
preClustering <- function(E.prep,
                          network.prep,
                          initial.number.of.clusters = 32,
                          network.annotation,
                          use.ICA = FALSE
                          ){
  
  E.prep <- E.prep[rownames(E.prep) %in% network.prep$gene, , drop = F]
  flog.info("%d metabolic genes from the analysed dataset mapped to this component.",
            dim(E.prep)[1],
            name = "stats.logger")

  ### gene.cor <- cor(t(E.prep), use="pairwise.complete.obs")
  # gene.cor <- (E.prep %*% t(E.prep)) / max(rowSums(E.prep**2)) # max(rowSums(E.prep**2)) = x, while x+1 samples
  # gene.cor.dist <- as.dist(1 - gene.cor)
  # gene.pam <- cluster::pam(gene.cor.dist, k=initial.number.of.clusters)
  # cur.centers <- E.prep[gene.pam$medoids,]
  # OR
  gene.kmeans <- kmeans(E.prep, centers=initial.number.of.clusters)
  cur.centers <- gene.kmeans$centers
  
  if(use.ICA == TRUE){
    if(all(grepl("PC", colnames(E.prep)))) {
      ica_result <- fastICA::fastICA(t(E.prep), n.comp = initial.number.of.clusters) 
      cur.centers <- t(ica_result$S) } else {
        message.string <- "To perform ICA, set `use.PCA = TRUE` in `prepareData()` function."
        flog.error(message.string, name = "stats.logger")
        stop(message.string)
      }
  }
  
  cur.centers
}

#' GAM-clustering analysis
#'
#' Performs the core GAM-clustering procedure to identify transcriptionally
#' coordinated metabolic modules. Starting from initial expression patterns,
#' the algorithm iteratively scores genes, solves SGMWCS optimization problems
#' on the metabolic network, and updates module-specific expression profiles
#' until convergence.
#' Additional refinement steps merge highly similar modules, remove
#' uninformative patterns, and control module size. The function returns
#' final metabolic modules, scored networks, module expression patterns,
#' and detailed iteration statistics.
#' 
#' @section Parameter tuning:
#' There is a set of parameters which determine the size and number of the
#' final modules. We recommend starting with the default settings, but the
#' following adjustments may be useful:
#'
#' \itemize{
#'   \item If you consider final modules to be too small or too big and it 
#'   complicates interpretation for you, you can either increase or reduce by 
#'   10 units the max.module.size parameter.
#'   \item If among final modules you consider presence of any modules with too 
#'   similar patterns, you can reduce by 0.1 units the cor.threshold parameter.
#'   \item If among final modules you consider presence of any uninformative 
#'   modules, you can reduce by 10 times the p.adj.val.threshold parameter.
#' }
#' 
#' @param E.prep Expression matrix after the `prepareData()` function.
#' @param E.prep Expression matrix after the `prepareData()` function.
#' @param network.prep Network edge table driven from `prepareNetwork()` function.
#' @param cur.centers Initial patterns produced by `preClustering()` function.
#' @param start.base The parameter which influences modules sizes.
#' @param base.dec The value controlling how strongly `base` parameter should be reduced if some module's size is bigger that `max.module.size`.
#'                 The update rule is: `base <- base * (1 - base.dec)`. Detaulf: `0.1`.
#' @param max.module.size Maximal number of unique genes in the final module.
#' @param cor.threshold Threshold for correlation between module patterns.
#' @param p.adj.val.threshold Padj threshold of geseca score for final modules.
#' @param batch.solver Solver for SGMWCS problem.
#' @param work.dir Working directory where results should be saved.
#' @param show.intermediate.clustering Whether to show or not heatmap of intermediate clusters.
#' @param verbose Verbose running.
#' @param collect.stats Whether to save or not running statistics.
#' @param reference.patterns Matrix of reference patterns to track correlation of centers against. 
#'     Pattern per row. Number of columns should be the sames as in E.prep.
#' @return results$modules -- Metabolic modules.
#' @return results$nets -- Scored networks.
#' @return results$patterns.pos -- Modules' patterns (genes with positive score only considered).
#' @return results$patterns.all -- Modules' patterns (all genes considered).
#' @return results$iter.stats -- Statistics from iterations.
#' @export
gamClustering <- function(E.prep,
                          network.prep,
                          cur.centers,
                          
                          start.base = 0.5,
                          base.dec = 0.1,
                          max.module.size = 50,
                          
                          cor.threshold = 0.8,
                          p.adj.val.threshold = 1e-5,
                          
                          batch.solver = seq_batch_solver(solver),
                          work.dir,
                          
                          show.intermediate.clustering = TRUE,
                          verbose = TRUE,
                          collect.stats = TRUE,
                          reference.patterns = NULL
                          ){
  
  flog.info("GAM-CLUSTERING starts here.", name = "stats.logger")
  
  iteration <- 1
  base <- start.base
  iter.stats <- list()
  
  while (T) {
    
    k <- 1
    revs <- list()
    
    while (T) {
      
      flog.info("[*] Iteration %s", iteration, name = "stats.logger")
      
      if (!is.null(reference.patterns)) {
        flog.info("correlations with references: %s", 
                  paste0(matrixStats::colMaxs(cosine(cur.centers, 
                                                          reference.patterns)),
                                              collapse =" "),
                  name = "stats.logger")
      }
      
      # 0. PREPARE ENVIRONMENT

      gK1 <- nrow(cur.centers)
      
      rev <- new.env()
      rev$modules <- list()
      rev$centers.pos <- matrix(nrow=gK1,
                                ncol=ncol(E.prep),
                                dimnames = list(
                                  paste0("c.pos", seq_len(gK1)),
                                  colnames(E.prep)))
      rev$centers.all <- matrix(nrow=gK1,
                                ncol=ncol(E.prep),
                                dimnames = list(
                                  paste0("c.all", seq_len(gK1)),
                                  colnames(E.prep)))
      
      # 1. CALCULATE CORRELATIONS -> DISTANCES -> SCORES

      # m <- corFromPrep(cur.centers, E.prep)
      m <- cosine(cur.centers, E.prep)
      
      dist.to.centers <- 1-m
      dist.to.centers[dist.to.centers < 1e-10] <- 0

      idxs <- seq_len(gK1)
      
      posScores_keeping_var <- c()
      
      # TODO: replace base to `correlation.threshold`
      nets <- lapply(idxs, function(i) {

        if(nrow(dist.to.centers) > 1) { 
          minOther <- pmin(apply(dist.to.centers[-i, , drop=F], 2, min), base) } else { 
            minOther <- base }
        
        score <- log2(minOther) - log2(dist.to.centers[i, ])
        score[score == Inf] <- 0
        score <- pmax(score, -1000)
        posScores_keeping_var <<- c(posScores_keeping_var, length(which(score>0)))

        EdgeTable <- data.table::as.data.table(data.table::copy(network.prep))
        EdgeTable[, score := score[gene]]
        EdgeTable[from > to, c("from", "to") := list(to, from)]
        EdgeTable <- EdgeTable[order(score, decreasing = T)]
        EdgeTable <- unique(EdgeTable, by=c("from", "to"))
        # we still keep loops here

        scored_graph <- igraph::graph_from_data_frame(EdgeTable, directed = F)
        igraph::V(scored_graph)$score <- 0
        scored_graph
      })

      nets_attr <- lapply(nets, mwcsr::normalize_sgmwcs_instance,
                          edges.weight.column = "score",
                          nodes.weight.column = "score",
                          edges.group.by = "gene",
                          nodes.group.by = NULL,
                          group.only.positive = TRUE)
      
      # 2. SOLVE SGMWCS TO GET MODULES
      
      # cat("Calling: batch.solver(nets)\n")
      ms <- batch.solver(nets_attr)
      # cat("Done: batch.solver(nets)\n")
      ms_mods <- lapply(ms, `[[`, "graph")
      
      # 2.a. COLLECT CORRESPONDING LOGS
      
      m.size.unique <- unlist(lapply(ms_mods, function(x) ulength(igraph::edge_attr(x)$gene)))

      if (collect.stats) {
        iter.stats_add <- data.frame(
          genes.n = dim(E.prep)[1],
          genes.pos.scored = posScores_keeping_var,
          base = base,
          m.size = unlist(lapply(ms_mods, igraph::gsize)),
          m.size.unique = m.size.unique,
          m.pos = unlist(lapply(ms_mods, function(x) sum(igraph::edge_attr(x)$score > 0))),
          m.non.neg = unlist(lapply(ms_mods, function(x) sum(igraph::edge_attr(x)$score >= 0)))
        )
        iter.stats[[iteration]] <- iter.stats_add
      }
      if (verbose) {
        flog.info(">> base was equal to: %.2g", base, 
                  name = "stats.logger")
        flog.info(">> number of modules was equal to: %s", length(ms_mods), 
                  name = "stats.logger")
        flog.info(">> sizes of modules (unique genes) were in range: %s-%s", 
                  min(m.size.unique), max(m.size.unique), 
                  name = "stats.logger")
      }

      # 2.b. RECORD MODULES AND CENTERS
      
      rev$modules <- ms_mods
      
      for (i in idxs) {

        module <- ms_mods[[i]]

        center.pos <- if (ulength(igraph::E(module)[score > 0]$gene) >= 3) {
          getCenter(E.prep, unique(igraph::E(module)[score > 0]$gene))
        } else {
          cur.centers[i, ]
        }
        center.all <- if (ulength(igraph::E(module)$gene) >= 3) {
          getCenter(E.prep, unique(igraph::E(module)$gene))
        } else {
          cur.centers[i, ]
        }
        rev$centers.pos[i, ] <- center.pos
        rev$centers.all[i, ] <- center.all
      }

      if (show.intermediate.clustering) {
        heatmapTable <- rbind(cur.centers, rev$centers.pos)[rbind(
          seq_len(gK1),
          seq_len(gK1) + gK1), ]
        pheatmap::pheatmap(
          normalize.rows(heatmapTable),
          cluster_rows=F, cluster_cols=F,
          show_rownames=T, show_colnames=F)
      }
      
      revs[[k]] <- rev
      
      # 3. DID MODULES CONVERGE, i.e. CAN WE LEAVE THE SECOND LOOP
      
      # previous iterations, in which there was the same number of modules:
      revsToCheck <- revs[sapply(revs[seq_len(k-1)], function(rev) nrow(rev$centers.pos)) 
                          == nrow(rev$centers.pos)]
      
      diff <- max(abs(rev$centers.pos - cur.centers))
      
      if (length(revsToCheck) > 0) {
        diff <- min(sapply(revsToCheck,
                           function(prevRev) max(abs(rev$centers.pos - prevRev$centers.pos))))
      }
      
      # 4. UPDATE PARAMETERS

      cur.centers <- rev$centers.pos
      
      if (!is.null(reference.patterns)) {
        flog.info("updated correlations with references (before potential merging): %s", 
                  paste0(matrixStats::colMaxs(cosine(cur.centers, 
                                                          reference.patterns)),
                         collapse =" "),
                  name = "stats.logger")
      }
      
      iteration <- iteration + 1
      
      if (verbose) {flog.info(">> max diff: %s", round(diff, 2), name = "stats.logger")}
      
      if (diff < 0.01) {break}

      k <- k + 1

    } # -------------------------------------------------------------------------------------- SECOND LOOP

    # 5. IF MODULES CONVERGED, WE CHECK THEM FOR PRESENCE OF 
    
    # (i) TOO BIG ONES:
    
    biggest.one <- max(sapply(ms_mods, function(m) ulength(igraph::E(m)$gene))) 
    
    if (biggest.one > max.module.size) {
      base <- base - base.dec * base
    }
    
    # (ii) CORRELATED ONES: 
    
    # centers.cors <- cor(t(cur.centers))
    centers.cors <- cosine(cur.centers)
    diag(centers.cors) <- 0
    correlation.max <- apply(centers.cors, 1, max, na.rm=T)
    
    if (any(correlation.max > cor.threshold)) {
      
      flog.info(">> max cor exceeded %s: %s", cor.threshold, round(max(correlation.max), 2),
                name = "stats.logger")
      max.cor.mod1 <- which.max(correlation.max) 
      max.cor.mod2 <- which.max(centers.cors[max.cor.mod1, ])
      cur.centers <- updCenters(cur.centers = cur.centers, 
                                m1 = max.cor.mod1, m2 = max.cor.mod2, 
                                E.prep = E.prep, ms_mods = ms_mods)
    } else {
      
      # (iii) or UNINFORMATIVE ONES:
      
      gesecaRes <- doGeseca(E.prep = E.prep,
                            network.prep = network.prep,
                            network.annotation = network.annotation,
                            modules = rev$modules,
                            scale = FALSE,
                            center = FALSE,
                            verbose = verbose,
                            gesecaSeed = 0)
      
      # TODO: geseca is randomized, there can be instabilities of comparing with a threshold
      # currently using sampleSize=1001 for more stability, but size-based thresholds
      # can be precomputed
      good <- gesecaRes$pathway[which(gesecaRes$padj < p.adj.val.threshold)]
      bad <- rownames(cur.centers)[!rownames(cur.centers) %in% good]

      if (length(bad) == 0 & biggest.one > max.module.size) {next} 
      if (length(bad) == 0 & biggest.one <= max.module.size) {break} 
      if (length(bad) != 0 & nrow(cur.centers) > 1) {
        
        m <- cosine(cur.centers, E.prep)
        is.positive <- (1 - m < base)
        # number of genes with potentially positive scores in two modules
        gene.overlaps <- crossprod(t(is.positive)) 
        diag(gene.overlaps) <- 0
        
        if (max(gene.overlaps[bad, ]) > 0) {
          # There is some overlap, we can try to merge a bad module with the good one,
          # without losing the genes from the bad one (and not destroying good one too much).
          # It also means that the centers are pretty similar.
          max.cor.mod1 <- as.integer(gsub("c.pos", "", bad[which.max(apply(gene.overlaps, 1, max, na.rm=T)[bad])]))
          max.cor.mod2 <- which.max(gene.overlaps[max.cor.mod1, ])
          cur.centers <- updCenters(cur.centers = cur.centers, 
                                    m1 = max.cor.mod1, m2 = max.cor.mod2, 
                                    E.prep = E.prep, ms_mods = ms_mods)
        } else {
          # No point in merging, let's remove module with the least number of positive genes. 
          to.remove <- names(which.min(rowSums(is.positive)[bad]))
          cur.centers <- cur.centers[rownames(cur.centers) != to.remove, ]
        }
      } else {
        saveStats(work.dir, rev, gesecaRes, iter.stats)
        message.string <- "[Attention!] No modules found.\n
          Try to tune method's parameters.\n
          Check '/stats' folder for the statistics of the run."
        flog.warn(message.string, name = "stats.logger")
        warning(message.string, call. = FALSE)
        return()
      }
    }

    # keep expressions devoted to sizes of modules:
    # m.sizes <- sapply(modules, function(m) ulength(igraph::E(m)$gene))
    # modules <- modules[m.sizes >= min.module.size] # add as param

  } # ---------------------------------------------------------------------------------------- FIRST LOOP
  
  # 9. FINAL ADJUSTMENTS OF MODULES 
  
  # (i) compactise
  modules_pre <- lapply(rev$modules, function(x) {
    igraph::graph.attributes(x)$signals[which(names(igraph::graph.attributes(x)$signals) %in%
                                                igraph::vertex_attr(x)$signal)] <- -0.001
    x
  })
  modules_set <- batch.solver(modules_pre)
  modules <- lapply(modules_set, `[[`, "graph")

  # (ii) recalc geseca & sort modules
  gesecaRes <- doGeseca(E.prep = E.prep,
                        network.prep = network.prep,
                        network.annotation = network.annotation,
                        modules = modules,
                        scale = FALSE,
                        center = FALSE,
                        verbose = verbose,
                        gesecaSeed = 0)
  
  order_idx <- as.numeric(gsub("c.pos", "", gesecaRes$pathway))
  
  modules <- modules[order_idx]
  nets_attr <- nets_attr[order_idx]
  rev$centers.pos <- rev$centers.pos[order_idx, , drop = FALSE]
  rownames(rev$centers.pos) <- paste0("c.pos", 1:nrow(rev$centers.pos)) 
  rev$centers.all <- rev$centers.all[order_idx, , drop = FALSE]
  rownames(rev$centers.all) <- paste0("c.all", 1:nrow(rev$centers.all))  
  
  gesecaRes$pathway <- paste0("m", 1:nrow(gesecaRes)) 

  saveStats(work.dir, rev, gesecaRes, iter.stats)
  
  flog.info("GAM-CLUSTERING ends here.", name = "stats.logger")

  return(list(
    modules = modules,
    nets = nets_attr,
    patterns.pos = rev$centers.pos,
    patterns.all = rev$centers.all,
    iter.stats = iter.stats
  ))
}
