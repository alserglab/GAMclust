#' Get files with modules' graphs
#'
#' Creates graphical representations of the identified metabolic modules.
#' For each module, the function annotates genes and metabolites, generates
#' network layouts, and exports the results in multiple formats suitable for
#' visualization and downstream analysis.
#' Supported outputs include PDF, PNG, SVG, GraphML, and XGMML files.
#' These visualizations facilitate biological interpretation of the detected
#' metabolic subnetworks.
#'
#' @param modules Metabolic modules.
#' @param network.annotation Metabolic network annotation.
#' @param metabolites.annotation Metabolites annotation.
#' @param work.dir Working directory where results should be saved.
#' @param seed.for.layout Random seed required for gatom to draw metabolic graph.
#' each atom to corresponding metabolite's name
#' @return Results of this function can be seen in work.dir (.pdf, .png and .xgmml files).
#' @export
getGraphs <- function(modules,
                      network.annotation,
                      metabolites.annotation,
                      seed.for.layout = 22,
                      work.dir = work.dir){

  for (i in seq_along(modules)) {

    dir.create(work.dir, showWarnings=F)
    m <- modules[[i]]
    name <- sprintf("%s.%s", "m", i)

    # expand edge/node attributes for drawing:
    # at the moment it has: from | to | gene | score | signal | index
    igraph::E(m)$label <- network.annotation$genes[igraph::E(m)$gene]$symbol
    igraph::E(m)$pval <- 0.01
    igraph::E(m)$log2FC <- igraph::E(m)$score
    # at the moment it has: name | score | signal | index
    # take only metabolite id without atom coord
    igraph::V(m)$label <- sapply(strsplit(igraph::V(m)$name, split = "_"), `[[`, 1)
    # fix 01.04.2022 + [, on = .(metabolite)]
    igraph::V(m)$label <- metabolites.annotation$metabolites[igraph::V(m)$label, on = .(metabolite)]$metabolite_name
    igraph::V(m)$logPval <- -5
    igraph::V(m)$log2FC <- NA

    file <- file.path(work.dir, paste0(name, ".dot"))
    gatom::saveModuleToDot(m, file=file, name=name)
    gatom::saveModuleToXgmml(m, name=name, file=file.path(work.dir, paste0(name, ".xgmml")))
    igraph::write_graph(m, file=file.path(work.dir, paste0(name, ".graphml")), format="graphml")
    # system2("neato", c("-Tpdf", "-o", file.path(work.dir, paste0(name, ".pdf")), file), stderr = F)
    # system2("neato", c("-Tpng", "-o", file.path(work.dir, paste0(name, ".png")), file), stderr = F)
    # use Rgraphviz instead of command line graphviz
    mm <- Rgraphviz::agread(file, layoutType="dot", layout=TRUE)
    Rgraphviz::toFile(mm, layoutType="neato", filename=file.path(work.dir, paste0(name, ".svg")), fileType="svg")
    tryCatch(
      rsvg::rsvg_pdf(file.path(work.dir, paste0(name, ".svg")), file.path(work.dir, paste0(name, ".pdf"))),
      error = function(e) NULL)
    tryCatch(
      rsvg::rsvg_png(file.path(work.dir, paste0(name, ".svg")), file.path(work.dir, paste0(name, ".png"))),
      error = function(e) NULL)
    tryCatch(
      gatom::saveModuleToPdf(module = m, seed = seed.for.layout,
                             file = file.path(work.dir, paste0(name, "_seed", seed.for.layout, ".pdf"))),
      error = function(e) NULL)
    messagef("Graphs for module %s are built", i)
  }
}


#' Get files with modules' genes
#'
#' Generates gene-level summaries for each identified metabolic module.
#' The function exports tables containing module genes, positively scored
#' genes that were not selected into the final module, and the top correlated
#' genes from the complete expression dataset.
#' Each table includes gene identifiers, symbols, module scores, and
#' correlations with module expression patterns. These outputs provide a
#' detailed view of the transcriptional composition of every module.
#'
#' @param modules Metabolic modules.
#' @param nets Scored networks.
#' @param patterns Patterns of metabolic modules.
#' @param gene.exprs Gene expression.
#' @param network.annotation Metabolic network annotation.
#' @param work.dir Working directory where results should be saved.
#' @return m.*.genes.tsv -- module genes.
#' @return m.*.notInModule.genes.tsv -- genes with positive score not included into module.
#' @return m.*.complete.genes.tsv -- top 300 of all genes sorted by correlation value.
#' @return Results of this function can be seen in work.dir (three .tsv files for each module with gene lists).
#' @export
getGeneTables <- function(modules, nets, patterns, gene.exprs,
                          network.annotation,
                          work.dir = work.dir){
  
  m.gene.list <- list()

  dir.create(work.dir, showWarnings=F)
  
  for (i in seq_along(modules)) {

    m <- modules[[i]]
    net <- nets[[i]]

    t <- igraph::as_data_frame(m)[, c("gene", "score")]
    t <- t[!duplicated(t$gene), ]
    colnames(t)[1] <- "Entrez"
    t$symbol <- network.annotation$genes[t$Entrez]$symbol
    t$cor <- cor(patterns[i, ], t(gene.exprs[t$Entrez, ]))[1,]
    t <- t[order(t$cor, decreasing = T), c("symbol", "Entrez", "score", "cor")]
    write.tsv(t, file=sprintf("%s/m.%s.genes.tsv", work.dir, i))
    m.gene.list[[i]] <- t$symbol
    names(m.gene.list[[i]]) <- t$Entrez

    notInModule <- data.table::data.table(Entrez=setdiff(igraph::E(net)$gene, igraph::E(m)$gene))
    notInModule[, score := igraph::E(net)[match(Entrez, gene)]$score]
    notInModule[, symbol := network.annotation$genes[match(Entrez, gene)]$symbol]
    notInModule[, cor := cor(patterns[i, ], t(gene.exprs[Entrez, ]))[1,]]
    notInModule <- notInModule[order(cor, decreasing=T), c("symbol", "Entrez", "score", "cor")]
    write.tsv(notInModule[score > 0], file=sprintf("%s/m.%s.notInModule.genes.tsv", work.dir, i))

    notInModule2 <- data.table::data.table(Entrez=rownames(gene.exprs))
    notInModule2[, symbol := network.annotation$genes[match(Entrez, gene)]$symbol]
    notInModule2[, cor := cor(patterns[i, ], t(gene.exprs[Entrez, ]))[1,]]
    notInModule2 <- notInModule2[order(cor, decreasing=T), c("symbol", "Entrez", "cor")]
    write.tsv(notInModule2[1:300], file=sprintf("%s/m.%s.complete.genes.tsv", work.dir, i))
    messagef("Gene tables for module %s are produced", i)
  }
  
  names(m.gene.list) <- paste0("m_", seq_along(modules))
  m.gene.list
}

#' Get files with modules' annotations
#'
#' Annotates metabolic modules using pathway enrichment analysis. For each
#' module, over-representation analysis is performed against the pathway
#' collection contained in the network annotation, followed by pathway
#' redundancy reduction.
#' The function exports pathway tables containing significantly enriched
#' biological processes and their associated module genes. These annotations
#' help connect discovered modules to known metabolic functions.
#'
#' @param network.annotation Metabolic network annotation.
#' @param nets Scored networks.
#' @param work.dir Working directory where files with module genes are (results will be saved here as well).
#' @param padj.threshold Threshold, adjusted p-value, for pathways.
#' @return Results of this function can be seen in work.dir (.tsv files).
#' @export
getAnnotationTables <- function(network.annotation, nets, work.dir,
                                 padj.threshold = Inf){

  g_files <- list.files(work.dir, "m\\.[0-9]+\\.genes\\.tsv", full.names = T)
  if(length(g_files) == 0){stop("Use `annotateModules()` after `getGeneTables()` only.")}

  for (i in 1:length(g_files)){

    file <- data.table::fread(g_files[i])
    name <- basename(g_files[i])
    name <- sapply(strsplit(name, ".", fixed = T),"[[", 2)
    net <- nets[[as.numeric(name)]]

    foraRes <- fgsea::fora(pathways=network.annotation$pathways,
                           genes=file$Entrez,
                           universe=unique(igraph::as_data_frame(net)$gene))
    foraResCllpsd <- fgsea::collapsePathwaysORA(foraRes = foraRes[order(pval)][padj < padj.threshold],
                                                pathways = network.annotation$pathways,
                                                genes = file$Entrez,
                                                universe = unique(igraph::as_data_frame(net)$gene))
    mainPathways <- foraRes[pathway %in% foraResCllpsd$mainPathways][order(-overlap), ]

    pz <- lapply(mainPathways$overlapGenes, function(x) unlist(strsplit(x, " ")))
    new_pz <- lapply(pz, function(set) network.annotation$genes[gene %in% set]$symbol)

    if(length(new_pz) != 0){

      for(i in 1:length(new_pz)){
        mainPathways[i, "overlapGenes"] <- paste(new_pz[[i]], collapse=" ")}
      mainPathways <- data.frame(lapply(mainPathways, as.character), stringsAsFactors=F)

      write.tsv(mainPathways, file=sprintf("%s/m.%s.pathways.tsv", work.dir, name))
      messagef("Pathway annotation for module %s is produced", name)

    } else {
      messagef("No pathways annotating module %s are found", name)
    }
  }
}


#' Heatmap for annotated pathways
#'
#' Builds a heatmap summarizing pathway annotations across all identified
#' modules. The visualization displays the proportion of module genes
#' associated with each enriched pathway, enabling rapid comparison of
#' functional themes between modules.
#' Filtering options allow users to restrict displayed pathways based on
#' adjusted p-values and gene overlap thresholds. The resulting heatmap
#' provides a compact overview of the functional landscape revealed by
#' GAM-clustering.
#'
#' @param work.dir Directory with gene and pathways tables
#' @param padj.threshold Set threshold for p-adjusted
#' @param threshold Minimal percent of genes in the pathway
#' @param file_name Name of file with heatmap
#' @return Results of this function can be seen in work.dir (png files)
#' @export
getAnnotationHeatmap <- function(work.dir,
                                 padj.threshold = Inf,
                                 threshold = 0.05,
                                 file_name="Modules_heatmap.png") {

  # files with pathways
  # reorder because if >9 default order is lexicographic: 1, 10, 11, 2, 3....
  m_files <- list.files(work.dir, "m\\.[0-9]+\\.pathways\\.tsv", full.names = T)
  m_files <- m_files[order(as.numeric(gsub(".*m\\.(\\d+).*", "\\1", m_files)))]

  # files with genes
  # TODO: should we enrich on positively-scored genes only?
  m_files_genes <- list.files(work.dir, "m\\.[0-9]+\\.genes\\.tsv", full.names = T)
  m_files_genes <- m_files_genes[order(as.numeric(gsub(".*m\\.(\\d+).*", "\\1", m_files_genes)))]

  # if exist
  if (length(m_files_genes) == 0) {
    stop("Use 'getGeneTables' first.", call. = FALSE)
  } else if (length(m_files) == 0) {
    stop("No pathways annotating any of the modules. Try to use 'getModulesAnnotation' first.", call. = FALSE)
  }

  m_files_genes <- m_files_genes[sort(as.numeric(gsub(".*m\\.(\\d+).*", "\\1", m_files)))]

  set <- list()

  # options(readr.num_columns = 0)
  for (i in 1:length(m_files)) {
    messagef("Processing module %s", as.numeric(gsub(".*m\\.(\\d+).*", "\\1", m_files[i])))
    file <- read.tsv(m_files[i])
    moduleSize <- nrow(read.tsv(m_files_genes[i]))
    messagef(sprintf("Module size: %s", moduleSize))
    if (length(file$overlap) == 0) {
      paths <- file[, c("pathway", "overlap")]
      paths$percentOfGenesInThePthw <- c(0)
      paths$PATHNAME <- c("No pathways")
    } else {
      paths <- file[which(file$padj < padj.threshold), ]
      paths <- paths[, c("pathway", "overlap")]
      paths <- paths[order(paths$overlap, decreasing = T), ]
      paths$percentOfGenesInThePthw <- paths$overlap / moduleSize
      paths <- paths[which(paths$percentOfGenesInThePthw > threshold), ]
    }
    set <- c(set, list(paths))
  }

  pathNames <-c()
  for(i in seq_along(set)){
    pathNames <- c(pathNames, set[[i]]$pathway)
  }
  pathNames <- unique(pathNames)

  set_df <- as.data.frame(matrix(NA, nrow = length(set), ncol=length(pathNames)))
  colnames(set_df) <- pathNames

  for(i in seq_along(set)){
    set_df[i, set[[i]]$pathway] <- set[[i]]$percentOfGenesInThePthw
  }

  if (nrow(set_df) == 0) {
    stop("No pathways were found", call. = FALSE)
  }

  set_df[is.na(set_df)] <- 0

  rownames(set_df) <- c(paste("Module", sort(as.numeric(gsub(".*m\\.(\\d+).*", "\\1", m_files)))))
  # for heatmap height calculation
  max_len <- max(nchar(colnames(set_df)))

  pheatmap::pheatmap(
    t(set_df),
    cluster_rows=F, cluster_cols=F,
    show_rownames=T, show_colnames=T,
    # width=ifelse(ncol(set_df) < 10, 7, ncol(set_df) * 0.35),
    # height=ifelse(max_len < 90, 8, ceiling(max_len/11)),
    cellwidth=30,
    # file=file.path(work.dir, file_name),
    breaks = seq(0, 1, 0.05),
    color = colors4heatmap)

}


#' Compare two GAM-clustering runs
#'
#' Compares the results of two independent GAM-clustering analyses. The
#' function evaluates module overlap, similarity indices, and, when
#' applicable, correlations between module expression patterns.
#' Comparisons can be performed using either all module genes or only
#' positively scored genes. A combined summary figure is generated to
#' facilitate assessment of module reproducibility and consistency across
#' datasets or parameter settings.
#'
#' @param dir1 Folder with GAM-clustering results of the first run.
#' @param dir2 Folder with GAM-clustering results of the second run.
#' @param name1 Name of the first run.
#' @param name2 Name of the second run.
#' @param same.data Whether two runs were made on the same data.
#' @param use.genes.with.pos.score Whether to build figure considering positive genes only.
#' @param work.dir Folder where final figure should be saved in.
#' @param file.name Name of the final figure.
#' @param plot.height Height of the resulting plot in inches.
#' @param plot.width Width of the plot in inches.
#' @return Results of this function can be seen in work.dir (.pdf, .png and .xgmml files).
#' @export
modulesSimilarity <- function(dir1, dir2,
                              name1 = basename(dir1),
                              name2 = basename(dir2),
                              same.data = FALSE,
                              use.genes.with.pos.score = FALSE,
                              work.dir = getwd(),
                              plot.height = 10,
                              plot.width = NULL,
                              file.name = sprintf("%s_VS_%s.png", name1, name2)) {
  # number of cols
  n_col <- ifelse(same.data, 3, 2)
  
  # intersection and index for all genes
  inter <- compareModules(dir1, dir2,
                          file.name = sprintf("%s/intersection.png", work.dir))
  index <- compareModulesIndex(dir1, dir2,
                               file.name = sprintf("%s/index.png", work.dir))
  
  # intersection and index for pos.genes
  if (use.genes.with.pos.score) {
    inter_pos <- compareModules(dir1, dir2,
                                use.genes.with.pos.score.only = use.genes.with.pos.score,
                                file.name = sprintf("%s/intersection_pos.png", work.dir))
    index_pos <- compareModulesIndex(dir1, dir2,
                                     use.genes.with.pos.score.only = use.genes.with.pos.score,
                                     file.name = sprintf("%s/index_pos.png", work.dir))
  }
  
  # correlation
  if (same.data) {
    cor_all <- compareModulesCor(paste0(dir1, "/stats/patterns_all.tsv"), paste0(dir2, "/stats/patterns_all.tsv"),
                                 file.name = sprintf("%s/correlation_all.png", work.dir))
    if (use.genes.with.pos.score) {
      cor_p <- compareModulesCor(paste0(dir1, "/stats/patterns_pos.tsv"), paste0(dir2, "/stats/patterns_pos.tsv"),
                                 file.name = sprintf("%s/correlation_pos.png", work.dir))
    }
  }
  
  # get final plots
  if (use.genes.with.pos.score) {
    if (same.data) {
      plots <- cowplot::plot_grid(inter$gtable, index$gtable, cor_all$gtable,
                                  inter_pos$gtable, index_pos$gtable, cor_p$gtable,
                                  ncol = n_col)
    } else {
      plots <- cowplot::plot_grid(inter$gtable, index$gtable,
                                  inter_pos$gtable, index_pos$gtable,
                                  ncol = n_col)
    }
  } else {
    if (same.data) {
      plots <- cowplot::plot_grid(inter$gtable, index$gtable, cor_all$gtable,
                                  ncol = n_col)
    } else {
      plots <- cowplot::plot_grid(inter$gtable, index$gtable,
                                  ncol = n_col)
    }
  }
  
  title <- cowplot::ggdraw() +
    cowplot::draw_label(
      sprintf("Rows: %s / Columns: %s", name1, name2),
      fontface = 'bold',
      x = 0,
      hjust = 0,
      size = 24
    ) + ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 7)
    )
  
  final <- cowplot::plot_grid(title, plots,
                              ncol = 1,
                              rel_heights = c(0.1, 1)) + ggplot2::theme_bw()
  
  cowplot::save_plot(plot = final, filename = sprintf("%s/%s", work.dir, file.name),
                     base_height = plot.height,
                     base_width = plot.width)
}


