#!/usr/bin/env Rscript

# This script is based on https://github.com/fml-fam/fmlr/blob/master/vignettes/rebuild.r

library(rmarkdown)

rmf <- function(f) {
  if (file.exists(f))
    file.remove(f)
}

clean <- function() {
  files = dir(pattern="*.Rmd", recursive=FALSE)
  for (f in files)
    rmf(f)
}

set_path <- function() {
  while (!file.exists("DESCRIPTION"))
  {
    setwd("..")
    if (getwd() == "/")
      stop("couldn't find package!")
  }
  
  setwd("vignettes")
}


build_vignette <- function(src, quiet=TRUE) {
  out_name <- sub("^_", "", tools::file_path_sans_ext(basename(src)))
  fig_dir  <- paste0(out_name, "_files/")
  
  vignettes_dir <- normalizePath(".")
  root_dir <- file.path(vignettes_dir, "work")
  dir.create(root_dir, showWarnings = FALSE)
  
  knitr::opts_knit$set(base.dir = vignettes_dir, base.url = "", root.dir = root_dir)
  
  # fig.cap = "" otherwise it renders captions in Rmd
  knitr::opts_chunk$set(fig.path = fig_dir, fig.width = 6, fig.height = 5, fig.cap = "")
  # but still want automatic alt text
  knitr::opts_hooks$set(fig.cap = function(options) {
    if (identical(options$fig.cap, "") && is.null(options$fig.alt))
      options$fig.alt <- paste("plot of chunk", options$label)
    options
  })
  
  
  # patch include_graphics to copy code-generated images into fig_dir
  local({
    orig <- get("include_graphics", envir = asNamespace("knitr"))
    fig_dir_abs <- file.path(vignettes_dir, fig_dir)
    assignInNamespace("include_graphics", function(path, ...) {
      abs <- normalizePath(path, mustWork = FALSE)
      # NB: this flattens all the files
      if (!startsWith(abs, paste0(fig_dir_abs, .Platform$file.sep)) && file.exists(abs)) {
        dir.create(fig_dir_abs, recursive = TRUE, showWarnings = FALSE)
        file.copy(abs, file.path(fig_dir_abs, basename(abs)), overwrite = TRUE)
        abs <- file.path(fig_dir_abs, basename(abs))
      }
      # Pass the path already relative to vignettes_dir so knitr's own
      # input_dir()-based relativisation and existence check don't interfere.
      orig(xfun::relative_path(abs, vignettes_dir), rel_path = FALSE, error = FALSE, ...)
    }, ns = "knitr")
  })
  
  knitr::knit(src, output = paste0(out_name, ".Rmd"), quiet = quiet)
  message("Built: ", out_name, ".Rmd  (figures -> ", fig_dir, ")")
}

set_path()
# clean()

build_vignette("./src/_Algorithm_overview.Rmd", quiet=FALSE)
# build_vignette("./src/_GAMclust_tutorial_BULK.Rmd", quiet=FALSE)
# build_vignette("./src/_GAMclust_tutorial_SC.Rmd", quiet=FALSE)
# build_vignette("./src/_GAMclust_tutorial_SPAT.Rmd", quiet=FALSE)
