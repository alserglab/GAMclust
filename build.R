# Pre-render vignettes
devtools::load_all(".")
knitr::opts_chunk$set(error = FALSE)
knitr::knit(
  input = "vignettes/src/_GAMclust_tutorial_BULK.Rmd",
  output = "vignettes/GAMclust_tutorial_BULK.Rmd"
)
knitr::knit(
  input = "vignettes/src/_GAMclust_tutorial_SC.Rmd", 
  output = "vignettes/GAMclust_tutorial_SC.Rmd"
)
knitr::knit(
  input = "vignettes/src/_GAMclust_tutorial_SPAT.Rmd", 
  output = "vignettes/GAMclust_tutorial_SPAT.Rmd"
)

devtools::document()

# Build by component
pkgdown::check_pkgdown()
pkgdown::build_home()
pkgdown::build_articles_index()
pkgdown::preview_site()

# Build ALL
pkgdown::build_site()
