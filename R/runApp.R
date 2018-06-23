require(shiny)

#' @import shiny Biostrings DT GenomicFeatures GenomicRanges dplyr magrittr plotly readr reshape2 rtracklayer shinythemes stringr GenomeInfoDb
#' @importFrom shinyjs hidden
#' @importFrom systemPipeR predORF
NULL

#' @export
runApp <- function() {
  exampleDir <- system.file("shinyapp", package = "repeatcraft")
  if (exampleDir == "") {
    stop("Could not find app directory. Try re-installing `repeatcraft`.", call. = FALSE)
  }

  if(!require(TnT)){
    print("Installing package TnT...")
    source("https://bioconductor.org/biocLite.R")
    biocLite("TnT")
  }

  appDir <- system.file("shinyapp",package = "repeatcraft")

  shiny::runApp(appDir, display.mode = "normal")
}
