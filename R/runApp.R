require(shiny)

#' @import shiny DT dplyr magrittr plotly readr reshape2 shinythemes stringr
#' @importFrom Biostrings DNAStringSet substr
#' @importFrom rtracklayer import.gff import.gff3 mcols start end export.gff3 
#' @importFrom GenomicRanges seqnames GRanges makeGRangesFromDataFrame
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
