require(shiny)

#' @import shiny Biostrings DT GenomicFeatures GenomicRanges TnT dplyr magrittr plotly readr reshape2 rtracklayer shinythemes stringr
#' @importFrom shinyjs hidden
#' @importFrom systemPipeR predORF
NULL

#' @export
runApp <- function() {
  exampleDir <- system.file("shinyapp/exampleData", package = "repeatcraft")
  if (exampleDir == "") {
    stop("Could not find example directory. Try re-installing `repeatcraft`.", call. = FALSE)
  }

  appDir <- system.file("shinyapp",package = "repeatcraft")

  shiny::runApp(appDir, display.mode = "normal")
}
