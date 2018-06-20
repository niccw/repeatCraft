require(shiny)

#' @export
runApp <- function() {
  exampleDir <- system.file("shinyapp/exampleData", package = "repeatcraft")
  if (example == "") {
    stop("Could not find example directory. Try re-installing `repeatcraft`.", call. = FALSE)
  }

  appDir <- paste0(system.file(package = "repeatcraft"),"/R")

  shiny::runApp(appDir, display.mode = "normal")
}
