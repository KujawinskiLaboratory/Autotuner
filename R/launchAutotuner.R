#' @title launchAutotuner
#'
#' @description This function deploys the Autotuner shiny app.
#'
#' @export
launchAutotuner <- function() {
    shiny::runApp(appDir = system.file("application/", package = "Autotuner"))
}
