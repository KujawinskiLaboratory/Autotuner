#' @title exportPeaksUI
#'
#' @description This function contains the UI component of the export Peak
#' portion of the shiny app used within Autotuner. This is the part where a
#' spreadsheet is saved as output.
#'
#' @param id - character string used to map UI and server within App
#'
#' @import shiny
#' @export
exportPeaksUI <- function(id) {
  ns <- NS(id)
  uiOutput(ns("exportPeaksCentralUI"))
}

#' @title exportPeaks
#'
#' @description This function contains the server component of the export Peak
#' portion of the shiny app used within Autotuner. This is the part where a
#' spreadsheet is saved as output.
#'
#' @param input - List of UI entered values.
#' @param output - Values exported from server to the UI after computation.
#' @param session - String used to match the namespace of the UI and server
#' @param filledPeaks - The xcms object to get exported as a feature table.
#'
#' @export
exportPeaks <- function(input, output, session, filledPeaks) {

  ns <- session$ns

  ## Central UI containing all information
  output$exportPeaksCentralUI <- renderUI({
    req(filledPeaks())
    tagList(
      titlePanel("Export Aligned Peak Table"),
      wellPanel(
        textInput(ns("fileName"), label = "Enter name of export file:", value = "alignedPeaks.csv"),
        actionButton(ns("exportPeaks"), "Export Peak Table")
      )
    )
  })

  observeEvent(input$exportPeaks, {
    fileName <- isolate(input$fileName)
    write.csv(x = xcms::peakTable(filledPeaks()), file = fileName, row.names = F)
  })

}
