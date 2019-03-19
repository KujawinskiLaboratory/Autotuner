#' @title fineParamsUI
#'
#' @description This function contains the UI component of the EIC based
#' parameter selection portion of the shiny app used within Autotuner.
#'
#' @param id - character string used to map UI and server within App
#'
#' @export

fineParamsUI <- function(id) {

  ns <- NS(id)
  tagList(
    titlePanel("Estimate parameters From Scan and EIC Data"),
        wellPanel(
          h2("Parameters determined with TIC data:"),
          tableOutput(ns("ticTable")),
          h2("Parameters determined with EIC data:"),
          tableOutput(ns("eicTable")),
          uiOutput(ns("parambutton")),
          br(),
          uiOutput(ns("exportParams"))
        )
  )
}

#' @title fineParams
#'
#' @description This function contains the server component of the EIC based
#' parameter selection portion of the shiny app used within Autotuner. This
#' module is specifically involved in presentation of recommended EIC peaks
#' following the completion of Autotuner. It can help the user export a
#' the calculated parameters as a spreadsheet.
#'
#' @param input - List of UI entered values.
#' @param output - Values exported from server to the UI after computation.
#' @param session - String used to match the namespace of the UI and server
#' @param params - The reactive object containing all calculated parameter
#' estimates from the data.
#' @param runfile - A data.frame containing metadata for all the samples within
#' the selected folder.
#' @param Autotuner - Data with specifics on samples path and metadata.
#'
#' @export
fineParams <- function(input, output, session, params, runfile, Autotuner) {

    ns <- session$ns

    # Rendering tables of calculated parameters -------------------------------
    # renders UI for table of recommended parameters
    output$ticTable <- renderTable({

        req(params())
        max_missing <- max(table(runfile[,Autotuner@factorCol]))
        x <- data.frame(t(params()[[2]]), max_missing)
        colnames(x) <- c("Max Peakwdith", "Min Peakwidth", "Group Difference",
                         "Max Missing")
        x

    }, include.rownames=FALSE)


    # renders UI for table of recommended parameters
    output$eicTable <- renderTable({

        req(params())
        params()[[1]]

    }, include.rownames=FALSE)

    # Exporting Table of Parameters -------------------------------------------
    ## Central UI containing all information
    output$exportParams <- renderUI({

        req(params())
        tagList(
            textInput(ns("fileName"),
                      label = "Enter name of export file:",
                      value = "autotuner_parameters.csv"),
            actionButton(ns("exportParams"), "Export Calculated Parameters")
        )
    })

    observeEvent(input$exportParams, {

        fileName <- isolate(input$fileName)
        eicTable <- params()[[1]]
        ticTable <- params()[[2]]
        allParamsEsts <- cbind(eicTable, t(ticTable))

        write.csv(x = allParamsEsts,
                  file = file.path(path.expand("~"), fileName), row.names = F)

    })

}
