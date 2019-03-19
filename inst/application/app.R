# UI ----------------------------------------------------------------------
ui <- fluidPage(
    navbarPage("Autotuner",
               tabPanel("Load Data",
                        loadDataUI("load")
               ),
               tabPanel("Find Peaks",
                        signalProcessingUI("signal")
               ),
               tabPanel("Visualize Peaks",
                        peakVisUI("peakVisualization")
               ),
               tabPanel("Estimated Parameters",
                        fineParamsUI("fine")
                        )
    )
) # end of fluid page


# Server ------------------------------------------------------------------
server <- function(input, output) {

    loadedContent <- callModule(loadData, "load")
    calcSignal <- callModule(module = signalProcessing,
                             id = "signal",
                             Autotuner = loadedContent()$Autotuner)

    params <- callModule(peakVis, "peakVisualization", calcSignal,
                         loadedContent()$Autotuner)


    paramList <- callModule(fineParams, "fine", params,
                            loadedContent()$runfile,
                            loadedContent()$Autotuner)
}

shinyApp(ui = ui, server = server)
