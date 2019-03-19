# Outer UI function --------------------------------------------------------
#' @title signalProcessingUI
#'
#' @description This function contains the UI component of the signal processing
#' portion of the shiny app used within Autotuner.
#'
#' @param id - character string used to map UI and server within App
#'
#' @export
signalProcessingUI <- function(id) {
  ns <- NS(id)
  tagList(

    titlePanel("Finding Peaks within Chromatograms
               Using Sliding Window Analysis"),
    column(3,
           wellPanel(
                uiOutput(ns("sampleSelector")),
                uiOutput(ns("signalParams")),
                uiOutput(ns("goodSignals"))
           )
    ),
    column(9,
           wellPanel(
                plotOutput(ns("signalPlot"))
           )
        )
    )
}

# Outer Server Function ---------------------------------------------------
#' @title signalProcessing
#'
#' @description This function contains the server component of the signal
#' processing portion of the shiny app used within Autotuner.
#'
#' @param input - List of UI entered values.
#' @param output - Values exported from server to the UI after computation.
#' @param session - String used to match the namespace of the UI and server
#' @param Autotuner - An Autotuner objected containing sample specific raw
#' data.
#'
#' @export
signalProcessing <- function(input, output, session, Autotuner) {

    ns <- session$ns

    ## Checking if load data is completed before moving forward
    inputCheck <- reactiveValues(miss = F)

    ## checks to make sure load data has been placed
    observe({

        req(Autotuner)
        inputCheck$miss <- T

    })


    ## UI to load samples into app
    output$sampleSelector <- renderUI({

        validate(need(inputCheck$miss,
                      "Files must first be loaded to use find peaks."))
        selectizeInput(ns("sample"),
                       "Choose samples to visualize:",
                       choices = Autotuner@metadata[,Autotuner@file_col],
                       multiple = TRUE)
    })

    ## generating UI for picking sliding window parameters
    output$signalParams <- renderUI({

        validate(need(inputCheck$miss, "\n"))
        tagList(
            sliderInput(ns("lag"), "Lag:", min = 5,
                        max = 100, value = 10, step = 1),
            sliderInput(ns("threshold"),
                        "Significance Threshold (sd):",
                        min = 0, max = 5, value = 3.5, step = 0.1),
            numericInput(ns("influence"), "Influence:",
                         value = 0.01, min = 0, max = 1,
                         step = 0.01)
        )
    })

    signalPar <- reactive({
        return(list(lag = input$lag,
                    threshold = input$threshold,
                    influence = input$influence))
    })

    ## function used to calculate signal across the trace
    signals <- reactive({

        req(inputCheck$miss)
        req(!is.null(input$lag))


        signals <- lapply(Autotuner@intensity,
                                   ThresholdingAlgo,
                                   lag = input$lag,
                                   threshold = input$threshold,
                                   influence = input$influence)

        return(signals)

    })

    ## UI button to output signals
    output$goodSignals <- renderUI({

        req(inputCheck$miss)
        actionButton(ns("findPeakAction"),
                   label = "Accept Observed Peaks")

    })

    ## Determining which samples to visualize
    sampleIndex <- reactive({

        req(input$sample)
        sampleCol <- grep(Autotuner@file_col, colnames(Autotuner@metadata))
        temp <- which(unlist(Autotuner@metadata[,sampleCol]) %in% input$sample)
        return(temp)

    })

    ## Generating UI to plot signals
    output$signalPlot <- renderPlot({

        ## the specific samples that are entered into the data
        validate(need(length(input$sample) > 0, "Please select samples"))
        req(signals())

        signalsOut <- plot_signals(Autotuner = Autotuner,
                                threshold = input$threshold,
                                sample_index = sampleIndex(),
                                signals = signals())
        return(signalsOut)
    })


    acceptedSignal <- eventReactive(input$findPeakAction, {
        return(signals())
    })

    return(acceptedSignal)
}
