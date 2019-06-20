#'  @title peakVisUI
#'
#'  @description This file contains all the code needed to run the peak
#'  visualization from the metabolomics data analysis shiny page. The main
#'  functions for constructing UI and server behavior are teh "peakVisUI" and
#'  "peakVis" respectively. All of the other functions depend on
#'
#' @param id - character string used to map UI and server within App
#'
#' @export
## outer peakvis page UI
peakVisUI <- function(id) {
  ns <- NS(id)
  tagList(
    titlePanel("Visualizing Peaks and Estimating Peak Widths"),
    column(3,
           wellPanel(
             uiOutput(ns("peakNumber")),
             uiOutput(ns("plotCode")),
             uiOutput(ns("bound")),
             uiOutput(ns("thresh")),
             uiOutput(ns("button"))
           )
    ),
    column(9,
           wellPanel(
               uiOutput(ns("peakPlot")),
               uiOutput(ns("status"))
           ))
  )
}

#' @title peakVis
#'
#' @description This function contains the server component of the peak
#' visualization portion of the shiny app used within Autotuner. In this section
#' the user selects if signal processing parameters were properly tuned to
#' indentify individual TIC peaks.
#'
#' @param input - List of UI entered values.
#' @param output - Values exported from server to the UI after computation.
#' @param session - String used to match the namespace of the UI and server
#' @param signalData - Information on peaks selected from sliding windown
#' analysis in the previous page of shiny app.
#' @param Autotuner - Data with specifics on samples path and metadata.
#'
#' @export
## outer peakvis server
peakVis <- function(input, output, session, signalData, Autotuner) {

    ns <- session$ns


    # checks for completion of signal processing steps ------------------------
    # checks to make sure signal processing is run before peak visualization
    check <- reactiveValues(miss = 0,
                            paramEst = F,
                            params = list())

    observe({
        req(signalData())
        check$miss <- 1
    })

    # all UI functions --------------------------------------------------------
    # Numeric input for max number of peaks to check
    output$peakNumber <- renderUI({

        validate(
            need(check$miss > 0,
               "Please accept the signal processing parameters on the 'Find Peaks' page.")
        )

        tagList(
            h3("You have accepted peak identifimessageion parameters!"),
            br(),
            numericInput(ns("peaks"),
                       label = "Select Number of Peaks to Search between Samples:",
                       value = 10,
                       min = 5,
                       max = 15)
        )
    })

    # UIs associated with the plot visualization
    output$plotCode <- renderUI({

        req(peakView)

        selectInput(ns("plot"),
                    "Choose a peak found within atleast 2 samples:",
                    1:length(totalPeaks()))
    })

    ## UI for picking the peak boundaries to be plotted
    output$bound <- renderUI({

        req(peakView, totalPeaks, check$miss > 0)
        numericInput(ns("boundary"),
                     "Select a boundary window to plot peak:",
                     value = 10,
                     min = 1,
                     step = 1)
    })

    ## UI for generous mz threshold
    output$thresh <- renderUI({

        req(peakView, totalPeaks, check$miss > 0)
        numericInput(ns("mzThreshold"),
                     "Input a generous mz threshold (min = 0.005, max = 0.1):",
                     min = 0.005,
                     max = 0.1,
                     value = 0.01)
    })

    ## UI to continue to next step of autotuner
    output$button <- renderUI({

        req(peakView, totalPeaks, check$miss > 0)
        actionButton(ns("acceptPeaks"),
                     label = "Accept Observed Peaks")
    })

    ## UI to plot peaks
    output$renderedPlot <- renderPlot({

        req(input$boundary, input$plot, peakView())

        peak_diff <- peakView()$peak_diff
        peak_table <- peakView()$peak_table

        plot_peaks(Autotuner,
                   boundary = input$boundary,
                   peak = input$plot,
                   peak_difference = peak_diff,
                   peak_table = peak_table)

    })

    output$peakPlot <- renderUI({
        plotOutput(ns("renderedPlot"))
    })


    # Calculating things to plot ---------------------------------------------
    # Functions that identify, extract, infer, and match peaks across samples
    peakView <- reactive({

        req(input$peaks)
        req(signalData())

        ## the extracted peak list isn't used further beyond this reactive call
        peakList <- extract_peaks(Autotuner,
                                  returned_peaks = input$peaks,
                                  signalData())



        peak_table <- peakwidth_table(Autotuner,
                                      peakList,
                                      returned_peaks = input$peaks)
        rm(peakList)
        peak_diff <- peak_time_difference(peak_table)

        return(list(peak_table = peak_table,
                    peak_diff = peak_diff))

    })




    totalPeaks <- reactive({

        req(peakView())
        temp <- unlist(peakView()$peak_diff$index)
        unique(temp)

    })

    # Actual algorithms that calculate parameters -----------------------------
    observeEvent(input$acceptPeaks, {

        message("~~~ Starting Autotuner Algorithm ~~~\n")
        # Organizing input for estimates ------------------------------------------
        peak_table <- peakView()[[1]]
        peak_diff <- peakView()[[2]]


        # Estimating TIC derived Parameters ---------------------------------------
        message("~~~ Estimating TIC Parameters ~~~\n")
        ticParams <- TIC_params(peak_table = peak_table,
                             peak_difference = peak_diff)
        message("~~~ TIC Parameter Estimation Complete ~~~\n")
        ticParams <- unlist(ticParams)



        # Estimating EIC derived Parameters ---------------------------------------
        ## estimating EIC parameters
        message("~~~ Estimating EIC Parameters ~~~\n")
        eicParamEsts <- EICparams(Autotuner,
                                  massThresh = input$mzThreshold,
                                  peak_table = peak_table)

        message("~~~ EIC Parameter Estimation Complete ~~~\n")

        colNameCheck <- all(c("ppm", "peakCount", "noiseThreshold", "prefilterI", "prefilterScan",
                              "TenPercentQuanSN","maxPw", "minPw") %in% colnames(eicParamEsts))

        assertthat::assert_that(colNameCheck,
                                msg = "Error in EICparams: some of the column names from the output are missing. Cannot complete parameter estimation.")

        ppmEst <- weighted.mean(eicParamEsts$ppm, eicParamEsts$peakCount)
        noiseEst <- min(eicParamEsts$noiseThreshold, na.rm = T)
        prefilterIEst <- min(eicParamEsts$prefilterI, na.rm = T)
        prefilterScanEst <- min(eicParamEsts$prefilterScan, na.rm = T)
        snEst <- min(eicParamEsts$TenPercentQuanSN, na.rm = T)

        ## added this heuristic when poor resoluti
        if(any(ticParams[1]*2 < eicParamEsts$maxPw)) {
            maxPw <- prod(eicParamEsts$maxPw)^(1/length(eicParamEsts$maxPw))
        } else {
            maxPw <- max(eicParamEsts$maxPw)
        }

        minPw <- min(eicParamEsts$minPw)

        estimates <- c(ppm = ppmEst,
                       noise = noiseEst,
                       preIntensity = prefilterIEst,
                       preScan = prefilterScanEst,
                       snThresh = snEst,
                       "Max Peakwidth" = maxPw,
                       "Min Peakwidth" = minPw)

        rm(ppmEst,noiseEst,prefilterIEst,prefilterScanEst,snEst,maxPw,minPw)

        ppmSd <- sd(eicParamEsts$ppm)
        noiseSd <- sd(eicParamEsts$noiseThreshold)
        prefilSd <- sd(eicParamEsts$prefilterI, na.rm = T)
        prefilScanSd <- sd(eicParamEsts$prefilterScan, na.rm = T)
        snEstSd <- sd(eicParamEsts$TenPercentQuanSN, na.rm = T)
        maxPwDist <- abs(diff(sort(eicParamEsts$maxPw, decreasing = T))[1])
        minPwDist <- abs(diff(sort(eicParamEsts$minPw, decreasing = T))[2])

        variability <- c(ppmSd, noiseSd, prefilSd, prefilScanSd, snEstSd,
                         maxPwDist,
                         minPwDist)
        rm(ppmSd, noiseSd, prefilSd, prefilScanSd, snEstSd, maxPwDist, minPwDist)
        description <- c("Standard Deviation of all PPM Estimates",
                         "Standard Deviation of all noise Estimates",
                         'Standard Deviation of all prefileter Intensity Estimates',
                         "Standard Deviation of all scan coung Estimates",
                         "Standard Deviation of all s/n threshold Estimates",
                         "Distance between two highest estimated peak widths",
                         "Distance between two lowest estimated peak widths")


        aggregatedEstimates <- data.frame(Parameters = names(estimates),
                                          estimates = signif(estimates, digits = 4),
                                          'Variability Measure' = signif(variability, digits = 4),
                                          "Measure" = description)




        # Returning Objects -------------------------------------------------------
        check$params <- list(aggregatedEstimates, ticParams)

    }) ## end of param estimation observer

    # Visual confirmation of end of algo --------------------------------------
    startUI <- eventReactive(input$acceptPeaks, {
        TRUE
    })

    observe({
        req(length(check$params) > 0)
        check$paramEst <- TRUE
    })

    output$status <- renderUI({

        req(startUI())
        if(check$paramEst == TRUE) {
            h3("Parameter Estimation is Complete!")
        } else {
            h3("Parameters are Being Estimated!")
        }
    })

    paramsOut <- reactive({
        req(length(check$params) > 0)
        check$params
    })

    return(paramsOut)
}

