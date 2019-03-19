#' @title loadDataUI
#'
#' @description This function contains the UI component of the data loading
#' portion of the shiny app used within Autotuner.
#'
#' @param id - character string used to map UI and server within App'
#'
#' @export

# Load data UI ------------------------------------------------------------

loadDataUI <- function(id) {
    ns <- NS(id)
    tagList(
        titlePanel("Load Data into Autotuner"),

# Left Hand Column --------------------------------------------------------
        column(3,
            wellPanel(

                 h3("Load Metadata File and Raw Data Files"),
                 shiny::fileInput(ns('datafile'),
                                  'Choose Metadata CSV file:',
                           accept=c('text/csv',
                                    'text/comma-separated-values,text/plain'),
                           placeholder = ""),
                 br(),

                 shiny::uiOutput(ns("sampleCols")),
                 br(),

                 h5("Choose Directory Containing Data:"),
                 shinyFiles::shinyDirButton(ns("dir"),
                                            label = "Go!",
                                            title = "Upload",
                                            buttonType = "default"),

                 uiOutput(ns("subsetMeta"))
        )
    ),

# Middle Column -----------------------------------------------------------
    column(6,
        wellPanel(
            h3("Entered Metadata"),
            DT::dataTableOutput(ns("table"))
        )
    ),
# Last Column -------------------------------------------------------------
    column(3,
        wellPanel(
            h3("Enter Information About Data Classes"),
            uiOutput(ns("classFactor")),
            uiOutput(ns("selectedSamples")),
            uiOutput(ns("nextStep")),
            uiOutput(ns("status"))
        )
    )
)

}

# Server Function for load data -------------------------------------------

#' @title loadData
#'
#' @description This function contains the server component of the data loading
#' portion of the shiny app used within Autotuner.
#'
#' @param input - List of UI entered values.
#' @param output - Values exported from server to the UI after computation.
#' @param session - String used to match the namespace of the UI and server
#'
#' @export

#This function is repsonsible for loading in the selected file
loadData <- function(input, output, session) {

    ns <- session$ns
    # Description of Stored Data ----------------------------------------------
    ## Data is ultimately what will leave this app:
    ## runfile - entire runfile subsetted by all samples w/in folder
    ## runfileSub - subset of runfile with files selected for visual analysis
    ## path - selected path to data folder saved as a string
    ## fileNames - string vector of all matching files within data folder
    ## files - path to each individual file selected
    ## filesSub - path to the subset of selected files for visual analysis
    ## Autotuner - Autotuner object loaded for further processing

    Data <- reactiveValues(runfile = NA,
                 runfileSub = NA,
                 runFilePath = "",
                 files = "",
                 filesSub = "",
                 path = "",
                 ready = FALSE,
                 fileNames = NA,
                 Autotuner = NA,
                 sampleID = "")

    # loading metadata --------------------------------------------------------
    # loading path to metadata
    observe({
        req(input$datafile != "")
        Data$runFilePath <- input$datafile
    })

    # reading in metadata into Data reactive value object
    observe({

        req(Data$runFilePath != "")
        infile <- Data$runFilePath

        # User has not uploaded a file yet
        if (is.null(infile)) {
            return(NULL)
        }

        Data$runfile <- utils::read.csv(infile$datapath, stringsAsFactors = F)
    })

    # generating UI for metadata spreadsheet
    output$table <- DT::renderDataTable({
        validate(need(!is.na(Data$runfile), "Enter a Metadata File"))
        validate(need(grepl("\\.csv", Data$runFilePath), "Metadata File Must be .csv File"))
        Data$runfile
    },
    options = list(scrollX = TRUE),
    rownames= FALSE)

    # UI for selecting column within metadata that has file names
    output$sampleCols <- renderUI({
        req(!is.na(Data$runfile))
        selectOptions <- colnames(Data$runfile)
        selectizeInput(ns("colChoice"),
            choices = selectOptions,
            label = "Select the column within the metadata cointaing file names.")
    })

    # loading in path to data -------------------------------------------------

    # shinyFiles linker function to search for local directory
    shinyFiles::shinyDirChoose(input,
                               id = 'dir',
                               session = session,
                               roots = c(home = path.expand("~")))

    # filling in data object with information on sample location
    observe({

        path <- shinyFiles::parseDirPath(roots = c(home = normalizePath("~/")),
                                         selection = input$dir)
        req(length(path) > 0)
        Data$path <- path
        dataFiles <- list.files(Data$path)
        Data$fileNames <- dataFiles
        Data$files <- file.path(Data$path, dataFiles)
    })

    # subsetting metadata with respect to loaded files ------------------------

    # UI narrowing down all metadata into a subset based on matched samples
    output$subsetMeta <- renderUI({

        req(Data$runFilePath != "")
        req(Data$files != "")
        req(input$colChoice)

        supportedFileTypes <- c("mzXML","mzML","mzData","NetCDF", "CDF")

        checkFileType <- vector(mode = "logical",
                                length = length(supportedFileTypes))
        for(fileTypeIndex in 1:length(supportedFileTypes)) {
            fileType <- supportedFileTypes[fileTypeIndex]
            checkFileType[fileTypeIndex] <- all(grepl(fileType, Data$files))
        }

        validate(
            need(any(checkFileType),
                 "This software curently only supports mzXML, mzML, mzData, or NetCDF files. The directory selected must only contain files of this type.")
        )

        tagList(
            h5("Subset Metadata by files:"),
            shiny::actionButton(ns("subsetMeta"), "Go!")
        )

    })

    # server narrowing down all metadata into a subset based on matched samples
    observeEvent(input$subsetMeta, {

        runfilePaths <- file.path(Data$path, Data$runfile[,input$colChoice])

        if(any(grepl(".mzML|.mzXML|.mzData|.NetCDF|.CDF",runfilePaths))) {
            runfilePaths <- sub(".mzML|.mzXML|.mzData|.NetCDF|.CDF", "",
                                runfilePaths)
        }

        selectedFiles <- basename(Data$files)

        noSuffix <- gsub(pattern = c(".mzML|.mzXML|.mzData|.NetCDF|.CDF"),
                 replacement = "", selectedFiles)

        runfile_subset <- which(basename(runfilePaths) %in% noSuffix)
        Data$runfile <- Data$runfile[runfile_subset,]

    })

    # Pick Class Factor Column ------------------------------------------------
    # selection UI to pick factor column within metadata
    output$classFactor <- renderUI({

        validate(need(input$subsetMeta > 0, "Please make sure to first subset the metadata with specific files being processed."))
        selectOptions <- colnames(Data$runfile)
        selectizeInput(ns("classFactor"),
               choices = selectOptions,
               label = "Pick the column in the metadata containing class specific information for sample classes:",
               selected = input$colChoice)

    })

    # subsetting data based on selected files  --------------------------------
    # replacing common string between sample to make selector input more
    # legible
    observe({

        req(!is.na(Data$fileNames))

        storeLCS <- vector()
        if(length(Data$fileNames) > 20) {
            randStrings <- base::sample(x = seq_along(Data$fileNames),
                                        size = 20)
        } else {
            randStrings <- base::sample(x = seq_along(Data$fileNames),
                                        size = length(Data$fileNames))
        }

        randOne <- randStrings[c(TRUE, FALSE)]
        randTwo <- randStrings[c(FALSE, TRUE)]

        for(rand in 1:min(length(randOne), length(randTwo))) {

            string1 <- randOne[rand]
            string2 <- randTwo[rand]

            lscOut <- PTXQC::LCS(Data$fileNames[string1],
                                 Data$fileNames[string2])

            # checking wether to store results
            if(rand == 1) {
                storeLCS <- lscOut
            } else if(nchar(lscOut) < nchar(storeLCS)) {
                storeLCS <- lscOut
            }

        }

        temp <- sub(storeLCS, "Sample_", Data$fileNames)
        Data$sampleID <- sub("\\..*", "", temp)

    })

    # UI for sample selection input
    output$selectedSamples <- renderUI({

        req(input$classFactor)
        req(Data$sampleID != "")
        validate(need(input$classFactor != input$colChoice,
            "Pick a column other than the sample name column."))

        shiny::selectizeInput(ns("selectedSamples"),
                              choices = Data$sampleID,
                              label = "Select samples to process further:",
                              selected = input$selectedSamples,
                              multiple = TRUE)

    })

    # loadind subsetted metadata file based on selected samples to data object
    observeEvent(input$nextStep ,{

        req(length(input$selectedSamples) > 2)

        filesSelected <- which(Data$sampleID %in% input$selectedSamples)
        data <- Data$runfile
        fileCol <- input$colChoice
        Data$runfileSub <- data[filesSelected,]

    })

    # UI button to tell R to load selected samples into an Autotuner object
    output$nextStep <- renderUI({

        req(input$classFactor != input$colChoice)
        validate(need(length(input$selectedSamples) > 2, "Load at least 3 samples."))
        tagList(
            h5("Tune XCMS with sample data subset:"),
            actionButton(ns("nextStep"), "Go!")
        )

    })

    # create Autotuner obj ---------------------------------------------------
    # filling in the Autotuner slot
    observe({

        req(!is.na(Data$runfileSub))
        filesSelected <- which(Data$sampleID %in% input$selectedSamples)
        Data$filesSub <- file.path(Data$path,
                                   Data$fileNames[filesSelected])
        Data$Autotuner <- createAutotuner(data_paths = Data$filesSub,
                                          runfile = Data$runfileSub,
                                          file_col = input$colChoice,
                                          factorCol = input$classFactor)
        Data$ready = TRUE

    })

    # return object  ----------------------------------------------------------
    # Visual check to let users known they are ready to go
    output$status <- renderUI({

        req(input$subsetMeta > 0)
        req(Data$filesSub != "")
        if(Data$ready == TRUE) {
            h3("You are good to go!")
        }
    })

    # Object leaving the function ---------------------------------------------
    # creating output object for the function
    loadReturn <- reactive({

        req(Data$ready == TRUE)
        enteredData <- list(
            Autotuner = Data$Autotuner,
            runfile = Data$runfile,
            dataPaths = Data$files)
        return(enteredData)
    })

    return(loadReturn)

}

