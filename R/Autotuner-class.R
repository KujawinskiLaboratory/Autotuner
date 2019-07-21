#' Autotuner
#'
#' @description This file contains the skeleton to the Autotuner class used
#' through out Autoutuner.
#'
#' @description This object is a generic object designed to run the different
#' functions of the ms2sweeper package. The slots represent content or data
#' that the package uses throughout the different functions.
#'
#' @slot time - A list containing vectors of scan time points
#' from each sample.
#' @slot intensity - A list containing vectors of scan intensity
#' points from each sample.
#' @slot peaks - Regions within each sample identified as
#' peaks by slidingwindow analysis.
#' @slot peak_table - A data.frame containing information on
#' each peak after further processing is done to the data.
#' @slot peak_difference - A data.frame containing information
#' on how peaks are eluted differently over time.
#' @slot metadata - A data.frame containing metadata for all
#' samples to be run on Autotuner.
#' @slot file_paths - A string path that leads to the samples
#' to be run on Autotuner.
#' @slot file_col - A string for the column name of the
#' column within the metadata that has specific sample names.
#' @slot factorCol - A string for the column name of the column
#' within the metadata that has specific sample class names.
#'
#' @importFrom MSnbase readMSData
#' @importFrom MSnbase filterFile
#' @importFrom MSnbase rtime
#' @importFrom MSnbase tic
#'
#' @export
Autotuner <- setClass(
    # Set the name for the class
    Class = "Autotuner",

    # Define the slots
    slots = c(
        time = "list",
        intensity = "list",
        peaks = "list",
        peak_table = "data.frame",
        peak_difference = "data.frame",
        metadata = "data.frame",
        file_paths = "character",
        file_col = "character",
        factorCol = "character"),

    prototype = prototype(
        time = list(),
        intensity = list(),
        peaks = list(),
        peak_table = data.frame(),
        peak_difference = data.frame(),
        metadata = data.frame(),
        file_paths = character(),
        file_col = character(),
        factorCol = character())
)


#' @title Update an \code{\linkS4class{Autotuner}} object
#'
#' @description This method updates an \code{\linkS4class{Autotuner}}
#'     object to the latest definition.
#'
#' @param .Object - the \code{\linkS4class{Autotuner}} object to update.
#' @param data_paths - A string path pointing at data files to
#' load in Autotuner.
#' @param runfile - a data.frame of sample metadata.
#' @param file_col - Character string of the column name of the
#' column within the runfile that contains sample names.
#' @param factorCol - Character string of the column name of the
#' column within the runfile that contains sample type factor.
#'
#' @return An updated \code{\linkS4class{Autotuner}} containing all data from
#' the input object.
#'
#' @author Craig McLean
setMethod(f = "initialize", signature = "Autotuner",
            function(.Object, data_paths, runfile, file_col, factorCol) {

                message("~~~ Autotuner: Initializator ~~~ \n")
                message("~~~ Parsing Raw Data into R ~~~ \n")
                raw <- suppressMessages(MSnbase::readMSData(data_paths,
                                                            msLevel. = 1,
                                                            mode = "onDisk"))


                # determining time and intensity data for each sample
                time <- list()
                intensity <- list()
                message(paste("~~~ Extracting the Raw Data from Individual",
                            "Samples ~~~ \n"))
                for(index in 1:nrow(runfile)) {

                    signal_data <- MSnbase::filterFile(raw, file = index)
                    time[[index]] <- MSnbase::rtime(signal_data)
                    intensity[[index]] <- MSnbase::tic(signal_data)

                    ## 2019-07-21 - Adding check for chromatograph data
                    if(sum(intensity[[index]]) == 0) {

                        if(index == 1) {
                            message("No TIC data was found.")
                            message("Using integrated intensity data instead.")
                        }

                        allInts <- MSnbase::intensity(signal_data)
                        storeInt <- list()
                        for(i in 1:length(allInts)) {
                            storeInt[[i]] <- sum(unlist(allInts[i]))
                        }
                        intensity[[index]] <- unlist(storeInt)
                        rm(allInts, storeInt)
                    }

                }

                message("~~~ Storing Everything in Autotuner Object ~~~ \n")
                .Object@time <- time
                .Object@intensity <- intensity
                .Object@metadata <- runfile
                .Object@file_paths <- data_paths
                .Object@file_col <- file_col
                .Object@factorCol <- factorCol

                message("~~~ The Autotuner Object has been Created ~~~ \n")
                return(.Object)

})

#' @title createAutotuner
#'
#' @description This function will create a Autotuner used to extract ms2s.
#'
#' @param data_paths - A string path pointing at data files to
#' load in Autotuner.
#' @param runfile - a data.frame of sample metadata.
#' @param file_col - Character string of the column name of the
#' column within the runfile that contains sample names.
#' @param factorCol - Character string of the column name of the
#' column within the runfile that contains sample type factor.
#'
#' @export
#'
#' @examples
#' library(devtools)
#' if(!require("mmetspData")) {
#'     install_github("crmclean/mmetspData")
#' }
#' library(mmetspData)
#' mmetspFiles <- c(system.file("mzMLs/mtab_mmetsp_ft_120815_24.mzML",
#' package = "mmetspData"), system.file("mzMLs/mtab_mmetsp_ft_120815_25.mzML",
#' package = "mmetspData"), system.file("mzMLs/mtab_mmetsp_ft_120815_26.mzML",
#' package = "mmetspData"))
#' metadata <- read.csv(system.file("mmetsp_metadata.csv",
#' package = "mmetspData"),stringsAsFactors = FALSE)
#' metadata <- metadata[metadata$File.Name %in%
#' sub(pattern = ".mzML", "",basename(mmetspFiles)),]
#' Autotuner <- createAutotuner(mmetspFiles, metadata,
#' file_col = "File.Name", factorCol = "Sample.Type")
#'
#' @return This function returns an Autotuner object
createAutotuner <- function(data_paths, runfile, file_col, factorCol) {

    if(nrow(runfile) == 0) {
        stop(paste('The entered metadata file does not contain any rows.',
                    'Check the input.'))
    }

    Autotuner <- methods::new(Class="Autotuner", data_paths,
                                runfile,
                                file_col,
                                factorCol)



    return(Autotuner)
}

#' @title getAutoIntensity
#'
#' @description This function is designed to return the list intensities
#' obtained from applying the sliding window analysis to the raw data stored
#' within an AutoTuner object.
#'
#' @param Autotuner - An AutoTuner object.
#'
#' @return The content of the intensity slot
#'
#' @export
#'
#' @examples
#' Autotuner <- readRDS(system.file("extdata/Autotuner.rds",
#' package="Autotuner"))
#' getAutoIntensity(Autotuner)
getAutoIntensity <- function(Autotuner) {

    return(Autotuner@intensity)

}

