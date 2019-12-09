#' Autotuner
#'
#' @description This file contains the skeleton to the Autotuner class used
#' through out Autoutuner.
#'
#' @description This object is a generic object designed to run the different
#' functions of the Autotuner package. The slots represent content or data
#' that the package uses throughout the different functions.
#'
#' @slot time A list containing vectors of scan time points
#' from each sample.
#' @slot intensity A list containing vectors of scan intensity
#' points from each sample.
#' @slot peaks Regions within each sample identified as
#' peaks by slidingwindow analysis.
#' @slot peak_table A data.frame containing information on
#' each peak after further processing is done to the data.
#' @slot peak_difference A data.frame containing information
#' on how peaks are eluted differently over time.
#' @slot metadata A data.frame containing metadata for all
#' samples to be run on Autotuner.
#' @slot file_paths A string path that leads to the samples
#' to be run on Autotuner.
#' @slot file_col A string for the column name of the
#' column within the metadata that has specific sample names.
#' @slot factorCol A string for the column name of the column
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
#' @param .Object The \code{\linkS4class{Autotuner}} object to update.
#' @param data_paths A string path pointing at data files to
#' load in Autotuner.
#' @param runfile A data.frame of sample metadata.
#' @param file_col Character string of the column name of the
#' column within the runfile that contains sample names.
#' @param factorCol Character string of the column name of the
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
                for(index in seq_len(nrow(runfile))) {

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
                        for(i in seq_along(allInts)) {
                            storeInt[[i]] <- sum(unlist(allInts[i]))
                        }
                        intensity[[index]] <- unlist(storeInt)
                        rm(allInts, storeInt)
                    }

                }

                message("~~~ Storing Everything in Autotuner Object ~~~ \n")
                .Object <- setAutoTime(time = time, Autotuner = .Object)
                .Object <- setAutoIntensity(intensity = intensity,
                                            Autotuner = .Object)
                .Object <- setAutoMetadata(metadata = runfile,
                                                    Autotuner = .Object)
                .Object <- setAutoFile_paths(file_paths = data_paths,
                                                        Autotuner = .Object)
                .Object <- setAutoFile_col(file_col = file_col,
                                                    Autotuner = .Object)
                .Object <- setAutoFactorCol(factorCol = factorCol,
                                            Autotuner = .Object)

                message("~~~ The Autotuner Object has been Created ~~~ \n")
                return(.Object)

})

#' @title createAutotuner
#'
#' @description This function will create a Autotuner used to extract ms2s.
#'
#' @param data_paths A string path pointing at data files to
#' load in Autotuner.
#' @param runfile A data.frame of sample metadata.
#' @param file_col Character string of the column name of the
#' column within the runfile that contains sample names.
#' @param factorCol Character string of the column name of the
#' column within the runfile that contains sample type factor.
#'
#' @export
#'
#' @examples
#' library(mtbls2)
#' rawPaths <- c(
#' system.file("mzData/MSpos-Ex2-cyp79-48h-Ag-1_1-B,3_01_9828.mzData",
#' package = "mtbls2"),
#' system.file("mzData/MSpos-Ex2-cyp79-48h-Ag-2_1-B,4_01_9830.mzData",
#' package = "mtbls2"),
#' system.file("mzData/MSpos-Ex2-cyp79-48h-Ag-4_1-B,4_01_9834.mzData",
#' package = "mtbls2"))
#'
#' metadata <- read.table(system.file(
#' "a_mtbl2_metabolite_profiling_mass_spectrometry.txt", package = "mtbls2"),
#' header = TRUE, stringsAsFactors = FALSE)
#' metadata <- metadata[sub("mzData/", "",
#' metadata$Raw.Spectral.Data.File) %in% basename(rawPaths),]
#'
#' Autotuner <- Autotuner::createAutotuner(rawPaths,
#' metadata,
#' file_col = "Raw.Spectral.Data.File",
#' factorCol = "Factor.Value.genotype.")
#'
#'
#' @return This function returns an Autotuner object
createAutotuner <- function(data_paths, runfile, file_col, factorCol) {

    if(nrow(runfile) == 0) {
        stop(paste('The entered metadata file does not contain any rows.',
                    'Check the input.'))
    }

    if(!all(file.exists(data_paths))) {
        stop("One or more of the file paths do not exist.")
    }

    if(nrow(runfile) != length(data_paths)) {
        stop("Number of file paths and metadata entries do not match. Check input.")
    }

    Autotuner <- methods::new(Class="Autotuner", data_paths,
                                runfile,
                                file_col,
                                factorCol)



    return(Autotuner)
}


# accesor methods for slots -----------------------------------------------


#' @title getAutoIntensity
#'
#' @description This function is designed to return the list intensities
#' obtained from applying the sliding window analysis to the raw data stored
#' within an AutoTuner object.
#'
#' @param Autotuner An AutoTuner object.
#'
#' @return The content of the intensity slot
#'
#' @export
#'
#' @examples
#' Autotuner <- readRDS(system.file("extdata/Autotuner.rds",
#' package="Autotuner"))
#' intensity <- getAutoIntensity(Autotuner)
getAutoIntensity <- function(Autotuner) {

    return(Autotuner@intensity)

}


#' @title getAutoTime
#'
#' @description This function returns the list of numerics stored within
#' the 'time' slot of the Autotuner Object
#'
#' @param Autotuner An AutoTuner object.
#'
#' @return The content of the time slot
#'
#' @export
#'
#' @examples
#' Autotuner <- readRDS(system.file("extdata/Autotuner.rds",
#' package="Autotuner"))
#' intensity <- getAutoTime(Autotuner)
getAutoTime <- function(Autotuner) {

    return(Autotuner@time)

}


#' @title getAutoPeaks
#'
#' @description This function returns the list of numerics stored within
#' the 'peaks' slot of the Autotuner Object
#'
#' @param Autotuner An AutoTuner object.
#'
#' @return The content of the peaks slot
#'
#' @export
#'
#' @examples
#' Autotuner <- readRDS(system.file("extdata/Autotuner.rds",
#' package="Autotuner"))
#' intensity <- getAutoPeaks(Autotuner)
getAutoPeaks <- function(Autotuner) {

    return(Autotuner@peaks)

}


#' @title getAutoPeak_table
#'
#' @description This function returns the data.frame stored within
#' the 'peak_table' slot of the Autotuner Object
#'
#' @param Autotuner An AutoTuner object.
#'
#' @return The content of the peak_table slot
#'
#' @export
#'
#' @examples
#' Autotuner <- readRDS(system.file("extdata/Autotuner.rds",
#' package="Autotuner"))
#' intensity <- getAutoPeak_table(Autotuner)
getAutoPeak_table <- function(Autotuner) {

    return(Autotuner@peak_table)

}

#' @title getAutoPeak_difference
#'
#' @description This function returns the data.frame stored within
#' the 'peak_difference' slot of the Autotuner Object
#'
#' @param Autotuner An AutoTuner object.
#'
#' @return The content of the peak_difference slot
#'
#' @export
#'
#' @examples
#' Autotuner <- readRDS(system.file("extdata/Autotuner.rds",
#' package="Autotuner"))
#' intensity <- getAutoPeak_difference(Autotuner)
getAutoPeak_difference <- function(Autotuner) {

    return(Autotuner@peak_difference)

}


#' @title getAutoMetadata
#'
#' @description This function returns the data.frame stored within
#' the 'meatadata' slot of the Autotuner Object
#'
#' @param Autotuner An AutoTuner object.
#'
#' @return The content of the meatadata slot
#'
#' @export
#'
#' @examples
#' Autotuner <- readRDS(system.file("extdata/Autotuner.rds",
#' package="Autotuner"))
#' intensity <- getAutoMetadata(Autotuner)
getAutoMetadata <- function(Autotuner) {

    return(Autotuner@metadata)

}

#' @title getAutoFile_paths
#'
#' @description This function returns the character string stored within
#' the 'file_paths' slot of the Autotuner Object
#'
#' @param Autotuner An AutoTuner object.
#'
#' @return The content of the file_paths slot
#'
#' @export
#'
#' @examples
#' Autotuner <- readRDS(system.file("extdata/Autotuner.rds",
#' package="Autotuner"))
#' intensity <- getAutoFile_paths(Autotuner)
getAutoFile_paths <- function(Autotuner) {

    return(Autotuner@file_paths)

}


#' @title getAutoFile_col
#'
#' @description This function returns the character string stored within
#' the 'file_col' slot of the Autotuner Object
#'
#' @param Autotuner An AutoTuner object.
#'
#' @return The content of the file_col slot
#'
#' @export
#'
#' @examples
#' Autotuner <- readRDS(system.file("extdata/Autotuner.rds",
#' package="Autotuner"))
#' intensity <- getAutoFile_col(Autotuner)
getAutoFile_col <- function(Autotuner) {

    return(Autotuner@file_col)

}

#' @title getAutoFactorCol
#'
#' @description This function returns the character string stored within
#' the 'factorCol' slot of the Autotuner Object
#'
#' @param Autotuner An AutoTuner object.
#'
#' @return The content of the factorCol slot
#'
#' @export
#'
#' @examples
#' Autotuner <- readRDS(system.file("extdata/Autotuner.rds",
#' package="Autotuner"))
#' intensity <- getAutoFactorCol(Autotuner)
getAutoFactorCol <- function(Autotuner) {

    return(Autotuner@factorCol)

}



# Set methods for slots ---------------------------------------------------

#' @title setAutoIntensity
#'
#' @description This function fills the "intensity" slot within an Autotuner
#' object.
#'
#' @param Autotuner An AutoTuner object.
#' @param intensity A list of numeric values representing intensity
#'
#' @return An Autotuner object with a filled intensity slot
#'
#' @export
#'
#' @examples
#' Autotuner <- readRDS(system.file("extdata/Autotuner.rds",
#' package="Autotuner"))
#' intensity <- getAutoIntensity(Autotuner)
#' Autotuner <- setAutoIntensity(intensity, Autotuner)
setAutoIntensity <- function(intensity, Autotuner) {

    Autotuner@intensity <- intensity
    return(Autotuner)

}


#' @title setAutoTime
#'
#' @description This function fills the "time" slot within an Autotuner object.
#'
#' @param Autotuner An AutoTuner object.
#' @param time A list of numeric values representing time
#'
#' @return An Autotuner object with a filled time slot
#'
#' @export
#'
#' @examples
#' Autotuner <- readRDS(system.file("extdata/Autotuner.rds",
#' package="Autotuner"))
#' time <- getAutoTime(Autotuner)
#' Autotuner <- setAutoTime(time, Autotuner)
setAutoTime <- function(time, Autotuner) {

    Autotuner@time <- time
    return(Autotuner)

}


#' @title setAutoPeaks
#'
#' @description This function fills the peaks slot within an Autotuner object.
#'
#' @param Autotuner An AutoTuner object.
#' @param peaks A list of numeric values representing peaks
#'
#' @return An Autotuner object with a filled peaks slot
#'
#' @export
#'
#' @examples
#' Autotuner <- readRDS(system.file("extdata/Autotuner.rds",
#' package="Autotuner"))
#' peaks <- getAutoPeaks(Autotuner)
#' Autotuner <- setAutoPeaks(peaks, Autotuner)
setAutoPeaks <- function(peaks, Autotuner) {

    Autotuner@peaks <- peaks
    return(Autotuner)

}


#' @title setAutoPeak_table
#'
#' @description This function fills the peak_table slot within an Autotuner
#' object.
#'
#' @param Autotuner An AutoTuner object.
#' @param peak_table A data.frame representing peak_table
#'
#' @return An Autotuner object with a filled peak_table slot
#'
#' @export
#'
#' @examples
#' Autotuner <- readRDS(system.file("extdata/Autotuner.rds",
#' package="Autotuner"))
#' peak_table <- getAutoPeak_table(Autotuner)
#' Autotuner <- setAutoPeak_table(peak_table, Autotuner)
setAutoPeak_table <- function(peak_table, Autotuner) {

    Autotuner@peak_table <- peak_table
    return(Autotuner)

}

#' @title setAutoPeak_difference
#'
#' @description This function fills the peak_difference slot within an
#' Autotuner object.
#'
#' @param Autotuner An AutoTuner object.
#' @param peak_difference A data.frame representing peak_difference
#'
#' @return An Autotuner object with a filled peak_difference slot
#'
#' @export
#'
#' @examples
#' Autotuner <- readRDS(system.file("extdata/Autotuner.rds",
#' package="Autotuner"))
#' peak_difference <- getAutoPeak_table(Autotuner)
#' Autotuner <- setAutoPeak_difference(peak_difference, Autotuner)
setAutoPeak_difference <- function(peak_difference, Autotuner) {

    Autotuner@peak_difference <- peak_difference
    return(Autotuner)

}


#' @title setAutoMetadata
#'
#' @description This function fills the metadata slot within an Autotuner
#' object.
#'
#' @param Autotuner An AutoTuner object.
#' @param metadata A data.frame representing metadata
#'
#' @return An Autotuner object with a filled metadata slot
#'
#' @export
#'
#' @examples
#' Autotuner <- readRDS(system.file("extdata/Autotuner.rds",
#' package="Autotuner"))
#' metadata <- getAutoMetadata(Autotuner)
#' Autotuner <- setAutoMetadata(metadata, Autotuner)
setAutoMetadata <- function(metadata, Autotuner) {

    Autotuner@metadata <- metadata
    return(Autotuner)

}

#' @title setAutoFile_paths
#'
#' @description This function fills the file_paths slot within an Autotuner
#' object.
#'
#' @param Autotuner An AutoTuner object.
#' @param file_paths A character vector representing file_paths
#'
#' @return An Autotuner object with a filled file_paths slot
#'
#' @export
#'
#' @examples
#' Autotuner <- readRDS(system.file("extdata/Autotuner.rds",
#' package="Autotuner"))
#' file_paths <- getAutoFile_paths(Autotuner)
#' Autotuner <- setAutoFile_paths(file_paths = file_paths, Autotuner)
setAutoFile_paths <- function(file_paths, Autotuner) {

    Autotuner@file_paths <- file_paths
    return(Autotuner)

}


#' @title setAutoFile_col
#'
#' @description This function fills the file_col slot within an Autotuner
#' object.
#'
#' @param Autotuner An AutoTuner object.
#' @param file_col A character vector representing file_col
#'
#' @return An Autotuner object with a filled file_col slot
#'
#' @export
#'
#' @examples
#' Autotuner <- readRDS(system.file("extdata/Autotuner.rds",
#' package="Autotuner"))
#' file_col <- getAutoFile_col(Autotuner)
#' Autotuner <- setAutoFile_col(file_col, Autotuner)
setAutoFile_col <- function(file_col, Autotuner) {

    Autotuner@file_col <- file_col
    return(Autotuner)

}

#' @title setAutoFactorCol
#'
#' @description This function fills the factorCol slot within an Autotuner
#' object.
#'
#' @param Autotuner An AutoTuner object.
#' @param factorCol A character vector representing factorCol
#'
#' @return An Autotuner object with a filled factorCol slot
#'
#' @export
#'
#' @examples
#' Autotuner <- readRDS(system.file("extdata/Autotuner.rds",
#' package="Autotuner"))
#' factorCol <- getAutoFactorCol(Autotuner)
#' Autotuner <- setAutoFactorCol(factorCol, Autotuner)
setAutoFactorCol <- function(factorCol, Autotuner) {

    Autotuner@factorCol <- factorCol
    return(Autotuner)

}




