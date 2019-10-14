#' @title dissectScans
#'
#' @description This function is designed to extract all MS1 scan features
#' observed within the bounds of the current TIC peak.
#'
#' @param mzDb This is a list of two column data frames containing information
#' on each mass spectra within the data.
#' @param observedPeak A list with 'start' and 'stop' boundaries of the
#' current peak.
#' @param header This is the header file containing all the metadata
#' for the currently loaded sample.
#'
#' @return A Peak Matrix from Scan Data
dissectScans <- function(mzDb, observedPeak, header) {


    scansOfPeak <- which(observedPeak$start < header$retentionTime &
                             header$retentionTime  < observedPeak$end)
    peakHead <- header[scansOfPeak,]
    ms1 <- peakHead$msLevel == 1L

    ## removing if statement here since everything is entered as MS1 spectra
    ## from the get go 2019-06-19
    scanID <- as.numeric(sub("(.* )?scan=|(.* )?scanId=", "", peakHead$spectrumId[ms1]))

    rm(peakHead,ms1)

    peakMassSpectras <- mzDb[scansOfPeak]
    for(i in seq_along(scansOfPeak)) {
        peakMassSpectras[[i]] <- cbind(peakMassSpectras[[i]],
                                       scansOfPeak[i],
                                       scanID[i])
    }
    peakMassSpectras <- Reduce(rbind, peakMassSpectras)
    colnames(peakMassSpectras) <- c("mz", "intensity", "scan", "scanID")
    peakMassSpectras <- data.frame(peakMassSpectras)

    peakMassSpectras <- peakMassSpectras[order(peakMassSpectras$mz),]
    peakMassSpectras <- peakMassSpectras[peakMassSpectras$intensity > 0,]

    return(peakMassSpectras)

}
