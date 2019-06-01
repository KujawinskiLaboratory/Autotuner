# Creating Peak Matrix from Scan Data -------------------------------------
#' @title dissectScans
#'
#' @description This function is designed to extract all MS1 scan features
#' observed within the bounds of the current TIC peak.
#'
#' @param mzDb - This is a list of two column data frames containing information
#' on each mass spectra within the data.
#' @param observedPeak - A list with 'start' and 'stop' boundaries of the
#' current peak.
#' @param header - This is the header file containing all the metadata
#' for the currently loaded sample.
#'
#' @export
dissectScans <- function(mzDb, observedPeak, header) {


    scansOfPeak <- which(observedPeak$start < header$retentionTime &
                             header$retentionTime  < observedPeak$end)
    peakHead <- header[scansOfPeak,]
    ms1 <- peakHead$msLevel == 1L

    # Added this on 2019-03-24 for cases where ms2 data is not within the
    # ms convert file
    if(!all(header$msLevel == 1L)) {
        scansOfPeak <- as.numeric(sub(".* scan=", "", peakHead$spectrumId[ms1]))
    }

    rm(peakHead,ms1)

    peakMassSpectras <- mzDb[scansOfPeak]
    for(i in 1:length(scansOfPeak)) {
        peakMassSpectras[[i]] <- cbind(peakMassSpectras[[i]], scansOfPeak[i])
    }
    peakMassSpectras <- Reduce(rbind, peakMassSpectras)
    colnames(peakMassSpectras) <- c("mz", "intensity", "scan")
    peakMassSpectras <- data.frame(peakMassSpectras)

    peakMassSpectras <- peakMassSpectras[order(peakMassSpectras$mz),]
    peakMassSpectras <- peakMassSpectras[peakMassSpectras$intensity > 0,]

    return(peakMassSpectras)

}
