# Creating Peak Matrix from Scan Data -------------------------------------
#' @title dissectScans
#'
#' @description This function is designed to extract all MS1 scan features
#' observed within the bounds of the current TIC peak.
#'
#' @param currentMsFile - This is an mzR object is linked to the current sample
#' being investigated.
#' @param observedPeak - A list with 'start' and 'stop' boundaries of the
#' current peak.
#' @param sampleChrom - This is the header file containing all the metadata
#' for the currently loaded sample.
#'
#' @export
dissectScans <- function(currentMsFile, observedPeak, sampleChrom) {


    scansOfPeak <- which(observedPeak$start < sampleChrom$retentionTime & sampleChrom$retentionTime  < observedPeak$end)
    peakHead <- sampleChrom[scansOfPeak,]
    ms1 <- peakHead$msLevel == 1L

    # Added this on 2019-03-24 for cases where ms2 data is not within the
    # ms convert file
    if(!all(sampleChrom$msLevel == 1L)) {
        scansOfPeak <- as.numeric(sub(".* scan=", "", peakHead$spectrumId[ms1]))
    }

    rm(peakHead,ms1)

    scanCounter <- 1
    allEIC <- list()
    ## potential solution to speed up computation: save everything to drive
    ## by not using mzR we save space on the drive.

    #cdf <- MSnbase::readMSData(files = data_paths[1], mode = "onDisk")
    #j <- MSnbase::mz(cdf)
    #k <- MSnbase::intensity(cdf)
    #jj <- j[scansOfPeak]
    #kk <- k[scansOfPeak]

    #w <- list()
    #for(i in seq_along(jj)) {
    #    w[[i]] <- cbind(mz = jj[[i]],
    #                    intensity = kk[[i]],
    #                    scan = scansOfPeak[i])
    #}
    #eic2 <- Reduce(rbind, w)


    for(i in 1:length(scansOfPeak)) {

        scan <- scansOfPeak[i]





        peakMatrix <- mzR::peaks(currentMsFile,scan)
        if(is.null(peakMatrix)) {
            next
        }

        colnames(peakMatrix) <- c("mz", "intensity")

        if(nrow(peakMatrix) == 0) {
            next()
        }

        peakStorage <- data.frame(peakMatrix, scanCounter,
                                  scan)
        allEIC[[scanCounter]] <- peakStorage
        scanCounter <- scanCounter + 1
    }

    eicDF <- Reduce(rbind, allEIC)
    returnEIC <- eicDF[order(eicDF$mz),]
    returnEIC <- returnEIC[returnEIC$intensity > 0,]

    return(returnEIC)

}
