#' @title checkBounds
#'
#' @description Recursive function used to find how far a binned feature might
#' extend beyond the boundary of the originally defined TIC peak.
#'
#' @param mass Specific mass being checked against adjacent scans.
#' @param upper A boolean value that tells the algorithm to check indices
#' greater than the entered one.  If false, it will check values less thant the
#' entered one.
#' @param mzDb A list of data.frames containing the m/z and intensity values
#' from each scan's mass spectra.
#' @param currentIndex Numerical index indicating which scan contains
#' feature specific information.
#' @param intensityStorage A vector used during recursion to store intensity
#' values as they are added to peak expansion.
#' @param ppmEst Scalar numerical value meant to represent the ppm of the
#' instrument.
#' @param origBound The original scan bound location of the peak.
#' @param scans Set of all possible ms1 scans for the sample.
#' @param header A data.fame containing metadata on the sample like
#' spectra type (MS1 vs MS2), retention time, and scan count.
#'
#' @import Biobase methods
#'
#' @return This function returns the last index the feature is detected.
checkBounds <- function(mass,
                        upper = TRUE,
                        mzDb,
                        currentIndex,
                        intensityStorage,
                        ppmEst,
                        scans,
                        origBound,
                        header) {


    # Check to make sure we havent reached the boundary -----------------------
    # introducing check to make sure this doesn't run forever
    runawayPeak <- FALSE
    if(upper) {
        toowide <- currentIndex  >  origBound + 500
        edge <- max(scans) == currentIndex
        if(toowide || edge) {
            runawayPeak <- TRUE
        }

    } else {

        toowide <- currentIndex < origBound - 500
        edge <- min(scans) == currentIndex
        if(toowide || edge) {
            runawayPeak <- TRUE
        }
    }

    if(runawayPeak) {
        return("Runaway Peak")
    }

    # initializing storage objects --------------------------------------------
    # intensitty storage
    if(missing(intensityStorage)) {
        intensityStorage <- vector(mode = "numeric")
    }

    # next scan

    scanIndex <- which(scans %in% currentIndex)
    if(upper) {
        adjIndex <- scanIndex+1
    } else {
        adjIndex <- scanIndex-1
    }
    rm(scanIndex)

    nextIndex <- scans[adjIndex]

    if(is.na(nextIndex)) {
        ## hack for netCDF files
        nextIndex <- suppressWarnings(as.numeric(
            sub("(.* )?scan=|(.* )?scanId=",
                "",header$spectrumId[adjIndex])))
    }

    peakMatrix <- data.frame(mzDb[[adjIndex]])
    if(ncol(peakMatrix) == 0) {
        return(0)
    }
    colnames(peakMatrix) <- c("mz", "intensity")


    # checking if there is a match in the boundary ----------------------------
    massSpectraMatch <- vector(mode = "numeric", length = length(peakMatrix$mz))
    for(i in seq_along(peakMatrix$mz)) {
        massSpectraMatch[i] <- estimatePPM(first = peakMatrix$mz[i],
                                            second = mass)
    }
    massSpectraMatch <- massSpectraMatch < ppmEst

    foundMass <- any(massSpectraMatch)

    # updating the function to check next scan --------------------------------
    if(foundMass) {

        # adding intensity to storage variable --------------------------------
        # correcting for multiple matches
        if(length(foundMass) > 1) {
            minMz <- which.min(peakMatrix$mz - peakMatrix$mz[massSpectraMatch])
            newBound <- peakMatrix$intensity[massSpectraMatch][minMz]
            intensityStorage[length(intensityStorage) + 1] <- newBound
        } else {
            intensityStorage[length(intensityStorage) + 1] <-
                peakMatrix$intensity[massSpectraMatch][1]
        }


        # correlating peak shape  ---------------------------------------------
        ## checking if peak is increasing or decreasing monotonically
        if(length(intensityStorage) >= 3) {

            fit <- lm(intensityStorage ~ seq_along(intensityStorage))
            slope <- stats::coef(fit)[2]
            r2 <- suppressWarnings(summary(fit)$r.squared)

            if(abs(slope) < 0.75 | r2 < .75) {
                return(currentIndex)
            }

        }

        rm(peakMatrix)

        bound <- checkBounds(mass = mass,
                            upper = upper,
                            mzDb = mzDb,
                            currentIndex = nextIndex,
                            intensityStorage = intensityStorage,
                            ppmEst = ppmEst,
                            scans = scans,
                            origBound = origBound,
                            header = header)
        return(bound)

    } else {

        if(length(intensityStorage) == 0) {
            intensityStorage <- NA
        }

        rm(peakMatrix)

        return(currentIndex)

    }

}
