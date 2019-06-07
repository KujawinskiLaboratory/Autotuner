#' @title checkEICPeaks
#'
#' @description This function is the outer most function used to check for
#' individual EIC peak specific parameters.
#'
#' @param mzDb - A list of data.frames containing the m/z and intensity values
#' from each scan's mass spectra.
#' @param header - A data.fame containing metadata on the sample like
#' spectra type (MS1 vs MS2), retention time, and scan count.
#' @param observedPeak - A list with names 'start' and 'end' containing
#' scalar values representing the calculated peak boundary points
#' @param massThresh - A generous exact mass error threshold used to estimate
#' PPM for features.
#' @param useGap - Parameter carried into checkEICPeaks that tells Autotuner
#' whether to use the gap statustic to determine the proper number of clusters
#' to use during ppm parameter estimation.
#' @param varExpThresh - Numeric value representing the variance explained
#' threshold to use if useGap is false.
#' @param returnPpmPlots - Boolean value that tells R to return plots for
#' ppm distributions.
#' @param plotDir - Path where to store plots.
#' @param filename - A string containing the name of the current data file being
#' analyzed.
#'
#' @export
checkEICPeaks <- function(mzDb,
                          header,
                          observedPeak,
                          massThresh = 0.005,
                          useGap,
                          varExpThresh,
                          returnPpmPlots,
                          plotDir,
                          filename) {

    # extracting ms1 information for current peak -----------------------------
    ## msnObj
    suppressWarnings(rm(no_match,approvedPeaks))

    scanDiff <- diff(header$retentionTime[header$msLevel == 1L])
    rate <- mean(scanDiff)
    rm(scanDiff)

    sortedAllEIC <- dissectScans(mzDb,
                                 observedPeak = observedPeak,
                                 header = header)


    assertthat::assert_that(nrow(sortedAllEIC) > 0,
                            msg = "Check dissectScans within checkEICPeaks. A table with 0 rows was returned. Make sure correct TIC peaks were selected.")
    boundaries <- range(sortedAllEIC$scan)


    # Checking if sorted scans pass mz error threshold ------------------------
    matchedMasses <- rle(diff(sortedAllEIC$mz) < massThresh)

    ### THIS COULD BE PLACE TO ADD NOISE FILTER TO MAKE SPEED FASTER
    noiseAndPeaks <- filterPeaksfromNoise(matchedMasses)
    no_match <- noiseAndPeaks[[1]]
    truePeaks <- noiseAndPeaks[[2]]
    rm(noiseAndPeaks)

    message("-------- Number of bins detected with absolute mass error threshold: " , length(truePeaks))
    approvedPeaks <- findTruePeaks(truePeaks, sortedAllEIC)

    message("-------- Number of bins retained after checking that features within bins come from consecutive scans: ",
            nrow(approvedPeaks))
    message("-------- ", signif(nrow(approvedPeaks)/length(truePeaks)*100, digits = 2),
            " % of bins retained after checking that features come from consecutive scans")

    overlappingScans <- sum(approvedPeaks$multipleInScan)
    message("-------- Number of bins with scans with 2+ mass observations: ", overlappingScans)

    if(nrow(approvedPeaks) == 0) {
        message("No observed m/z value met was observed across adjacent scans below an error of massThresh.")
        return(NULL)
    }

    # Filtering data by variability and ppm checks ----------------------------
    ppmEst <- filterPpmError(approvedPeaks, useGap, varExpThresh,
                             returnPpmPlots, plotDir, observedPeak,
                             filename)

    assertthat::assert_that(!is.na(ppmEst),
                            msg = "Output of filterPpmError function was NA. Something may have gone wrong here with input.")

    ppmObs <- approvedPeaks$meanPPM
    ppmObs <- strsplit(split = ";", x = as.character(ppmObs)) %>%
        sapply(as.numeric)

    noisyBin <- lapply(ppmObs, function(ppm) {
        any(ppm > ppmEst)
    }) %>% unlist()
    approvScorePeaks <- approvedPeaks[!noisyBin,]

    # Estimating PeakPicking Parameters ---------------------------------------
    SNest <- estimateSNThresh(no_match,
                              sortedAllEIC, approvScorePeaks)
    SNest <- min(SNest)

    assertthat::assert_that(!is.na(SNest),
                            msg = "Output of estimateSNThresh within checkEICPeaks was NA. Something went wrong here.")

    scanEst <- min(approvScorePeaks$scanCount)

    ### Noise Intensity Estimate
    noiseEst <- min(approvScorePeaks$minIntensity) - 100
    if(noiseEst < 0) {
        noiseEst <- min(approvScorePeaks$minIntensity) + 10
    }

    ### Prefilter Intensity Estimate
    intensityEst <- min(approvScorePeaks$Intensity)/sqrt(2)

    ### peakWidth Estimate
    maxPw <- findPeakWidth(approvScorePeaks = approvScorePeaks,
                           mzDb = mzDb,
                           header = header,
                           sortedAllEIC = sortedAllEIC,
                           boundaries = boundaries,
                           ppmEst = ppmEst)

    minPw <- scanEst * rate
    if(max(maxPw) < 10*min(maxPw)) {

        message("Expanding Max Peakwidth")
        maxPw <- maxPw*2

    }


    # Returning Data ----------------------------------------------------------
    estimatedPeakData <- data.frame(ppm = ppmEst,
                                    noiseThreshold = noiseEst,
                                    peakCount = nrow(approvedPeaks),
                                    prefilterI = intensityEst,
                                    prefilterScan = scanEst,
                                    TenPercentQuanSN = unname(SNest),
                                    maxPw, minPw)

    return(estimatedPeakData)
}



