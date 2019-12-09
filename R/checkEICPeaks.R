#' @title checkEICPeaks
#'
#' @description This function is the outer most function used to check for
#' individual EIC peak specific parameters.
#'
#' @param mzDb A list of data.frames containing the m/z and intensity values
#' from each scan's mass spectra.
#' @param header A data.fame containing metadata on the sample like
#' spectra type, retention time, and scan count.
#' @param observedPeak A list with names 'start' and 'end' containing
#' scalar values representing the calculated peak boundary points
#' @param massThresh A generous exact mass error threshold used to estimate
#' PPM for features.
#' @param useGap Parameter carried into checkEICPeaks that tells Autotuner
#' whether to use the gap statustic to determine the proper number of clusters
#' to use during ppm parameter estimation.
#' @param varExpThresh Numeric value representing the variance explained
#' threshold to use if useGap is false.
#' @param returnPpmPlots Boolean value that tells R to return plots for
#' ppm distributions.
#' @param plotDir Path where to store plots.
#' @param filename A string containing the name of the current data file being
#' analyzed.
#'
#' @return This function returns a peak specific set of processign parameters.
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
                            msg = paste("Check dissectScans within",
                                        "checkEICPeaks. A table with 0 rows",
                                        "was returned. Make sure correct TIC",
                                        "peaks were selected."))
    assertthat::assert_that(!all(is.na(sortedAllEIC$scanID)),
                            msg = paste("There was a problem finding spectrum",
                                        "IDs within header file for this",
                                        "data. Error occured after function",
                                        "'dissectScans'."))

    boundaries <- range(sortedAllEIC$scanID)


    # Checking if sorted scans pass mz error threshold ------------------------
    matchedMasses <- rle(diff(sortedAllEIC$mz) < massThresh)

    ### THIS COULD BE PLACE TO ADD NOISE FILTER TO MAKE SPEED FASTER
    noiseAndPeaks <- filterPeaksfromNoise(matchedMasses)
    no_match <- noiseAndPeaks[[1]]
    truePeaks <- noiseAndPeaks[[2]]
    rm(noiseAndPeaks)

    message(paste("-------- Number of bins detected with absolute mass error",
            "threshold: "),
            length(truePeaks))
    approvedPeaks <- findTruePeaks(truePeaks, sortedAllEIC)

    message(paste("-------- Number of bins retained after checking that",
                  "features within bins come from consecutive scans: "),
            nrow(approvedPeaks))
    message("-------- ", signif(nrow(approvedPeaks)/length(truePeaks)*100,
                                digits = 2),
            paste(" % of bins retained after checking that features",
                  "come from consecutive scans"))

    overlappingScans <- sum(approvedPeaks$multipleInScan)
    message("-------- Number of bins with scans with 2+ mass observations: ",
            overlappingScans)

    if(nrow(approvedPeaks) == 0) {
        message(paste("No observed m/z value met was observed across",
                      "adjacent scans below an error of massThresh."))
        return(NULL)
    }

    # Filtering data by variability and ppm checks ----------------------------
    ppmEst <- filterPpmError(approvedPeaks, useGap, varExpThresh,
                             returnPpmPlots, plotDir, observedPeak,
                             filename)

    assertthat::assert_that(!is.na(ppmEst),
                            msg = paste("Output of filterPpmError function",
                                        "was NA. Something may have gone",
                                        "wrong here with input."))

    ppmObs <- approvedPeaks$meanPPM
    ppmObs <- strsplit(split = ";", x = as.character(ppmObs))
    ppmObs <- lapply(ppmObs, as.numeric)


    noisyBin <- lapply(ppmObs, function(ppm) {
        any(ppm > ppmEst)
    })
    noisyBin <- unlist(noisyBin)
    approvScorePeaks <- approvedPeaks[!noisyBin,]

    # Estimating PeakPicking Parameters ---------------------------------------
    SNest <- estimateSNThresh(no_match,
                              sortedAllEIC, approvScorePeaks)
    SNest <- min(SNest)

    if(is.infinite(SNest)) {
        stop("There was an issue estimating the s/n threshold")
    }

    assertthat::assert_that(!is.na(SNest),
                            msg = paste("Output of estimateSNThresh within",
                                        "checkEICPeaks was NA. Something",
                                        "went wrong here."))

    scanEst <- min(approvScorePeaks$scanCount)

    ### Noise Intensity Estimate
    noiseEst <- min(approvScorePeaks$minIntensity) -
        min(approvScorePeaks$minIntensity)*.1

    ### Prefilter Intensity Estimate
    intensityEst <- min(approvScorePeaks$Intensity) -
        min(approvScorePeaks$Intensity)*0.1

    ### peakWidth Estimate
    maxPw <- findPeakWidth(approvScorePeaks = approvScorePeaks,
                           mzDb = mzDb,
                           header = header,
                           sortedAllEIC = sortedAllEIC,
                           boundaries = boundaries,
                           ppmEst = ppmEst)


    if(length(maxPw) == 0) {

        stop("Something went wrong with findPeakWidth Function...",
             " Please report this error to github issues.")
    }

    minPw <- scanEst * rate

    ## 2019-06-19
    ## fixed this since peakwidth estimate functions should return scalars
    if(maxPw < 5*minPw) {

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



