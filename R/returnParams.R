#' @title returnParams
#'
#' @description This function is designed to return a list of data.frames
#' containing parameter estimates obtained from the EIC and TIC parameter
#' estimation.
#'
#' @param eicParamEsts - The objection containing all parameter estimates
#' obtained from running Autotuner's EICparam function.
#' @param Autotuner - An Autotuner object used to return the TIC estimated
#' parameters
#'
#' @return A list of data.frames with parameter estimates.
#'
#' @export

returnParams <- function(eicParamEsts, Autotuner) {

    params <- TIC_params(Autotuner@peak_table, Autotuner@peak_difference)
    params <- data.frame(descriptions = names(params), estimates = unlist(params))

    colNameCheck <- all(c("ppm", "peakCount",
                          "noiseThreshold", "prefilterI",
                          "prefilterScan",
                          "TenPercentQuanSN","maxPw", "minPw") %in% colnames(eicParamEsts))

    assertthat::assert_that(colNameCheck,
                            msg = "Error in EICparams: some of the column names from the output are missing. Cannot complete parameter estimation.")

    ppmEst <- weighted.mean(eicParamEsts$ppm, eicParamEsts$peakCount)
    noiseEst <- min(eicParamEsts$noiseThreshold, na.rm = T)
    prefilterIEst <- min(eicParamEsts$prefilterI, na.rm = T)
    prefilterScanEst <- min(eicParamEsts$prefilterScan, na.rm = T)
    snEst <- min(eicParamEsts$TenPercentQuanSN, na.rm = T)

    maxPw <- split(eicParamEsts$maxPw, eicParamEsts$sampleID)
    maxPw <- sapply(maxPw, max)
    maxPw <- mean(maxPw)

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
                                      estimates = signif(estimates,
                                                         digits = 4),
                                      'Variability Measure' = signif(variability,
                                                                     digits = 4),
                                      "Measure" = description)

    output <- list(eicParams = aggregatedEstimates, ticParams = params)

    return(output)

}
