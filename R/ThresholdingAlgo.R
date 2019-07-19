#' @title ThresholdingAlgo
#'
#' @description This function performs a sliding window analysis on the
#' chromatograms in order to identify peaks within the data. I would recommend
#' to keep influence low in order to use adjacent peak lengths as a measure of
#' peak width.
#'
#' @param y - numerical vector of measured chromatographic intensity values
#' @param lag - scalar value of number of observations to calculate intensity
#' prior to peak selection.
#' @param threshold - number of standard deviations above chromatogram. Used
#' to detect significantly observed peaks.
#' @param influence - scalar values between 0-1 that describes how much the
#' value of a peak (measured index value above threshold) should contribute to
#' the sliding window analysis of downstream peaks.
#'
#'
#' @examples
#' Autotuner <- readRDS(system.file("extdata/Autotuner.rds",
#' package="Autotuner"))
#' lapply(getAutoIntensity(Autotuner), ThresholdingAlgo,
#' lag, threshold, influence)
#'
#' @return A list of calculated sliding window values.
#'
#' @export
ThresholdingAlgo <- function(y, lag, threshold, influence) {

    # correct for possible NA values in data
    assertthat::assert_that(length(y) > 0,
                        msg = paste("Error: Intensity slot within Autotuner",
                                        "Object is zero length."))

    signals <- rep(0,length(y))
    filteredY <- y[1:lag]
    avgFilter <- NULL
    stdFilter <- NULL
    avgFilter[lag] <- mean(y[1:lag])
    stdFilter[lag] <- sd(y[1:lag])

    for (i in (lag+1):length(y)){
        if (abs(y[i]-avgFilter[i-1]) > threshold*stdFilter[i-1]) {
            if (y[i] > avgFilter[i-1]) {
                signals[i] <- 1;
            } else {
                signals[i] <- -1;
            }
            filteredY[i] <- influence*y[i]+(1-influence)*filteredY[i-1]
        } else {
            signals[i] <- 0
            filteredY[i] <- y[i]
        }
        avgFilter[i] <- mean(filteredY[(i-lag):i])
        stdFilter[i] <- sd(filteredY[(i-lag):i])

    }
    return(list("signals"=signals,"avgFilter"=avgFilter,"stdFilter"=stdFilter))
}
