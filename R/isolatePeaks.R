#' @title isolatePeaks
#'
#' @description This function is designed to handle the isolation of TIC peak
#' regions to throw into AutoTuner.
#'
#' @param Autotuner An Autotuner object. Ideally generated right after signal
#' processing peak identification is complete.
#' @param returned_peaks A numerical value representing the number of peaks
#' that are to be returned to the user for downstream parameter optimization.
#' @param signals A list of list containing data from sliding window analysis.
#'
#' @return Returns an Autotuner object with selected TIC regions.
#'
#'
#' @examples
#' Autotuner <- readRDS(system.file("extdata/Autotuner.rds",
#' package="Autotuner"))
#' lag <- 25
#' threshold<- 3.1
#' influence <- 0.1
#' signals <- lapply(getAutoIntensity(Autotuner), ThresholdingAlgo,
#' lag, threshold, influence)
#' isolatePeaks(Autotuner, returned_peaks = 10, signals)
#'
#' @export
isolatePeaks <- function(Autotuner, returned_peaks, signals) {

    for(i in seq_along(signals)) {
        nameCheck <- all(unique(names(signals[[i]])) %in%
                             c("signals","avgFilter","stdFilter"))
        if(!nameCheck) {
            stop("Parameter signals input is incorrect. Check slinding window.")
        }

    }

    peaks <- extract_peaks(Autotuner,returned_peaks,signals)
    Autotuner <- setAutoPeaks(peaks = peaks,Autotuner = Autotuner)

    peak_table <- peakwidth_table(Autotuner, returned_peaks)
    Autotuner <- setAutoPeak_table(peak_table = peak_table,
                                   Autotuner = Autotuner)

    peak_difference <- peak_time_difference(Autotuner)
    Autotuner <- setAutoPeak_difference(peak_difference = peak_difference,
                           Autotuner = Autotuner)

    return(Autotuner)

}
