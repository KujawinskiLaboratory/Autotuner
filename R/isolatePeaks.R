#' @title isolatePeaks
#'
#' @description This function is designed to handle the isolation of TIC peak
#' regions to throw into AutoTuner.
#'
#' @param Autotuner - Autotuner object. Ideally generated right after signal
#' processing peak identification is complete.
#' @param returned_peaks - numerical value representing the number of peaks
#' that are to be returned to the user for downstream parameter optimization.
#' @param signals - List of list containing data from sliding window analysis.
#'
#' @return Returns an Autotuner object with selected TIC regions.
#'
#' @export
isolatePeaks <- function(Autotuner, returned_peaks, signals) {

    Autotuner@peaks <- extract_peaks(Autotuner,
                                     returned_peaks,
                                     signals)

    Autotuner@peak_table <- peakwidth_table(Autotuner, returned_peaks)

    Autotuner@peak_difference <- peak_time_difference(Autotuner@peak_table)

    return(Autotuner)

}
