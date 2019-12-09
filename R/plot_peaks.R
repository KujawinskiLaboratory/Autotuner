#'  @title plot_peaks
#'
#'  @description This funciton plots the peak identified within
#'  chromatography.
#'
#' @details This function plots individual peaks selected by signal
#' processing and expanded with a regression to allow the user to validate the
#' selected signal processing parameters.
#'
#' @param Autotuner An Autotuner objected containing sample specific raw
#' data.
#' @param boundary UI input value that defines the boundary around the peak to
#' visualize it.
#' @param peak A Numeric index obtained from UI that indicates which peak
#' should be visualized.
#'
#' @import mzR
#'
#' @return This function outputs plots that are meant to go into the peakVis
#' UI.
#'
#' @examples
#' Autotuner <- readRDS(system.file("extdata/Autotuner.rds",
#' package="Autotuner"))
#' plot_peaks(Autotuner = Autotuner, boundary = 100, peak = 1)
#'
#' @export
plot_peaks <- function(Autotuner, boundary = 10, peak) {


    # extracting relevant values from input args ------------------------------
    factorCol <- getAutoFactorCol(Autotuner)
    metadata <- getAutoMetadata(Autotuner)
    peak_difference <- getAutoPeak_difference(Autotuner)
    peak_table <- getAutoPeak_table(Autotuner)


    peak_difference <- peak_difference
    peak_difference <- peak_difference[peak_difference$index == peak,]
    #peak_table <- peak_table
    sample_names <- unlist(metadata[,factorCol])
    sample_names <- paste(sample_names, seq_along(sample_names))

    colors <- seq_along(unique(peak_table$Sample))

    row <- 1
    peak_counter <- 1

    checkPeaks <- c(peak_difference$cur_row[1], peak_difference$next_row)
    sample_index <- peak_table$Sample[checkPeaks]

    ## ploting the peaks
    lapply(seq_along(checkPeaks), function(row_index) {

        ## renaming info for clarity
        current_row <- checkPeaks[row_index]
        sample_index <- peak_table$Sample[current_row]

        ## extracting relevant info to plot figures
        bdd_names <- peak_table[current_row,]
        bdd_names <- bdd_names[,grep("name", colnames(bdd_names))]
        time <- getAutoTime(Autotuner)[[sample_index]]
        intensity <- getAutoIntensity(Autotuner)[[sample_index]]
        bdd_points <- which(names(time) %in% unlist(bdd_names))


        peak_interval <- c(0,0)
        if((bdd_points[1] - boundary) < 1) {
            peak_interval[1] <- 1
        } else {
            peak_interval[1] <- bdd_points[1]-boundary
        }

        if((bdd_points[2]+boundary) > length(time)) {
            peak_interval[2] <- length(time)
        } else {
            peak_interval[2] <- bdd_points[2]+boundary
        }

        peak_interval <- peak_interval[1]:peak_interval[2]

        if(row_index == 1) {
            upper_bdd <- peak_difference$max_intensity[row] +
                peak_difference$max_intensity[row]/
                log(peak_difference$max_intensity[row])
            plot(x = time[peak_interval], y = intensity[peak_interval],
                 type = "l",
                 xlab = "Time (s)",
                 ylab = "Intensity",
                 main = paste("Max Peak Width:",
                              signif(max(peak_difference$Max_width)), "(s)"),
                 col = colors[peak_table$Sample[current_row]],
                 ylim = c(0, upper_bdd))
            abline(v = time[bdd_points[1]], lty = 5,
                   col = colors[peak_table$Sample[current_row]])
            abline(v = time[bdd_points[2]], lty = 5,
                   col = colors[peak_table$Sample[current_row]])
        } else {
            lines(time[peak_interval],
                  intensity[peak_interval],
                  col = colors[peak_table$Sample[current_row]])
            abline(v = time[bdd_points[1]], lty = 5,
                   col = colors[peak_table$Sample[current_row]])
            abline(v = time[bdd_points[2]], lty = 5,
                   col = colors[peak_table$Sample[current_row]])
        }
    }) # end of plotting function

    sampleIds <- peak_table$Sample[checkPeaks]

    legend("topleft",
           legend = sample_names[sampleIds],
           col = colors[sampleIds],
           cex = 0.75,
           fill = which(sample_names %in%
                            sample_names[peak_table$Sample[checkPeaks]]))


}
