#'  @title plot_peaks - this funciton plots the peak identified within chromatography.
#'
#' @description This function plots individual peaks selected by signal
#' processing and expanded with a regression to allow the user to validate the
#' selected signal processing parameters.
#'
#' @param Autotuner - An Autotuner objected containing sample specific raw
#' data.
#' @param boundary - UI input value that defines the boundary around the peak to
#' visualize it.
#' @param peak - Numeric index obtained from UI that indicates which peak
#' should be visualized.
#' @param peak_difference - A data.frame containing peak overlap information.
#' @param peak_table - A data.frame containing the sample specific peak
#' information.
#'
#' @importFrom rlang .data
#'
#' @return This function outputs plots that are meant to go into the peakVis
#' UI.
#'
#' @export
plot_peaks <- function(Autotuner, boundary = 10, peak, peak_difference,
                       peak_table) {


    # extracting relevant values from input args ------------------------------
    factorCol <- Autotuner@factorCol
    metadata <- Autotuner@metadata

    peak_difference <- peak_difference %>%
        dplyr::filter(.data$index == peak)
    peak_table <- peak_table
    sample_names <- unlist(metadata[,factorCol])
    sample_names <- paste(sample_names, 1:length(sample_names))


    row <- 1
    peak_counter <- 1
    while(row <= dim(peak_difference)[1]) {


        # Updating the itterator for the while loop -------------------------------
        if(row < dim(peak_difference)[1] && identical(peak_difference$cur_row[row],
                                                      peak_difference$cur_row[row+1])) {

            same_peak <- which(peak_difference$cur_row ==
                                   peak_difference$cur_row[row])
            match_rows <- c(peak_difference$cur_row[row],
                            peak_difference$next_row[same_peak])
            row_update <- row + 2

        } else {
            match_rows <- c(peak_difference[row,"cur_row"],
            peak_difference[row,"next_row"])
            row_update <- row + 1
        } # end of update to itterator


        colors <- 1:length(unique(peak_table$Sample))
        if(length(boundary) == 0) {
            boundary <- 1
        }


        # making plots ------------------------------------------------------------
        lapply(1:length(match_rows), function(row_index) {

            ## renaming info for clarity
            current_row <- match_rows[row_index]
            sample_index <- peak_table$Sample[current_row]

            ## extracting relevant info to plot figures
            bdd_names <- peak_table[current_row,] %>%
                dplyr::select(dplyr::contains("name"))
            time <- Autotuner@time[[sample_index]]
            intensity <- Autotuner@intensity[[sample_index]]
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
                peak_difference$max_intensity[row]/ log(peak_difference$max_intensity[row])
                plot(time[peak_interval],
                intensity[peak_interval],
                type = "l",
                xlab = "Time (s)",
                ylab = "Intensity",
                main = paste("Max Peak Width:", signif(max(peak_difference$Max_width)), "(s)"),
                col = colors[peak_table$Sample[current_row]],
                ylim = c(0, upper_bdd))
                abline(v = time[bdd_points[1]], lty = 5, col = colors[peak_table$Sample[current_row]])
                abline(v = time[bdd_points[2]], lty = 5, col = colors[peak_table$Sample[current_row]])
            } else {
                lines(time[peak_interval],
                intensity[peak_interval],
                col = colors[peak_table$Sample[current_row]])
                abline(v = time[bdd_points[1]], lty = 5, col = colors[peak_table$Sample[current_row]])
                abline(v = time[bdd_points[2]], lty = 5, col = colors[peak_table$Sample[current_row]])
            }
        }) # end of plotting function

        legend("topleft",
        legend = sample_names[peak_table$Sample[match_rows]],
        col = colors,
        cex = 0.75,
        fill = which(sample_names %in% sample_names[peak_table$Sample[match_rows]]))

        peak_counter <- 1 + peak_counter
        row <- row_update
    } # end of while loop

}
